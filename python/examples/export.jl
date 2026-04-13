# Export graph + simulation snapshots
#
# Run from the CO2BatchFill repo root:
# julia --project=examples python/examples/export.jl
# Can then load the exported data in Python with:
# cd python && uv run -m examples.load_and_train

using CO2BatchFill
using SurfaceWaterIntegratedModeling
using GaussianRandomFields
using Random
using NPZ, JSON  # triggers the ExportExt extension

# Setup
Random.seed!(101)
const NX, NY = 100, 100
const DX, DY = 10.0, 10.0
const LAYER_THICK = 10.0

cov = CovarianceFunction(2, Matern(200, 2, σ=3.0))
pts = range(0.0, stop=(NX - 1) * DX, step=DX)
sample_surface(base) = sample(GaussianRandomField(cov, CirculantEmbedding(), pts, pts, minpadding=113)) .+ base

surf_L1 = sample_surface(950.0)
surf_L2 = sample_surface(900.0)
surf_L3 = sample_surface(850.0)

sand_layers = [
    Dict{String,Any}("name" => "Storage layer 3", "top" => surf_L3, "base" => surf_L3 .+ LAYER_THICK),
    Dict{String,Any}("name" => "Storage layer 2", "top" => surf_L2, "base" => surf_L2 .+ LAYER_THICK),
    Dict{String,Any}("name" => "Storage layer 1", "top" => surf_L1, "base" => surf_L1 .+ LAYER_THICK),
]
topo = GenericTopography(sand_layers, NX, NY, DX, DY, minimum(surf_L3), maximum(surf_L1) + LAYER_THICK)
domain = create_domain(topo, 0.1)
layers = analyze_base_surfaces(topo; boundary_condition=:closed)

rp = ReservoirProperties(0.3, 0.4, 0.1, 25_000.0, 5.0)
rp_cap = ReservoirProperties(0.3, 0.4, 0.1, Inf, 5.0)
all_rp = [rp, rp, rp_cap]

# Run simulation
rate_L1 = create_injection_rate(layers, (div(NX, 2), div(NY, 2)), 80_000.0)
null_rate = create_no_injection(layers)

injection_events = [
    [InjectionEvent(0.0, rate_L1), InjectionEvent(10.0, create_injection_rate(layers, (1, 1), 0.0))],
    null_rate,
    null_rate,
]

seqs, ls, we = fill_layers(layers, domain, all_rp, injection_events; verbose=false)

# Generate snapshots
timepoints = collect(range(0.0, 15.0, length=50))
snapshots = generate_multi_layer_snapshots(layers, seqs, ls, we, timepoints)

# Export everything
export_dir = joinpath(@__DIR__, "..", "exported_graph_data")
export_graph_data(layers, domain, all_rp; path=export_dir)
export_snapshots(snapshots; path=export_dir)
