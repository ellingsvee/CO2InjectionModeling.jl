# =============================================================================
# Synthetic Multi-Layer CO2 Storage Example
# =============================================================================
#
# This example demonstrates the core workflow of CO2BatchFill using a small,
# fully synthetic reservoir. No external data files are needed.
#
# Scenario: 3 sand layers with structural traps created by Gaussian domes.
# CO2 is injected into the deepest layer and can migrate upward through
# caprock leakage points when CO2 column height exceeds a threshold.
#
# Workflow:
#   1. Build synthetic topography (depth surfaces with dome features)
#   2. Analyze trap structures in each layer
#   3. Define reservoir properties and injection schedule
#   4. Run multi-layer simulation
#   5. Generate snapshots and inspect results

using CO2BatchFill
using SurfaceWaterIntegratedModeling

# =============================================================================
# 1. Create synthetic topography
# =============================================================================
#
# We define a small 30x30 grid with 50 m cell spacing (1500 m x 1500 m domain).
# Each sand layer's top surface has a gentle regional dip plus Gaussian domes
# that create structural highs — these act as CO2 traps.
#
# In the depth convention used here, smaller values = shallower = structural high.
# SWIM's spillanalysis sees these as topographic depressions where "water" (CO2)
# accumulates.

nx, ny = 30, 30
dx, dy = 50.0, 50.0  # meters

# Coordinate arrays (cell centers)
xs = [(i - 0.5) * dx for i in 1:nx]
ys = [(j - 0.5) * dy for j in 1:ny]

"""
    make_surface(base_depth, slope_x, slope_y, domes, xs, ys)

Create a depth surface with regional dip and Gaussian dome features.

Each dome is specified as (center_x, center_y, sigma, amplitude).
Domes are subtracted from the background to create structural highs (depth minima).
"""
function make_surface(base_depth, slope_x, slope_y, domes, xs, ys)
    nx, ny = length(xs), length(ys)
    surface = zeros(nx, ny)
    for i in 1:nx, j in 1:ny
        surface[i, j] = base_depth + slope_x * xs[i] + slope_y * ys[j]
        for (cx, cy, sigma, amp) in domes
            r2 = (xs[i] - cx)^2 + (ys[j] - cy)^2
            surface[i, j] -= amp * exp(-r2 / (2 * sigma^2))
        end
    end
    return surface
end

# Layer 1 (deepest, ~900 m) — two domes
surface_1 = make_surface(900.0, 0.002, 0.001,
    [(500.0, 400.0, 200.0, 4.0),     # Dome A
     (1000.0, 800.0, 200.0, 3.0)],   # Dome B
    xs, ys)

# Layer 2 (middle, ~850 m) — two domes
surface_2 = make_surface(850.0, 0.002, 0.001,
    [(600.0, 600.0, 200.0, 3.5),     # Dome C
     (900.0, 400.0, 150.0, 2.5)],    # Dome D
    xs, ys)

# Layer 3 (shallowest, ~800 m) — two domes
surface_3 = make_surface(800.0, 0.002, 0.001,
    [(400.0, 700.0, 250.0, 4.0),     # Dome E
     (1100.0, 500.0, 180.0, 3.0)],   # Dome F
    xs, ys)

# Sand layer thickness (top to base of each sand body)
thickness = 10.0  # meters

# Build layer dictionaries
# NOTE: sand_layers must be ordered shallowest-first, deepest-last.
# analyze_base_surfaces() reverses this to produce layers ordered deepest-first,
# which is the convention used by fill_layers().
sand_layers = [
    Dict{String,Any}("name" => "Sand 3", "top" => surface_3, "base" => surface_3 .+ thickness),
    Dict{String,Any}("name" => "Sand 2", "top" => surface_2, "base" => surface_2 .+ thickness),
    Dict{String,Any}("name" => "Sand 1", "top" => surface_1, "base" => surface_1 .+ thickness),
]

# Create topography object
topography = GenericTopography(
    sand_layers, nx, ny, dx, dy,
    minimum(surface_3),                # shallowest depth
    maximum(surface_1 .+ thickness);   # deepest depth
)

# =============================================================================
# 2. Analyze layers and create domain
# =============================================================================
#
# analyze_base_surfaces runs SWIM's spillanalysis on each layer's top surface
# to identify structural traps, their volumes, and spill relationships.
#
# With boundary_condition=:closed, walls are added around the domain edges
# so CO2 cannot escape laterally.

boundary_condition = :closed
domain = create_domain(topography, 1.0)
layers = analyze_base_surfaces(topography; boundary_condition=boundary_condition)

println("Domain: $(domain.nx) x $(domain.ny), dx=$(domain.dx) m, dy=$(domain.dy) m")
println("Number of layers: $(length(layers))")
for (i, layer) in enumerate(layers)
    nt = SurfaceWaterIntegratedModeling.numtraps(layer.trap_structure)
    println("  Layer $i ($(layer.name)): $nt traps")
end

# =============================================================================
# 3. Define reservoir properties
# =============================================================================
#
# ReservoirProperties defines petrophysical parameters and caprock strength.
# The leakage_height is automatically computed from the pressure threshold:
#   h = P / ((rho_brine - rho_co2) * g)
#
# For the top layer, we use Inf pressure threshold (sealed caprock, no leakage
# out of the reservoir).

rp_leaky = ReservoirProperties(
    0.30,       # sand_porosity
    0.20,       # sand_residual_co2_saturation
    0.10,       # sand_irreducible_water_saturation
    20000.0,    # shale_pressure_threshold (Pa) — ~3.6 m leakage height
    5.0;        # residual_leakage_time (years)
)

rp_sealed = ReservoirProperties(
    0.30, 0.20, 0.10,
    Inf,        # shale_pressure_threshold — impermeable caprock
    5.0;
)

# Per-layer properties: bottom two layers can leak, top layer is sealed
reservoir_properties = [rp_leaky, rp_leaky, rp_sealed]

println("\nReservoir properties:")
println("  Porosity: $(rp_leaky.sand_porosity)")
println("  Leakage height (layers 1-2): $(round(rp_leaky.leakage_height, digits=2)) m")
println("  Leakage height (layer 3): sealed (Inf)")

# =============================================================================
# 4. Define injection events
# =============================================================================
#
# CO2 is injected into the deepest layer (layers[1] = "Sand 1") for 10 years
# at a single cell near Dome A. Other layers receive CO2 only through leakage
# from the layer below.
#
# The injection rate matrix must match the padded topography size (the padding
# is added by analyze_base_surfaces for closed boundary conditions).

injection_events = Vector{Vector{InjectionEvent}}()

for (i, layer) in enumerate(layers)
    topo_size = size(layer.trap_structure.topography)
    pad = layer.boundary_padding

    if i == 1  # Deepest layer — inject here
        injection_rate = zeros(topo_size)
        # Inject near Dome A center, adjusting for boundary padding
        inj_i = round(Int, 500.0 / dx) + pad
        inj_j = round(Int, 400.0 / dy) + pad
        injection_rate[inj_i, inj_j] = 50000.0  # m³/year

        push!(injection_events, [
            InjectionEvent(0.0, injection_rate),          # Start injection
            InjectionEvent(10.0, zeros(topo_size)),       # Stop at t=10 years
        ])
    else  # Upper layers — no direct injection
        push!(injection_events, [InjectionEvent(0.0, zeros(topo_size))])
    end
end

println("\nInjection: 50,000 m³/year into layer 1 for 10 years (500,000 m³ total)")

# =============================================================================
# 5. Run simulation
# =============================================================================

println("\nRunning multi-layer simulation...")
seqs, leakage_states = fill_layers(
    layers, domain, reservoir_properties, injection_events;
    verbose=true
)
println("Simulation complete.\n")

# Print leakage summary
for (i, ls) in enumerate(leakage_states)
    n_leaking = count(ls.leaking)
    if n_leaking > 0
        println("Layer $i: $n_leaking trap(s) leaking through caprock")
    else
        println("Layer $i: no leakage")
    end
end

# =============================================================================
# 6. Generate and inspect snapshots
# =============================================================================

snapshots = generate_reservoir_snapshots(
    layers, seqs, leakage_states, domain,
    reservoir_properties, injection_events;
    num_snapshots=10,
    start_time=0.0,
    end_time=15.0,
    verbose=true
)

# Print summary at final time
print_snapshot_summary(snapshots[end])

# Print per-layer details
for ls in snapshots[end].layer_snapshots
    print_layer_snapshot_summary(ls)
end
