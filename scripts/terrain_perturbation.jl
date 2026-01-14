"""
Terrain perturbation functions for Monte Carlo uncertainty analysis.

Provides functions for:
- Sampling points from a 2D surface
- Adding spatially correlated noise to surface elevations
- Interpolating perturbed surfaces
- Creating perturbed SleipnerTopography objects
"""

using Random
using Statistics
using Interpolations

"""
    perturb_surface(
        surface::Matrix{Float64};
        noise_std::Float64=1.0,
        sample_fraction::Float64=0.5,
        interpolation_method::Symbol=:linear
    )

Perturb a 2D topography surface by sampling points, adding noise, and interpolating.

# Arguments
- `surface::Matrix{Float64}`: Original 2D surface (depth values)
- `noise_std::Float64=1.0`: Standard deviation of Gaussian noise to add (in same units as surface)
- `sample_fraction::Float64=0.5`: Fraction of grid points to sample (0.0-1.0)
- `interpolation_method::Symbol=:linear`: Interpolation method (:linear or :cubic)

# Returns
- `Matrix{Float64}`: Perturbed surface with same dimensions as input

# Notes
- Higher sample_fraction preserves more detail but increases computation
- Lower sample_fraction leads to smoother perturbations
- Noise is added independently to each sampled point
"""
function perturb_surface(
    surface::Matrix{Float64};
    noise_std::Float64=1.0,
    sample_fraction::Float64=0.5,
    interpolation_method::Symbol=:linear
)
    nx, ny = size(surface)

    # Compute step size from sample fraction
    # sample_fraction ≈ 1/step^2, so step ≈ 1/sqrt(sample_fraction)
    step = max(1, Int(round(1.0 / sqrt(sample_fraction))))

    # Create subsampled grid
    sub_x = 1:step:nx
    sub_y = 1:step:ny

    # Generate noise on subsampled points
    sub_noise = randn(length(sub_x), length(sub_y)) .* noise_std

    # Create interpolation object
    if interpolation_method == :cubic && length(sub_x) >= 4 && length(sub_y) >= 4
        itp = interpolate(sub_noise, BSpline(Cubic(Natural(OnGrid()))))
    else
        itp = interpolate(sub_noise, BSpline(Linear()))
    end

    # Scale interpolation to original grid
    scaled_itp = Interpolations.scale(itp, sub_x, sub_y)
    extrap_itp = extrapolate(scaled_itp, Flat())

    # Apply interpolated noise to entire surface
    perturbed = similar(surface)
    for i in 1:nx
        for j in 1:ny
            perturbed[i, j] = surface[i, j] + extrap_itp(i, j)
        end
    end

    return perturbed
end


"""
    perturb_surface_regular_grid(
        surface::Matrix{Float64};
        noise_std::Float64=1.0,
        sample_spacing::Int=2,
        interpolation_method::Symbol=:linear
    )

Perturb a 2D topography surface using regular grid sampling.

# Arguments
- `surface::Matrix{Float64}`: Original 2D surface (depth values)
- `noise_std::Float64=1.0`: Standard deviation of Gaussian noise (in surface units)
- `sample_spacing::Int=2`: Grid spacing for sampling (1 = every point, 2 = every 2nd, etc.)
- `interpolation_method::Symbol=:linear`: Interpolation method (:linear or :cubic)

# Returns
- `Matrix{Float64}`: Perturbed surface with same dimensions as input

# Notes
- Lower sample_spacing (more samples) preserves more high-frequency variation
- Higher sample_spacing leads to smoother, longer-wavelength perturbations
- sample_spacing=1 means every point gets independent noise (no smoothing)
"""
function perturb_surface_regular_grid(
    surface::Matrix{Float64};
    noise_std::Float64=1.0,
    sample_spacing::Int=2,
    interpolation_method::Symbol=:linear
)
    nx, ny = size(surface)

    # Handle edge case: if spacing is 1, just add noise directly
    if sample_spacing == 1
        return surface .+ randn(nx, ny) .* noise_std
    end

    # Create subsampled grid indices
    sub_x = 1:sample_spacing:nx
    sub_y = 1:sample_spacing:ny

    # Need at least 2 points in each dimension for interpolation
    if length(sub_x) < 2 || length(sub_y) < 2
        # Fall back to direct noise addition
        return surface .+ randn(nx, ny) .* noise_std
    end

    # Generate noise on subsampled grid
    sub_noise = randn(length(sub_x), length(sub_y)) .* noise_std

    # Create interpolation object
    if interpolation_method == :cubic && length(sub_x) >= 4 && length(sub_y) >= 4
        itp = interpolate(sub_noise, BSpline(Cubic(Natural(OnGrid()))))
    else
        itp = interpolate(sub_noise, BSpline(Linear()))
    end

    # Scale interpolation to subsampled coordinates
    scaled_itp = Interpolations.scale(itp, sub_x, sub_y)
    extrap_itp = extrapolate(scaled_itp, Flat())

    # Create perturbed surface
    perturbed = similar(surface)
    for i in 1:nx
        for j in 1:ny
            perturbed[i, j] = surface[i, j] + extrap_itp(i, j)
        end
    end

    return perturbed
end


"""
    create_perturbed_topography(
        base_topography::SleipnerTopography;
        noise_std::Float64=1.0,
        sample_spacing::Int=3,
        perturb_all_surfaces::Bool=true,
        correlation_between_layers::Float64=0.0
    )

Create a new SleipnerTopography with perturbed depth surfaces.

# Arguments
- `base_topography::SleipnerTopography`: Original topography to perturb
- `noise_std::Float64=1.0`: Standard deviation of depth perturbation (meters)
- `sample_spacing::Int=3`: Grid spacing for noise sampling (controls smoothness)
- `perturb_all_surfaces::Bool=true`: If true, perturb all surfaces; if false, only sand layer tops
- `correlation_between_layers::Float64=0.0`: Correlation of noise between adjacent layers (0-1)

# Returns
- `SleipnerTopography`: New topography with perturbed surfaces

# Notes
- Each surface is perturbed independently unless correlation_between_layers > 0
- The perturbations preserve the relative ordering of surfaces (no layer inversions)
- Setting correlation_between_layers > 0 ensures similar perturbations propagate vertically
"""
function create_perturbed_topography(
    base_topography::SleipnerTopography;
    noise_std::Float64=1.0,
    sample_spacing::Int=3,
    perturb_all_surfaces::Bool=true,
    correlation_between_layers::Float64=0.0
)
    nx, ny = base_topography.nx, base_topography.ny

    # Generate base noise field that can be shared across layers
    base_noise = nothing
    if correlation_between_layers > 0
        # Create a base noise field
        sub_x = 1:sample_spacing:nx
        sub_y = 1:sample_spacing:ny
        base_noise_grid = randn(length(sub_x), length(sub_y))
        itp = interpolate(base_noise_grid, BSpline(Linear()))
        scaled_itp = Interpolations.scale(itp, sub_x, sub_y)
        extrap_itp = extrapolate(scaled_itp, Flat())

        base_noise = zeros(nx, ny)
        for i in 1:nx
            for j in 1:ny
                base_noise[i, j] = extrap_itp(i, j)
            end
        end
    end

    # Helper function to perturb a single surface with optional correlation
    function perturb_with_correlation(surface::Matrix{Float64})
        if correlation_between_layers > 0 && base_noise !== nothing
            # Mix correlated and independent noise
            independent_noise = perturb_surface_regular_grid(
                zeros(nx, ny);
                noise_std=1.0,
                sample_spacing=sample_spacing
            )
            correlated_part = base_noise .* sqrt(correlation_between_layers)
            independent_part = independent_noise .* sqrt(1 - correlation_between_layers)
            noise = (correlated_part .+ independent_part) .* noise_std
            return surface .+ noise
        else
            return perturb_surface_regular_grid(
                surface;
                noise_std=noise_std,
                sample_spacing=sample_spacing
            )
        end
    end

    # Create new surfaces dictionary
    new_surfaces = Dict{String, Any}()

    # Perturb all named surfaces if requested
    if perturb_all_surfaces
        for (name, surface) in base_topography.surfaces
            if surface isa Matrix{Float64}
                new_surfaces[name] = perturb_with_correlation(surface)
            else
                new_surfaces[name] = surface
            end
        end
    else
        # Copy surfaces without perturbation
        for (name, surface) in base_topography.surfaces
            new_surfaces[name] = surface isa Matrix{Float64} ? copy(surface) : surface
        end
    end

    # Perturb top caprock
    new_top_caprock = perturb_all_surfaces ?
        perturb_with_correlation(base_topography.top_caprock) :
        copy(base_topography.top_caprock)

    # Create new sand_layers with perturbed surfaces
    new_sand_layers = Vector{Dict{String, Any}}()

    for layer in base_topography.sand_layers
        new_layer = Dict{String, Any}()
        new_layer["id"] = layer["id"]
        new_layer["name"] = layer["name"]

        # Perturb top and base surfaces
        new_layer["top"] = perturb_with_correlation(layer["top"])
        new_layer["base"] = perturb_with_correlation(layer["base"])

        push!(new_sand_layers, new_layer)
    end

    # Compute new depth range
    all_depths = Float64[]
    push!(all_depths, minimum(new_top_caprock), maximum(new_top_caprock))
    for layer in new_sand_layers
        push!(all_depths, minimum(layer["base"]), maximum(layer["base"]))
    end

    return SleipnerTopography(
        new_surfaces,
        new_top_caprock,
        new_sand_layers,
        nx,
        ny,
        base_topography.dx,
        base_topography.dy,
        minimum(all_depths),
        maximum(all_depths)
    )
end
