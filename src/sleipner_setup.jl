export generate_reservoir_properties_for_sleipner_layers, generate_sleipner_injection_events


function generate_reservoir_properties_for_sleipner_layers()::Vector{ReservoirProperties}

    n_layers = 9

    # From L1 up to L9. Using values from paper.
    brine_density = 1020
    co2_density = [570.0, 542.5, 515.0, 487.5, 460.0, 432.5, 405.0, 377.5, 350.0]
    brine_co2_density_difference = brine_density .- co2_density

    reservoir_properties = Vector{ReservoirProperties}(undef, n_layers)
    for i in 1:n_layers
        reservoir_properties[i] = ReservoirProperties(
            0.4,
            0.2,
            0.3,
            98000.0,
            brine_co2_density_difference[i],
            1.0
        )
    end
    return reservoir_properties
end


"""
    generate_sleipner_injection_events(layers, injection_cell::CartesianIndex)

Generate injection events for Sleipner simulation based on historical injection rates (1996-2010).
Injection occurs only in the bottom layer (L1) at the specified cell location.

# Arguments
- `layers`: Vector of layer structures
- `injection_cell`: CartesianIndex specifying the injection location in the bottom layer

# Returns
- `Vector{Vector{InjectionEvent}}`: Injection events for each layer (only L1 has non-zero injection)

# Notes
- Annual injection rates from Sleipner 2019 Benchmark (total 12.18 Mt over 14 years)
- Rates converted from Mt/year to m³/year using CO₂ density of 570 kg/m³
- Injection location is approximately at grid center (I=32, J=59) in the model
"""
function generate_sleipner_injection_events(
    layers,
    injection_cell::CartesianIndex
)::Vector{Vector{InjectionEvent}}

    # Historical annual injection rates from 1996-2010 (Mt/year)
    # Source: Sleipner 2019 Benchmark model
    annual_rates_mt = [
        0.07,  # 1996
        0.67,  # 1997
        0.85,  # 1998
        0.94,  # 1999
        0.94,  # 2000
        1.02,  # 2001
        0.96,  # 2002
        0.92,  # 2003
        0.76,  # 2004
        0.87,  # 2005
        0.83,  # 2006
        0.93,  # 2007
        0.82,  # 2008
        0.86,  # 2009
        0.76   # 2010
    ]

    # Convert Mt/year to m³/year
    # 1 Mt = 1e9 kg
    # CO₂ density at L1 (bottom layer) = 570 kg/m³
    co2_density_l1 = 570.0  # kg/m³
    annual_rates_m3_per_year = annual_rates_mt .* 1e9 ./ co2_density_l1

    # Create injection events for bottom layer (L1)
    # Time points are cumulative: start at year 0, events mark end of each year
    n_events = length(annual_rates_mt)
    bottom_layer_events = Vector{InjectionEvent}(undef, n_events)

    # Get the grid size from the bottom layer
    grid_size = size(layers[1].trap_structure.topography)

    for (i, rate) in enumerate(annual_rates_m3_per_year)
        # Time in years (0, 1, 2, ..., 14)
        time = float(i - 1)

        # Create injection rate field (only inject at specified cell)
        injection_rate = zeros(grid_size)
        injection_rate[injection_cell] = rate

        bottom_layer_events[i] = InjectionEvent(time, injection_rate)
    end

    # Create zero injection events for all other layers
    zero_injection = zeros(grid_size)
    zero_event = [InjectionEvent(0.0, zero_injection)]

    # Assemble injection events for all layers
    n_layers = length(layers)
    injection_events = Vector{Vector{InjectionEvent}}(undef, n_layers)
    injection_events[1] = bottom_layer_events  # Bottom layer (L1) has actual injection
    for i in 2:n_layers
        injection_events[i] = zero_event  # All other layers have zero injection
    end

    return injection_events
end