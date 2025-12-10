
function generate_injection_events(layers::Vector{Layer})
    # Define injection events
    # Make xy in the center of the domain
    trap_topo = layers[1].trap_structure.topography
    xy = CartesianIndex(div(size(trap_topo, 1), 2), div(size(trap_topo, 2), 2))

    # Some constant rate for testing
    annual_rates_mt = fill(0.86, 15) * 0.5

    # Convert Mt/year to m³/year
    co2_density_l1 = 425.0  # The same as used in reservoir properties
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
        injection_rate[xy] = rate
        injection_rate[xy[1] + 20, xy[2] + 20] = rate

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