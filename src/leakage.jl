
export compute_leakage_height

"""
Compute the leakage height based on the reservoir properties.
Uses the relation: cappilary_pressure_threshold = density_difference_between_brine_and_co2 * g * column_height
"""
function compute_leakage_height(reservoir_properties::ReservoirProperties) 
    g = 9.81 # m/s^2
    height = reservoir_properties.shale_pressure_threshold / (reservoir_properties.brine_co2_density_difference * g)
    # height = 15.0 # TEMPORARY FIX FOR TESTING. Know leakage should occur here.
    return height
end