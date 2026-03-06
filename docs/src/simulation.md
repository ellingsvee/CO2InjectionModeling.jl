# Simulation

## Multi-layer simulation

[`fill_layers`](@ref) runs the multi-layer simulation.  It fills each layer in
depth order (deepest first), converts leakage from one layer into injection
events for the layer above, and returns the complete fill sequences and leakage
states.

```@docs
CO2BatchFill.fill_layers
```

## Single-layer simulation

For single-layer use cases, [`fill_sequence_with_leakage`](@ref)
simulates one layer directly.

```@docs
CO2BatchFill.fill_sequence_with_leakage
```
