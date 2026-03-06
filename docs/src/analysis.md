# Analysis and Snapshots

After running [`fill_layers`](@ref), the raw output (fill sequences and leakage
states) must be post-processed into human-readable snapshots for inspection and
plotting.

## Snapshot types

```@docs
CO2BatchFill.LayerSnapshot
CO2BatchFill.MultiLayerSnapshot
```

## Generating snapshots

```@docs
CO2BatchFill.generate_multi_layer_snapshots
CO2BatchFill.generate_multi_layer_snapshot
CO2BatchFill.generate_layer_snapshots
CO2BatchFill.generate_layer_snapshot
```

## Mass-balance helpers

These functions extract specific quantities from a [`LayerSnapshot`](@ref).

```@docs
CO2BatchFill.total_to_next_layer
CO2BatchFill.total_upward_leakage
CO2BatchFill.total_passthrough
```

## Summary printing

```@docs
CO2BatchFill.print_summary
```
