"""Load exported graph data and do a simple plot."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from co2graph import build_multi_layer_graph, load_graph_data
from co2graph.plot import plot_co2_heights

PACKAGE_PATH = Path(__file__).parent.parent
DATA_DIR = PACKAGE_PATH / "exported_graph_data"


def main():
    # Load exported data
    data = load_graph_data(DATA_DIR)
    print(f"Loaded {data.n_layers} layers, grid {data.nx}x{data.ny}")
    for layer in data.layers:
        print(f"- Layer {layer.index} ({layer.name}): {layer.n_traps} traps")

    # Build multi-layer graph (all layers + leakage edges)
    mg = build_multi_layer_graph(data, include_leakage_edges=True)
    print("\nMulti-layer graph:")
    print(f"- nodes: {mg['n_nodes']}, edges: {mg['n_edges']}")
    print(f"- node_features: {mg['node_features'].shape}")
    print(f"- edge_features: {mg['edge_features'].shape}")
    print(f"- layer_offsets: {mg['layer_offsets']}")

    # Plot trap volumes on the spatial grid
    # Use dummy volumes (half capacity) to demonstrate plotting
    layer = data.layers[0]
    dummy_volumes = np.array(
        [node.net_capacity * 0.5 for node in layer.nodes], dtype=np.float64
    )
    plot_co2_heights(data, layer_idx=1, volumes=dummy_volumes, max_co2_height=5.0)
    plt.savefig(PACKAGE_PATH / "example_co2_heights.svg", bbox_inches="tight")
    print(f"\nSaved example plot to {PACKAGE_PATH / 'example_co2_heights.svg'}")


if __name__ == "__main__":
    main()
