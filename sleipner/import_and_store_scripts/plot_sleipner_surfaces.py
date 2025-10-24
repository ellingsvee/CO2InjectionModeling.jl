"""
Visualization script for Sleipner reservoir depth surfaces.
Creates multiple plots showing the geometry of the reservoir layers.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


DATA_FILE = Path(__file__).parent.parent / "sleipner_depth_surfaces.npz"
OUTPUT_PATH = Path(__file__).parent.parent / "plots"


def load_surfaces(data_file):
    """Load the surface data from npz file."""
    return np.load(data_file)


def plot_surface_map(ax, surface_data, title, cmap="terrain"):
    """Plot a single surface as a 2D map."""
    im = ax.imshow(surface_data, cmap=cmap, aspect="auto", origin="lower")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.set_xlabel("X grid index")
    ax.set_ylabel("Y grid index")
    plt.colorbar(im, ax=ax, label="Depth (m)")
    return im


def plot_all_surfaces_overview():
    """Create an overview plot showing all main surfaces."""
    # Load data
    data_file = DATA_FILE
    surfaces = load_surfaces(data_file)

    # Create figure with subplots
    fig, axes = plt.subplots(3, 3, figsize=(16, 14))
    fig.suptitle(
        "Sleipner Reservoir Depth Surfaces Overview", fontsize=16, fontweight="bold"
    )

    # Main surfaces to plot
    surface_list = [
        ("Top_Caprock", "Top of Caprock"),
        ("TopSW", "Top of Sand Wedge"),
        ("TopUtsiraFm", "Top of Utsira Formation"),
        ("Reflector7", "Reflector 7 (Top)"),
        ("Reflector5", "Reflector 5 (Top)"),
        ("Reflector3", "Reflector 3 (Top)"),
        ("Reflector1", "Reflector 1 (Top)"),
        ("BaseUtsiraFm", "Base of Utsira Formation"),
        ("ThickShale", "Thick Shale"),
    ]

    for idx, (key, title) in enumerate(surface_list):
        ax = axes.flat[idx]
        if key in surfaces.files:
            plot_surface_map(ax, surfaces[key], title)
        else:
            ax.set_visible(False)

    plt.tight_layout()
    output_file = OUTPUT_PATH / "sleipner_surfaces_overview.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Saved overview plot to: {output_file}")
    plt.close()


def plot_layer_thicknesses():
    """Plot the thicknesses of various layers."""
    data_file = DATA_FILE
    surfaces = load_surfaces(data_file)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle("Sleipner Reservoir Layer Thicknesses", fontsize=16, fontweight="bold")

    # 1. Utsira Formation thickness
    ax = axes[0, 0]
    utsira_thickness = surfaces["BaseUtsiraFm"] - surfaces["TopUtsiraFm"]
    im = ax.imshow(utsira_thickness, cmap="viridis", aspect="auto", origin="lower")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Utsira Formation Thickness", fontweight="bold")
    ax.set_xlabel("X grid index")
    ax.set_ylabel("Y grid index")
    plt.colorbar(im, ax=ax, label="Thickness (m)")

    # 2. Average intrashale thickness
    ax = axes[0, 1]
    avg_intrashale = np.zeros_like(surfaces["Reflector1"])
    for i in range(1, 8):
        thickness = surfaces[f"Base_Reflector{i}"] - surfaces[f"Reflector{i}"]
        avg_intrashale += thickness
    avg_intrashale /= 7
    im = ax.imshow(avg_intrashale, cmap="plasma", aspect="auto", origin="lower")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Average Intrashale Thickness", fontweight="bold")
    ax.set_xlabel("X grid index")
    ax.set_ylabel("Y grid index")
    plt.colorbar(im, ax=ax, label="Thickness (m)")

    # 3. Caprock to Top Utsira spacing
    ax = axes[1, 0]
    spacing = surfaces["TopUtsiraFm"] - surfaces["Top_Caprock"]
    im = ax.imshow(spacing, cmap="coolwarm", aspect="auto", origin="lower")
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Caprock to Top Utsira Spacing", fontweight="bold")
    ax.set_xlabel("X grid index")
    ax.set_ylabel("Y grid index")
    plt.colorbar(im, ax=ax, label="Spacing (m)")

    # 4. Histogram of Utsira thickness
    ax = axes[1, 1]
    valid_thickness = utsira_thickness[~np.isnan(utsira_thickness)]
    ax.hist(valid_thickness, bins=50, color="steelblue", edgecolor="black", alpha=0.7)
    ax.set_xlabel("Thickness (m)", fontweight="bold")
    ax.set_ylabel("Frequency", fontweight="bold")
    ax.set_title("Distribution of Utsira Formation Thickness", fontweight="bold")
    ax.axvline(
        np.mean(valid_thickness),
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"Mean: {np.mean(valid_thickness):.1f} m",
    )
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    output_file = OUTPUT_PATH / "sleipner_layer_thicknesses.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Saved thickness plot to: {output_file}")
    plt.close()


def plot_cross_section():
    """Create cross-sectional views of the reservoir."""
    data_file = DATA_FILE
    surfaces = load_surfaces(data_file)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    fig.suptitle("Sleipner Reservoir Cross-Sections", fontsize=16, fontweight="bold")

    # Define the layers in stratigraphic order (top to bottom)
    layer_order = [
        ("Top_Caprock", "Top Caprock", "brown"),
        ("TopSW", "Top Sand Wedge", "orange"),
        ("TopUtsiraFm", "Top Utsira", "blue"),
        ("Reflector7", "Reflector 7", "cyan"),
        ("Reflector6", "Reflector 6", "lightblue"),
        ("Reflector5", "Reflector 5", "skyblue"),
        ("Reflector4", "Reflector 4", "lightcyan"),
        ("Reflector3", "Reflector 3", "paleturquoise"),
        ("Reflector2", "Reflector 2", "powderblue"),
        ("Reflector1", "Reflector 1", "lightsteelblue"),
        ("BaseUtsiraFm", "Base Utsira", "darkblue"),
    ]

    # Cross-section 1: Middle row (X direction)
    j_slice = 32  # Middle of the Y dimension
    ax1.set_title(f"West-East Cross-Section (Y index = {j_slice})", fontweight="bold")
    x_coords = np.arange(surfaces["TopUtsiraFm"].shape[1])

    for key, label, color in layer_order:
        if key in surfaces.files:
            depth = surfaces[key][j_slice, :]
            ax1.plot(x_coords, depth, label=label, linewidth=2, color=color)

    ax1.set_xlabel("X grid index", fontweight="bold")
    ax1.set_ylabel("Depth (m)", fontweight="bold")
    ax1.legend(loc="best", fontsize=8, ncol=2)
    ax1.grid(alpha=0.3)
    ax1.invert_yaxis()  # Depth increases downward

    # Cross-section 2: Middle column (Y direction)
    i_slice = 59  # Middle of the X dimension
    ax2.set_title(f"South-North Cross-Section (X index = {i_slice})", fontweight="bold")
    y_coords = np.arange(surfaces["TopUtsiraFm"].shape[0])

    for key, label, color in layer_order:
        if key in surfaces.files:
            depth = surfaces[key][:, i_slice]
            ax2.plot(y_coords, depth, label=label, linewidth=2, color=color)

    ax2.set_xlabel("Y grid index", fontweight="bold")
    ax2.set_ylabel("Depth (m)", fontweight="bold")
    ax2.legend(loc="best", fontsize=8, ncol=2)
    ax2.grid(alpha=0.3)
    ax2.invert_yaxis()  # Depth increases downward

    plt.tight_layout()
    output_file = OUTPUT_PATH / "sleipner_cross_sections.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Saved cross-section plot to: {output_file}")
    plt.close()


def plot_3d_surface():
    """Create a 3D visualization of selected surfaces."""
    data_file = DATA_FILE
    surfaces = load_surfaces(data_file)

    fig = plt.figure(figsize=(16, 12))

    # Create 3D plot for Top Caprock
    ax1 = fig.add_subplot(2, 2, 1, projection="3d")
    surface = surfaces["Top_Caprock"]
    ny, nx = surface.shape
    X, Y = np.meshgrid(np.arange(nx), np.arange(ny))
    surf = ax1.plot_surface(X, Y, surface, cmap="terrain", alpha=0.9, edgecolor="none")
    ax1.set_title("Top Caprock Surface (3D)", fontweight="bold")
    ax1.set_xlabel("X")
    ax1.set_ylabel("Y")
    ax1.set_zlabel("Depth (m)")
    ax1.invert_zaxis()
    fig.colorbar(surf, ax=ax1, shrink=0.5)

    # Create 3D plot for Top Utsira
    ax2 = fig.add_subplot(2, 2, 2, projection="3d")
    surface = surfaces["TopUtsiraFm"]
    surf = ax2.plot_surface(X, Y, surface, cmap="viridis", alpha=0.9, edgecolor="none")
    ax2.set_title("Top Utsira Formation (3D)", fontweight="bold")
    ax2.set_xlabel("X")
    ax2.set_ylabel("Y")
    ax2.set_zlabel("Depth (m)")
    ax2.invert_zaxis()
    fig.colorbar(surf, ax=ax2, shrink=0.5)

    # Create 3D plot for Base Utsira
    ax3 = fig.add_subplot(2, 2, 3, projection="3d")
    surface = surfaces["BaseUtsiraFm"]
    surf = ax3.plot_surface(X, Y, surface, cmap="plasma", alpha=0.9, edgecolor="none")
    ax3.set_title("Base Utsira Formation (3D)", fontweight="bold")
    ax3.set_xlabel("X")
    ax3.set_ylabel("Y")
    ax3.set_zlabel("Depth (m)")
    ax3.invert_zaxis()
    fig.colorbar(surf, ax=ax3, shrink=0.5)

    # Create 3D plot showing Utsira thickness
    ax4 = fig.add_subplot(2, 2, 4, projection="3d")
    thickness = surfaces["BaseUtsiraFm"] - surfaces["TopUtsiraFm"]
    surf = ax4.plot_surface(
        X, Y, thickness, cmap="coolwarm", alpha=0.9, edgecolor="none"
    )
    ax4.set_title("Utsira Formation Thickness (3D)", fontweight="bold")
    ax4.set_xlabel("X")
    ax4.set_ylabel("Y")
    ax4.set_zlabel("Thickness (m)")
    fig.colorbar(surf, ax=ax4, shrink=0.5)

    plt.tight_layout()
    output_file = OUTPUT_PATH / "sleipner_3d_surfaces.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"Saved 3D surface plot to: {output_file}")
    plt.close()


def main():
    """Generate all visualizations."""
    print("=" * 70)
    print("Generating Sleipner Reservoir Visualizations")
    print("=" * 70)
    print()

    print("1. Creating overview of all surfaces...")
    plot_all_surfaces_overview()

    print("2. Creating layer thickness plots...")
    plot_layer_thicknesses()

    print("3. Creating cross-section plots...")
    plot_cross_section()

    print("4. Creating 3D surface visualizations...")
    plot_3d_surface()

    print()
    print("=" * 70)
    print("All visualizations completed successfully!")
    print("Check the 'sleipner' directory for output files.")
    print("=" * 70)


if __name__ == "__main__":
    main()
