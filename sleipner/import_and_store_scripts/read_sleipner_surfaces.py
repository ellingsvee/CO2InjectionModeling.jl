"""
Script to read Sleipner reservoir depth surfaces from RMS grid format files.
Extracts geometry of all layers and stores them as numpy arrays.
"""

import numpy as np
from pathlib import Path
from typing import Dict, Tuple


def parse_rms_grid(filepath: Path) -> Tuple[np.ndarray, Dict]:
    """
    Parse an RMS grid format file.

    Parameters
    ----------
    filepath : Path
        Path to the RMS grid file

    Returns
    -------
    grid_data : np.ndarray
        2D array of depth values with shape (ny, nx)
    metadata : dict
        Dictionary containing grid metadata (dimensions, bounding box, etc.)
    """
    with open(filepath, "r") as f:
        lines = f.readlines()

    # Parse header information
    metadata = {}
    data_start_idx = None

    for i, line in enumerate(lines):
        # Look for grid dimensions and bounding box (line with nx, ny, xmin, xmax, ymin, ymax)
        if (
            "," in line
            and not line.strip().startswith("!")
            and not line.strip().startswith("@")
        ):
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 2:
                # Check if this line contains the grid dimensions
                nums = []
                for part in parts:
                    try:
                        nums.extend([float(x) for x in part.split()])
                    except ValueError:
                        continue

                if (
                    len(nums) >= 6 and i < 15
                ):  # Grid dimension line typically in first 15 lines
                    metadata["nx"] = int(nums[0])
                    metadata["ny"] = int(nums[1])
                    metadata["xmin"] = nums[2]
                    metadata["xmax"] = nums[3]
                    metadata["ymin"] = nums[4]
                    metadata["ymax"] = nums[5]

        # Find where actual data starts (after the '@' and '+ Grid data starts' lines)
        if "+ Grid data starts after this line" in line or "+Grid data starts" in line:
            data_start_idx = i + 1
            break

    if data_start_idx is None:
        raise ValueError(f"Could not find data start marker in {filepath}")

    # Read all numeric data after the header
    data_values = []
    for line in lines[data_start_idx:]:
        line = line.strip()
        if line and not line.startswith("!") and not line.startswith("@"):
            try:
                values = [float(x) for x in line.split()]
                data_values.extend(values)
            except ValueError:
                continue

    # Reshape into 2D grid (nx columns, ny rows)
    nx = metadata["nx"]
    ny = metadata["ny"]

    # The data is stored row-by-row (row-major order)
    grid_data = np.array(data_values[: nx * ny]).reshape(ny, nx)

    # Replace undefined values with NaN
    undefined_value = -99999.0
    grid_data[grid_data == undefined_value] = np.nan

    # Extract surface name from filepath
    metadata["surface_name"] = filepath.name

    return grid_data, metadata


def read_all_sleipner_surfaces(data_dir: Path) -> Dict[str, Tuple[np.ndarray, Dict]]:
    """
    Read all Sleipner depth surface files from the DepthSurfaces_Grid directory.

    Parameters
    ----------
    data_dir : Path
        Path to the DepthSurfaces_Grid directory

    Returns
    -------
    surfaces : dict
        Dictionary mapping surface name to (grid_data, metadata) tuple
    """
    surfaces = {}

    # List of expected surface files based on the documentation
    surface_files = [
        # Top of intrashales (Reflectors 1-7)
        "Reflector1",
        "Reflector2",
        "Reflector3",
        "Reflector4",
        "Reflector5",
        "Reflector6",
        "Reflector7",
        # Base of intrashales (Base reflectors 1-7)
        "Base_Reflector1",
        "Base_Reflector2",
        "Base_Reflector3",
        "Base_Reflector4",
        "Base_Reflector5",
        "Base_Reflector6",
        "Base_Reflector7",
        # Top and base of Utsira Formation
        "TopUtsiraFm",
        "BaseUtsiraFm",
        # Top of Sand Wedge
        "TopSW",
        # Top of caprock
        "Top_Caprock",
        # Additional surface
        "ThickShale",
    ]

    for surface_name in surface_files:
        filepath = data_dir / surface_name
        if filepath.exists():
            try:
                grid_data, metadata = parse_rms_grid(filepath)
                surfaces[surface_name] = (grid_data, metadata)
                print(
                    f"✓ Loaded {surface_name}: shape {grid_data.shape}, "
                    f"depth range [{np.nanmin(grid_data):.2f}, {np.nanmax(grid_data):.2f}] m"
                )
            except Exception as e:
                print(f"✗ Error loading {surface_name}: {e}")
        else:
            print(f"✗ File not found: {surface_name}")

    return surfaces


def save_surfaces_to_npz(
    surfaces: Dict[str, Tuple[np.ndarray, Dict]], output_file: Path
):
    """
    Save all surfaces to a single compressed numpy file.

    Parameters
    ----------
    surfaces : dict
        Dictionary of surfaces from read_all_sleipner_surfaces()
    output_file : Path
        Output .npz file path
    """
    # Prepare data for saving
    arrays_dict = {}
    metadata_dict = {}

    for name, (grid_data, metadata) in surfaces.items():
        arrays_dict[name] = grid_data
        # Store metadata as string (npz doesn't handle nested dicts well)
        metadata_dict[name] = str(metadata)

    # Save all arrays and metadata
    np.savez_compressed(output_file, **arrays_dict, _metadata=str(metadata_dict))
    print(f"\nSaved all surfaces to {output_file}")


def load_surfaces_from_npz(input_file: Path) -> Dict[str, np.ndarray]:
    """
    Load surfaces from a compressed numpy file.

    Parameters
    ----------
    input_file : Path
        Input .npz file path

    Returns
    -------
    surfaces : dict
        Dictionary mapping surface name to grid data array
    """
    data = np.load(input_file)
    surfaces = {key: data[key] for key in data.files if key != "_metadata"}
    return surfaces


def load_sleipner_9_layers(
    input_file: Path,
) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Load the 9 permeable sand layers of the Sleipner Utsira Formation for CO2 storage modeling.

    The Sleipner reservoir consists of 9 PERMEABLE SAND layers where CO2 is stored and flows,
    separated by 7 thin INTRASHALE BARRIERS (shale layers ~1-2m thick) that act as seals.

    CO2 is stored WITHIN the sand layers (between the shale barriers), NOT in the barriers.
    The layers are numbered from bottom (Layer 1) to top (Layer 9).

    Parameters
    ----------
    input_file : Path
        Path to the .npz file containing all surfaces

    Returns
    -------
    layers : dict
        Dictionary with keys 'Layer_1' through 'Layer_9', each containing:
        - 'top': 2D numpy array of top surface depths (m)
        - 'base': 2D numpy array of base surface depths (m)
        - 'thickness': 2D numpy array of layer thickness (m)
        - 'description': String describing the permeable sand layer
    """
    surfaces = load_surfaces_from_npz(input_file)

    layers = {}

    # Layer 9: Topmost sand layer (PERMEABLE - CO2 is injected here)
    layers["Layer_9"] = {
        "top": surfaces["TopSW"],
        "base": surfaces["ThickShale"],
        "thickness": surfaces["ThickShale"] - surfaces["TopSW"],
        "description": "Topmost permeable sand (injection layer)",
    }

    # Layer 9: Topmost sand layer (PERMEABLE - CO2 is injected here)
    layers["Layer_8"] = {
        "top": surfaces["TopUtsiraFm"],
        "base": surfaces["Reflector7"],
        "thickness": surfaces["Reflector7"] - surfaces["TopUtsiraFm"],
        "description": "Topmost permeable sand (injection layer)",
    }

    # Layer 8: Sand between intrashale barriers 7 and 6 (PERMEABLE)
    layers["Layer_7"] = {
        "top": surfaces["Base_Reflector7"],
        "base": surfaces["Reflector6"],
        "thickness": surfaces["Reflector6"] - surfaces["Base_Reflector7"],
        "description": "Permeable sand between barriers 7 and 6",
    }

    # Layer 7: Sand between intrashale barriers 6 and 5 (PERMEABLE)
    layers["Layer_6"] = {
        "top": surfaces["Base_Reflector6"],
        "base": surfaces["Reflector5"],
        "thickness": surfaces["Reflector5"] - surfaces["Base_Reflector6"],
        "description": "Permeable sand between barriers 6 and 5",
    }

    # Layer 6: Sand between intrashale barriers 5 and 4 (PERMEABLE)
    layers["Layer_5"] = {
        "top": surfaces["Base_Reflector5"],
        "base": surfaces["Reflector4"],
        "thickness": surfaces["Reflector4"] - surfaces["Base_Reflector5"],
        "description": "Permeable sand between barriers 5 and 4",
    }

    # Layer 5: Sand between intrashale barriers 4 and 3 (PERMEABLE)
    layers["Layer_4"] = {
        "top": surfaces["Base_Reflector4"],
        "base": surfaces["Reflector3"],
        "thickness": surfaces["Reflector3"] - surfaces["Base_Reflector4"],
        "description": "Permeable sand between barriers 4 and 3",
    }

    # Layer 4: Sand between intrashale barriers 3 and 2 (PERMEABLE)
    layers["Layer_3"] = {
        "top": surfaces["Base_Reflector3"],
        "base": surfaces["Reflector2"],
        "thickness": surfaces["Reflector2"] - surfaces["Base_Reflector3"],
        "description": "Permeable sand between barriers 3 and 2",
    }

    # Layer 3: Sand between intrashale barriers 2 and 1 (PERMEABLE)
    layers["Layer_2"] = {
        "top": surfaces["Base_Reflector2"],
        "base": surfaces["Reflector1"],
        "thickness": surfaces["Reflector1"] - surfaces["Base_Reflector2"],
        "description": "Permeable sand between barriers 2 and 1",
    }

    # Layer 1: Bottommost sand layer (PERMEABLE)
    layers["Layer_1"] = {
        "top": surfaces["Base_Reflector1"],
        "base": surfaces["BaseUtsiraFm"],
        "thickness": surfaces["BaseUtsiraFm"] - surfaces["Base_Reflector1"],
        "description": "Bottommost permeable sand layer",
    }

    return layers


def load_sleipner_topographies(
    input_file: Path,
) -> np.ndarray:
    layers = load_sleipner_9_layers(input_file)
    topographies = np.zeros(
        (9, layers["Layer_1"]["top"].shape[0], layers["Layer_1"]["top"].shape[1])
    )

    # Topography is the top surface of each layer
    for i in range(9):
        layer_key = f"Layer_{i + 1}"
        topographies[i, :, :] = layers[layer_key]["top"]

    return topographies


def save_topographies_to_npz(topographies: np.ndarray, output_file: Path):
    """
    Save topographies to a compressed numpy file.

    Parameters
    ----------
    topographies : np.ndarray
        3D array of shape (9, ny, nx) containing top surface depths for each layer
    output_file : Path
        Output .npz file path
    """
    np.savez_compressed(output_file, topographies=topographies)
    print(f"\nSaved topographies to {output_file}")


def load_sleipner_topographies_from_npz(input_file: Path) -> np.ndarray:
    """
    Load topographies from a compressed numpy file.

    Parameters
    ----------
    input_file : Path
        Input .npz file path

    Returns
    -------
    topographies : np.ndarray
        3D array of shape (9, ny, nx) containing top surface depths for each layer
    """
    data = np.load(input_file)
    topographies = data["topographies"]
    return topographies


def main():
    """Main execution function."""
    # Define paths
    sleipner_dir = Path(__file__).parent.parent
    print(sleipner_dir)
    data_dir = (
        sleipner_dir / "velocities_trends_surfaces" / "data" / "DepthSurfaces_Grid"
    )
    output_dir = sleipner_dir
    output_dir.mkdir(exist_ok=True)

    print("=" * 70)
    print("Reading Sleipner Reservoir Depth Surfaces")
    print("=" * 70)
    print(f"Data directory: {data_dir}")
    print()

    # Read all surfaces
    surfaces = read_all_sleipner_surfaces(data_dir)

    print()
    print("=" * 70)
    print(f"Successfully loaded {len(surfaces)} surfaces")
    print("=" * 70)

    # Save to compressed numpy format
    output_file = output_dir / "sleipner_depth_surfaces.npz"
    save_surfaces_to_npz(surfaces, output_file)

    # Print summary statistics
    print("\nSummary of loaded surfaces:")
    print("-" * 70)

    # Group surfaces by type
    print("\nIntrashale layers (top surfaces):")
    for i in range(1, 8):
        name = f"Reflector{i}"
        if name in surfaces:
            grid, meta = surfaces[name]
            print(
                f"  {name:20s}: {grid.shape}, "
                f"depth = [{np.nanmin(grid):7.2f}, {np.nanmax(grid):7.2f}] m"
            )

    print("\nIntrashale layers (base surfaces):")
    for i in range(1, 8):
        name = f"Base_Reflector{i}"
        if name in surfaces:
            grid, meta = surfaces[name]
            print(
                f"  {name:20s}: {grid.shape}, "
                f"depth = [{np.nanmin(grid):7.2f}, {np.nanmax(grid):7.2f}] m"
            )

    print("\nMain formation boundaries:")
    for name in ["Top_Caprock", "TopSW", "TopUtsiraFm", "BaseUtsiraFm", "ThickShale"]:
        if name in surfaces:
            grid, meta = surfaces[name]
            print(
                f"  {name:20s}: {grid.shape}, "
                f"depth = [{np.nanmin(grid):7.2f}, {np.nanmax(grid):7.2f}] m"
            )

    # Example: demonstrate how to load the data back
    print("\n" + "=" * 70)
    print("Testing data reload...")
    print("=" * 70)
    loaded_surfaces = load_surfaces_from_npz(output_file)
    print(f"Loaded {len(loaded_surfaces)} surfaces from {output_file}")
    print(f"Available surfaces: {list(loaded_surfaces.keys())}")

    # Load topographies
    topographies = load_sleipner_topographies(output_file)
    print(f"Loaded topographies with shape: {topographies.shape}")
    # Save topographies to a separate file
    topography_output_file = output_dir / "sleipner_topographies.npz"
    save_topographies_to_npz(topographies, topography_output_file)


if __name__ == "__main__":
    main()
