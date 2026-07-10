from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np


def load_samples(csv_path: Path):
    xs = []
    ys = []
    zs = []

    with csv_path.open("r", encoding="utf-8") as f:
        next(f)  # Skip header: i,x,y,z
        for line in f:
            parts = line.strip().split(",")
            if len(parts) != 4:
                continue
            _, x, y, z = parts
            xs.append(float(x))
            ys.append(float(y))
            zs.append(float(z))

    return xs, ys, zs


def load_ellipsoid_params(path: Path):
    offset = None
    transform_rows = []

    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split(",")
            key = parts[0]
            values = [float(v) for v in parts[1:]]

            if key == "offset" and len(values) == 3:
                offset = np.array(values, dtype=float)
            elif key.startswith("transform_row") and len(values) == 3:
                transform_rows.append(values)

    if offset is None or len(transform_rows) != 3:
        return None, None

    transform = np.array(transform_rows, dtype=float)
    return offset, transform


def create_ellipsoid_surface(offset, transform, resolution=40):
    # Calibration does: v_calibrated = T * (v_raw - offset)
    # For plotting in raw space, map unit sphere points by inverse(T):
    # v_raw = offset + inv(T) * v_unit
    inv_transform = np.linalg.inv(transform)

    u = np.linspace(0.0, 2.0 * np.pi, resolution)
    v = np.linspace(0.0, np.pi, resolution)

    x = np.outer(np.cos(u), np.sin(v))
    y = np.outer(np.sin(u), np.sin(v))
    z = np.outer(np.ones_like(u), np.cos(v))

    sphere_points = np.stack([x, y, z], axis=-1)
    raw_points = sphere_points @ inv_transform.T + offset

    ex = raw_points[:, :, 0]
    ey = raw_points[:, :, 1]
    ez = raw_points[:, :, 2]
    return ex, ey, ez


def set_axes_equal(ax):
    """
    Wyrównuje skale osi X, Y, Z na wykresie 3D.
    Dzięki temu kula będzie wyglądać jak kula, a nie elipsoida.
    """
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = x_limits[1] - x_limits[0]
    x_middle = sum(x_limits) / 2
    
    y_range = y_limits[1] - y_limits[0]
    y_middle = sum(y_limits) / 2
    
    z_range = z_limits[1] - z_limits[0]
    z_middle = sum(z_limits) / 2

    # Znajdź maksymalny zakres (promień) ze wszystkich osi
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    # Ustaw nowe, symetryczne limity
    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    
    # Wymuś, aby ramka wykresu była sześcianem (wymaga matplotlib >= 3.3.0)
    ax.set_box_aspect([1, 1, 1])


def main():
    script_dir = Path(__file__).resolve().parent
    csv_path = script_dir / "samples.csv"
    calib_path = script_dir / "calibrated.csv"
    ellipsoid_path = script_dir / "ellipsoid.csv"

    if not csv_path.exists():
        raise FileNotFoundError(
            f"Missing {csv_path}. Run the C program first to generate samples."
        )

    xs, ys, zs = load_samples(csv_path)
    has_calib = calib_path.exists()
    has_ellipsoid = ellipsoid_path.exists()

    if has_calib:
        fig = plt.figure(figsize=(16, 7))
        ax1 = fig.add_subplot(121, projection="3d")
        ax2 = fig.add_subplot(122, projection="3d")
    else:
        fig = plt.figure(figsize=(9, 7))
        ax1 = fig.add_subplot(111, projection="3d")
        ax2 = None

    # Wykres 1: Surowe dane
    scatter1 = ax1.scatter(xs, ys, zs, c=zs, cmap="viridis", s=20, alpha=0.9)
    ax1.set_title("Magnetometer Samples (Uncalibrated)")
    ax1.set_xlabel("X")
    ax1.set_ylabel("Y")
    ax1.set_zlabel("Z")
    fig.colorbar(scatter1, ax=ax1, pad=0.1, label="Z value")

    if has_ellipsoid:
        offset, transform = load_ellipsoid_params(ellipsoid_path)
        if offset is not None and transform is not None:
            ex, ey, ez = create_ellipsoid_surface(offset, transform)
            ax1.plot_surface(ex, ey, ez, color="orange", alpha=0.22, linewidth=0)
    
    # Wyrównanie osi dla pierwszego wykresu
    set_axes_equal(ax1)

    # Wykres 2: Skalibrowane dane
    if has_calib and ax2 is not None:
        c_xs, c_ys, c_zs = load_samples(calib_path)
        
        scatter2 = ax2.scatter(c_xs, c_ys, c_zs, c=c_zs, cmap="plasma", s=20, alpha=0.9)
        ax2.set_title("Magnetometer Samples (Calibrated)")
        ax2.set_xlabel("X")
        ax2.set_ylabel("Y")
        ax2.set_zlabel("Z")
        fig.colorbar(scatter2, ax=ax2, pad=0.1, label="Z value")
        
        # Wyrównanie osi dla drugiego wykresu
        set_axes_equal(ax2)

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()