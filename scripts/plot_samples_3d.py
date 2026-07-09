from pathlib import Path

import matplotlib.pyplot as plt


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


def main():
    script_dir = Path(__file__).resolve().parent
    csv_path = script_dir / "samples.csv"

    if not csv_path.exists():
        raise FileNotFoundError(
            f"Missing {csv_path}. Run the C program first to generate samples."
        )

    xs, ys, zs = load_samples(csv_path)

    fig = plt.figure(figsize=(9, 7))
    ax = fig.add_subplot(111, projection="3d")

    scatter = ax.scatter(xs, ys, zs, c=zs, cmap="viridis", s=20, alpha=0.9)
    ax.set_title("Magnetometer Samples (3D)")
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")

    fig.colorbar(scatter, ax=ax, pad=0.1, label="Z value")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
