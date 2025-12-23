import os
import sys
import re
import glob
import imageio.v2 as imageio
import matplotlib.pyplot as plt

def safe_filename(name: str) -> str:
    name = name.strip()
    name = re.sub(r"[\x00-\x1f\x7f]", "", name)      # control chars
    name = re.sub(r'[<>:"/\\|?*]', "_", name)        # illegal chars
    return name

def infer_dt_id(*paths) -> str:
    text = " ".join(paths).replace("\\", "/")
    m = re.search(r"\bdt0*([0-9]+)\b", text, flags=re.IGNORECASE)
    return f"dt{m.group(1)}" if m else "dt"

def read_points(path):
    pts = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cid, x, y = line.split()
            pts[int(cid)] = (float(x), float(y))
    return pts

def read_path(path):
    order = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.lower().startswith("distance"):
                continue
            order.append(int(line))
    return order

def draw(points, order, out_png, title=""):
    xs = [points[i][0] for i in order] + [points[order[0]][0]]
    ys = [points[i][1] for i in order] + [points[order[0]][1]]

    plt.figure(figsize=(6, 6))
    plt.plot(xs, ys, "-o")
    plt.scatter(xs[:-1], ys[:-1])

    for i in order:
        plt.text(points[i][0], points[i][1], str(i), fontsize=9)

    plt.title(title)
    plt.axis("equal")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

def main():
    # Usage:
    #   python plot_bonus_gif.py dataset/dt4/point.txt output_bonus/dt4_frames
    # Optional:
    #   python plot_bonus_gif.py point.txt frames_dir output.gif
    if len(sys.argv) < 3:
        print("Usage: python plot_bonus_gif.py <point.txt> <frames_dir> [output_gif]")
        return

    point_file = sys.argv[1]
    frames_dir = sys.argv[2]

    if not os.path.exists(point_file):
        print(f"[ERROR] point file not found: {point_file}")
        return
    if not os.path.isdir(frames_dir):
        print(f"[ERROR] frames dir not found: {frames_dir}")
        return

    dt = infer_dt_id(point_file, frames_dir)
    frames_dir_name = os.path.basename(os.path.normpath(frames_dir))

    # default gif name (unique per dataset & frames folder)
    default_gif = safe_filename(f"{frames_dir_name}.gif")
    out_gif = sys.argv[3] if len(sys.argv) >= 4 else default_gif
    out_gif = safe_filename(out_gif)

    # store intermediate pngs in a dedicated folder
    png_dir = safe_filename(f"gif_frames_{dt}")
    os.makedirs(png_dir, exist_ok=True)

    points = read_points(point_file)
    path_files = sorted(glob.glob(os.path.join(frames_dir, "path_*.txt")))

    if not path_files:
        print(f"[ERROR] No path_*.txt found in: {frames_dir}")
        return

    images = []
    for p in path_files:
        step = os.path.basename(p).replace(".txt", "")  # e.g. path_000100
        png_path = os.path.join(png_dir, safe_filename(f"{step}.png"))

        order = read_path(p)
        draw(points, order, png_path, title=f"{dt}  {step}")

        images.append(imageio.imread(png_path))

    imageio.mimsave(out_gif, images, duration=0.5)
    print(f"[OK] GIF saved: {os.path.abspath(out_gif)}")
    print(f"[OK] PNG frames in: {os.path.abspath(png_dir)}")

if __name__ == "__main__":
    main()
