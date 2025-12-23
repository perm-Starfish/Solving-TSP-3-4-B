import os
import sys
import re
import matplotlib.pyplot as plt

def safe_filename(name: str) -> str:
    name = name.strip()
    name = re.sub(r"[\x00-\x1f\x7f]", "", name)      # control chars
    name = re.sub(r'[<>:"/\\|?*]', "_", name)        # illegal chars
    return name

def read_points(point_path):
    pts = {}
    with open(point_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            cid, x, y = line.split()
            pts[int(cid)] = (float(x), float(y))
    return pts

def read_ans_hw4(ans_path):
    order = []
    with open(ans_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            low = line.lower()
            if low.startswith("mean") or low.startswith("distance"):
                continue
            order.append(int(line))
    return order

def draw(points, order, out_png, title):
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
    # python plot_hw4_batch.py dataset output_hw4 fig_hw4
    dataset_root = sys.argv[1] if len(sys.argv) >= 2 else "dataset"
    ans_dir      = sys.argv[2] if len(sys.argv) >= 3 else "output_hw4"
    fig_dir      = sys.argv[3] if len(sys.argv) >= 4 else "fig_hw4"

    os.makedirs(fig_dir, exist_ok=True)

    # HW4: dt1-dt4
    for dt in ["dt1", "dt2", "dt3", "dt4"]:
        point_path = os.path.join(dataset_root, dt, "point.txt")
        ans_path   = os.path.join(ans_dir, f"ans_{dt}.txt")

        if not os.path.exists(point_path):
            print(f"[SKIP] missing point: {point_path}")
            continue
        if not os.path.exists(ans_path):
            print(f"[SKIP] missing ans:   {ans_path}")
            continue

        points = read_points(point_path)
        order  = read_ans_hw4(ans_path)

        out_png = os.path.join(fig_dir, safe_filename(f"fig_{dt}_hw4.png"))
        draw(points, order, out_png, title=f"HW4 - ACO TSP ({dt})")

        print(f"[OK] {dt} -> {os.path.abspath(out_png)}")

if __name__ == "__main__":
    main()
