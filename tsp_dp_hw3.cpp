#include <bits/stdc++.h>
using namespace std;

struct City {
    int id;
    double x, y;
};

static inline double dist2d(const City& a, const City& b) {
    double dx = a.x - b.x;
    double dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}

// Read point file: each line -> id x y
static bool read_points(const string& path, vector<City>& cities) {
    ifstream fin(path);
    if (!fin) return false;
    cities.clear();
    int id; double x, y;
    while (fin >> id >> x >> y) {
        cities.push_back({id, x, y});
    }
    return !cities.empty();
}

// Held–Karp DP solve; start is index 0
static bool solve_tsp_dp(const vector<City>& cities,
                         double& bestCost,
                         vector<int>& orderIndex /* indices in cities */) {
    int n = (int)cities.size();
    if (n == 0) return false;
    if (n == 1) {
        bestCost = 0.0;
        orderIndex = {0};
        return true;
    }

    // Precompute distances
    vector<vector<double>> d(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            d[i][j] = dist2d(cities[i], cities[j]);

    const double INF = 1e100;
    int FULL = 1 << n;
    // dp[mask][i] = min cost from 0 to i visiting mask; mask represents visited set; i represents last visited
    vector<vector<double>> dp(FULL, vector<double>(n, INF));
    vector<vector<int>> parent(FULL, vector<int>(n, -1));

    // Base case : start at 0
    dp[1 << 0][0] = 0.0;

    // (mask, i) -> unvisited j, update new mask. If shorter, update dp and record parent.
    for (int mask = 0; mask < FULL; mask++) {
        if (!(mask & (1 << 0))) continue;
        for (int i = 0; i < n; i++) {
            if (!(mask & (1 << i))) continue;
            double cur = dp[mask][i];
            if (cur >= INF / 2) continue;
            for (int j = 0; j < n; j++) {
                if (mask & (1 << j)) continue;
                int nmask = mask | (1 << j);
                double cand = cur + d[i][j];
                if (cand < dp[nmask][j]) {
                    dp[nmask][j] = cand;
                    parent[nmask][j] = i;
                }
            }
        }
    }

    int fullMask = FULL - 1;
    bestCost = INF;
    int bestLast = -1;
    for (int i = 1; i < n; i++) {
        double cand = dp[fullMask][i] + d[i][0]; // dp[fullMask][i] shorteset path to i + return to start
        if (cand < bestCost) {
            bestCost = cand;
            bestLast = i; // find min best last node
        }
    }
    if (bestLast == -1) return false;

    // Reconstruct path indices: 0 -> ... -> bestLast
    // backtrack from bestLast
    vector<int> rev;
    int mask = fullMask;
    int cur = bestLast;
    rev.push_back(cur);
    while (cur != 0) {
        int p = parent[mask][cur];
        if (p == -1) return false;
        mask ^= (1 << cur);
        cur = p;
        rev.push_back(cur);
    }
    reverse(rev.begin(), rev.end()); // starts at 0
    orderIndex = rev;
    return true;
}

// Write output
static bool write_answer(const string& path,
                         double bestCost,
                         const vector<City>& cities,
                         const vector<int>& orderIndex) {
    ofstream fout(path);
    if (!fout) return false;
    fout << fixed << setprecision(3);
    fout << "distance: " << bestCost << "\n";
    for (int idx : orderIndex) fout << cities[idx].id << "\n";
    return true;
}

static string join_path(const string& a, const string& b) {
#ifdef _WIN32
    const char sep = '\\';
#else
    const char sep = '/';
#endif
    if (a.empty()) return b;
    if (a.back() == '/' || a.back() == '\\') return a + b;
    return a + sep + b;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Modes:
    // 1) Single-file
    // 2) Batch

    if (argc >= 2 && string(argv[1]) == "--batch") {
        string datasetRoot, outDir;

        if (argc >= 3) datasetRoot = argv[2];
        else {
            cout << "Enter dataset root dir (e.g., dataset): ";
            cin >> datasetRoot;
        }

        if (argc >= 4) outDir = argv[3];
        else {
            // default output dir name
            outDir = "output_hw3";
        }

        // create output directory (portable-ish)
        // For simplicity: rely on system mkdir
#ifdef _WIN32
        string cmd = "mkdir " + outDir + " >nul 2>nul";
#else
        string cmd = "mkdir -p " + outDir + " >/dev/null 2>&1";
#endif
        system(cmd.c_str());

        vector<string> dts = {"dt1", "dt2", "dt3"}; // dt4 is skipped for HW3
        for (int k = 0; k < (int)dts.size(); k++) {
            string dt = dts[k];
            string pointPath = join_path(join_path(datasetRoot, dt), "point.txt");
            string outPath   = join_path(outDir, "ans_" + dt + ".txt");

            vector<City> cities;
            if (!read_points(pointPath, cities)) {
                cerr << "[HW3][SKIP] Cannot read " << pointPath << "\n";
                continue;
            }

            // DP warning if too large
            if ((int)cities.size() > 22) {
                cerr << "[HW3][WARN] " << dt << ": n=" << cities.size()
                     << " may be too large for DP.\n";
            }

            double bestCost;
            vector<int> orderIdx;
            bool ok = solve_tsp_dp(cities, bestCost, orderIdx);
            if (!ok) {
                cerr << "[HW3][FAIL] " << dt << ": DP failed.\n";
                continue;
            }

            if (!write_answer(outPath, bestCost, cities, orderIdx)) {
                cerr << "[HW3][FAIL] Cannot write " << outPath << "\n";
                continue;
            }

            cout << "[HW3][OK] " << dt
                 << "  n=" << cities.size()
                 << "  best=" << fixed << setprecision(3) << bestCost
                 << "  -> " << outPath << "\n";
        }

        cout << "[HW3] Batch done. (dt4 is skipped for HW3)\n";
        return 0;
    }

    // ---- Single-file mode ----
    string inputFile, outputFile;
    if (argc >= 2) inputFile = argv[1];
    else {
        cout << "Enter input point file path (e.g., dataset/dt1/point.txt): ";
        cin >> inputFile;
    }
    if (argc >= 3) outputFile = argv[2];
    else {
        cout << "Enter output file path (e.g., ans_dt01.txt): ";
        cin >> outputFile;
    }

    vector<City> cities;
    if (!read_points(inputFile, cities)) {
        cerr << "ERROR: Cannot read input file: " << inputFile << "\n";
        return 1;
    }

    if ((int)cities.size() > 22) {
        cerr << "WARNING: n=" << cities.size()
             << " may be too large for DP (Held–Karp).\n";
    }

    double bestCost;
    vector<int> orderIdx;
    if (!solve_tsp_dp(cities, bestCost, orderIdx)) {
        cerr << "ERROR: DP failed.\n";
        return 1;
    }

    if (!write_answer(outputFile, bestCost, cities, orderIdx)) {
        cerr << "ERROR: Cannot write output file: " << outputFile << "\n";
        return 1;
    }

    cout << "Solved with DP. Best distance: " << fixed << setprecision(3) << bestCost << "\n";
    cout << "Output written to: " << outputFile << "\n";
    return 0;
}
