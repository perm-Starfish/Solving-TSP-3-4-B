#include <bits/stdc++.h>
using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct City {
    int id;
    double x, y;
};

static inline double dist2d(double ax, double ay, double bx, double by) {
    double dx = ax - bx, dy = ay - by;
    return sqrt(dx*dx + dy*dy);
}

static inline double dist2d_city(const City& a, const City& b) {
    return dist2d(a.x, a.y, b.x, b.y);
}

static bool read_points(const string& path, vector<City>& cities) {
    ifstream fin(path);
    if (!fin) return false;
    cities.clear();
    int id; double x, y;
    while (fin >> id >> x >> y) cities.push_back({id, x, y});
    return !cities.empty();
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

static double tour_length_by_city_idx(const vector<int>& orderCityIdx, const vector<City>& cities) {
    int n = (int)orderCityIdx.size();
    if (n <= 1) return 0.0;
    double L = 0.0;
    for (int i = 0; i < n - 1; i++) {
        const City& a = cities[orderCityIdx[i]];
        const City& b = cities[orderCityIdx[i+1]];
        L += dist2d_city(a, b);
    }
    // close
    L += dist2d_city(cities[orderCityIdx.back()], cities[orderCityIdx[0]]);
    return L;
}

// rotate tour to start at city with id==1 if exists
static void rotate_start_city1(vector<int>& orderCityIdx, const vector<City>& cities) {
    int n = (int)orderCityIdx.size();
    int pos = -1;
    for (int i = 0; i < n; i++) {
        if (cities[orderCityIdx[i]].id == 1) { pos = i; break; }
    }
    if (pos <= 0) return;
    rotate(orderCityIdx.begin(), orderCityIdx.begin() + pos, orderCityIdx.end());
}

// Convert elastic net (ring nodes W) into a TSP city order.
// Rule: assign each city to its nearest ring node; then sort cities by assigned node index.
// If multiple cities share the same node index, break ties by distance to that node.
static vector<int> net_to_tour_city_order(
    const vector<pair<double,double>>& W, // size M
    const vector<City>& cities
) {
    int M = (int)W.size();
    int n = (int)cities.size();

    vector<int> assign(n, 0);
    vector<double> bestd(n, 1e100);

    for (int c = 0; c < n; c++) {
        double cx = cities[c].x, cy = cities[c].y;
        int best_i = 0;
        double bd = 1e100;
        for (int i = 0; i < M; i++) {
            double d = dist2d(cx, cy, W[i].first, W[i].second);
            if (d < bd) { bd = d; best_i = i; }
        }
        assign[c] = best_i;
        bestd[c] = bd;
    }

    vector<int> idx(n);
    iota(idx.begin(), idx.end(), 0);

    stable_sort(idx.begin(), idx.end(), [&](int a, int b){
        if (assign[a] != assign[b]) return assign[a] < assign[b];
        return bestd[a] < bestd[b];
    });

    // ensure unique cities (already unique); return as order of city indices
    return idx;
}

struct ENParams {
    int run_times = 30;                  // required
    long long eval_max = -1;             // default 10000*n
    int record_every = 100;              // required: record best path every 100 eval

    // Elastic Net hyperparameters
    int M_mul = 8;                       // M = M_mul * n (ring nodes count)
    int M_cap = 2000;                    // cap for M
    double lr0 = 0.2;                    // initial learning rate
    double lr_end = 0.02;                // final learning rate

    double sigma0 = 3.0;                 // initial neighborhood width (in ring index distance)
    double sigma_end = 0.5;              // final sigma

    double lambda0 = 0.02;               // initial elastic term weight
    double lambda_end = 0.2;             // final elastic term weight

    unsigned seed = 0;                   // 0 => random_device
};

// circular distance on ring indices
static inline int ring_dist(int i, int j, int M) {
    int d = abs(i - j);
    return min(d, M - d);
}

// initialize ring nodes in a circle around centroid, radius slightly bigger than max distance
static vector<pair<double,double>> init_ring(const vector<City>& cities, int M, mt19937& rng) {
    int n = (int)cities.size();
    double cx = 0.0, cy = 0.0;
    for (auto& c : cities) { cx += c.x; cy += c.y; }
    cx /= n; cy /= n;

    double r = 0.0;
    for (auto& c : cities) {
        r = max(r, dist2d(cx, cy, c.x, c.y));
    }
    r *= 1.2;
    if (r < 1e-6) r = 1.0;

    uniform_real_distribution<double> jit(-0.01, 0.01);

    vector<pair<double,double>> W(M);
    for (int i = 0; i < M; i++) {
        double ang = 2.0 * M_PI * (double)i / (double)M;
        W[i].first  = cx + r * cos(ang) + jit(rng);
        W[i].second = cy + r * sin(ang) + jit(rng);
    }
    return W;
}

// one evaluation step (paper-aligned):
// soft assignment eta_ij = exp(-||X - Wj||^2 / (2K^2)) / sum_k exp(...)
// then gradient step:
//   Wj += lr * [ alpha * eta_j * (X - Wj) + lambda * (W_{j-1}+W_{j+1}-2Wj) ]
static void en_step(
    vector<pair<double,double>>& W,
    const City& c,
    double lr,
    double sigma,   // NOTE: reinterpret as K (paper's scale parameter)
    double lambda
) {
    const double alpha = 1.0; // paper's attraction weight (keep fixed; lr scales overall step)

    int M = (int)W.size();
    if (M == 0) return;

    // K^2 with floor for numerical safety
    double K2 = max(sigma * sigma, 1e-12);

    // Copy for unbiased (simultaneous) update
    vector<pair<double,double>> W0 = W;

    // Compute unnormalized weights w_j = exp(-||X-Wj||^2 / (2K^2))
    // then normalize to get eta_j.
    vector<double> w(M);

    // Use max-trick for numerical stability
    double maxLog = -1e300;
    for (int j = 0; j < M; j++) {
        double dx = c.x - W0[j].first;
        double dy = c.y - W0[j].second;
        double d2 = dx*dx + dy*dy;                // squared distance
        double logwj = -d2 / (2.0 * K2);
        if (logwj > maxLog) maxLog = logwj;
        w[j] = logwj; // temporarily store log-weight
    }

    double denom = 0.0;
    for (int j = 0; j < M; j++) {
        double ej = exp(w[j] - maxLog);
        w[j] = ej;          // now store exp(logw - maxLog)
        denom += ej;
    }
    if (denom <= 0.0) return;

    // Update all nodes:
    // attraction term uses eta_j = w[j]/denom
    // elastic term uses discrete Laplacian on ring
    for (int j = 0; j < M; j++) {
        int jp = (j + 1) % M;
        int jm = (j - 1 + M) % M;

        double eta = w[j] / denom;

        // attraction (paper): alpha * eta * (X - Wj)
        double ax = alpha * eta * (c.x - W0[j].first);
        double ay = alpha * eta * (c.y - W0[j].second);

        // elastic smoothing (paper): lambda * (W_{j+1}+W_{j-1}-2Wj)
        double ex = (W0[jm].first + W0[jp].first - 2.0*W0[j].first);
        double ey = (W0[jm].second + W0[jp].second - 2.0*W0[j].second);

        W[j].first  = W0[j].first  + lr * (ax + lambda * ex);
        W[j].second = W0[j].second + lr * (ay + lambda * ey);
    }
}


static bool write_summary_ans(
    const string& outPath,
    double mean_best,
    double best_overall,
    const vector<City>& cities,
    const vector<int>& bestOrderCityIdx
) {
    ofstream fout(outPath);
    if (!fout) return false;
    fout << fixed << setprecision(3);
    fout << "mean distance: " << mean_best << "\n";
    fout << "distance: " << best_overall << "\n";
    for (int ci : bestOrderCityIdx) fout << cities[ci].id << "\n";
    return true;
}

// Write one snapshot record for GIF plotting.
// output the city visiting order (IDs) for best-so-far at that evaluation.
static bool write_snapshot_path(
    const string& path,
    const vector<City>& cities,
    const vector<int>& orderCityIdx,
    double bestDist
) {
    ofstream fout(path);
    if (!fout) return false;
    fout << fixed << setprecision(3);
    fout << "distance: " << bestDist << "\n";
    for (int ci : orderCityIdx) fout << cities[ci].id << "\n";
    return true;
}

static bool solve_elastic_net_bonus(
    const vector<City>& cities,
    const ENParams& P,
    double& mean_best_30,
    double& best_overall,
    vector<int>& best_order_city_idx,
    const string& snapshotDirForBestRun // if non-empty, write snapshots of the best run overall
) {
    int n = (int)cities.size();
    if (n < 2) return false;

    long long eval_max = P.eval_max;
    if (eval_max <= 0) eval_max = 10000LL * n;

    mt19937 rng;
    if (P.seed == 0) rng.seed(random_device{}());
    else rng.seed(P.seed);

    int M = min(P.M_mul * n, P.M_cap);
    M = max(M, n); // at least n

    vector<double> best_each_run;
    best_each_run.reserve(P.run_times);

    best_overall = 1e100;
    best_order_city_idx.clear();

    // Keep snapshots only for the run that finally becomes the best_overall for this dataset
    // store its snapshots in a temp list first, and write when confirmed.
    vector<pair<long long, vector<int>>> bestRunSnapshots; // (eval, order)
    vector<pair<long long, double>> bestRunSnapshotDist;   // (eval, bestDist)

    for (int run = 0; run < P.run_times; run++) {
        vector<pair<double,double>> W = init_ring(cities, M, rng);

        // choose random city order sampling for SGD steps
        uniform_int_distribution<int> pickCity(0, n - 1);

        double best_run = 1e100;
        vector<int> best_run_order;

        // store snapshots for this run if it becomes best overall
        vector<pair<long long, vector<int>>> runSnaps;
        vector<pair<long long, double>> runSnapDist;

        for (long long ev = 1; ev <= eval_max; ev++) {
            // linear schedules (simple + stable)
            double t = (double)(ev - 1) / (double)max(1LL, (eval_max - 1));
            double lr     = P.lr0     + (P.lr_end     - P.lr0)     * t;
            double sigma  = P.sigma0  + (P.sigma_end  - P.sigma0)  * t;
            double lambda = P.lambda0 + (P.lambda_end - P.lambda0) * t;

            const City& c = cities[pickCity(rng)];
            en_step(W, c, lr, sigma, lambda);

            // record every 100 eval: best-so-far path
            if (ev % P.record_every == 0) {

                // 1. 用目前 Elastic Net 形狀轉成一條 TSP 路徑
                vector<int> current_order = net_to_tour_city_order(W, cities);
                rotate_start_city1(current_order, cities);
                double current_L = tour_length_by_city_idx(current_order, cities);

                // 2. snapshot：存「當下路徑」（GIF 用）
                runSnaps.push_back({ev, current_order});
                runSnapDist.push_back({ev, current_L});

                // 3. best：只用來記錄最佳成績
                if (current_L < best_run) {
                    best_run = current_L;
                    best_run_order = current_order;
                }
            }
        }

        best_each_run.push_back(best_run);

        // update overall best and snapshot set
        if (best_run < best_overall) {
            best_overall = best_run;
            best_order_city_idx = best_run_order;
            bestRunSnapshots = move(runSnaps);
            bestRunSnapshotDist = move(runSnapDist);
        }
    }

    // mean best
    double sum = 0.0;
    for (double v : best_each_run) sum += v;
    mean_best_30 = sum / (double)best_each_run.size();

    // write snapshots for best run overall (for GIF)
    if (!snapshotDirForBestRun.empty()) {
#ifdef _WIN32
        string cmd = "mkdir " + snapshotDirForBestRun + " >nul 2>nul";
#else
        string cmd = "mkdir -p " + snapshotDirForBestRun + " >/dev/null 2>&1";
#endif
        system(cmd.c_str());

        // output files: path_000100.txt, path_000200.txt, ...
        for (size_t i = 0; i < bestRunSnapshots.size(); i++) {
            long long ev = bestRunSnapshots[i].first;
            const auto& order = bestRunSnapshots[i].second;
            double dist = bestRunSnapshotDist[i].second;

            ostringstream oss;
            oss << "path_" << setw(6) << setfill('0') << ev << ".txt";
            string fpath = join_path(snapshotDirForBestRun, oss.str());
            write_snapshot_path(fpath, cities, order, dist);
        }
    }

    return true;
}

// CLI arg helpers
static bool get_arg(int argc, char** argv, const string& key, string& val) {
    for (int i = 1; i + 1 < argc; i++) {
        if (string(argv[i]) == key) { val = argv[i + 1]; return true; }
    }
    return false;
}
static bool get_arg_int(int argc, char** argv, const string& key, int& val) {
    string s; if (!get_arg(argc, argv, key, s)) return false;
    val = stoi(s); return true;
}
static bool get_arg_ll(int argc, char** argv, const string& key, long long& val) {
    string s; if (!get_arg(argc, argv, key, s)) return false;
    val = stoll(s); return true;
}
static bool get_arg_double(int argc, char** argv, const string& key, double& val) {
    string s; if (!get_arg(argc, argv, key, s)) return false;
    val = stod(s); return true;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    // Modes:
    // 1) Single:
    // 2) Batch:

    // Snapshot files for GIF:
    //   output_bonus/dtX_frames/path_000100.txt, path_000200.txt, ...

    ENParams P;
    get_arg_int(argc, argv, "--runs", P.run_times);
    get_arg_ll(argc, argv, "--eval", P.eval_max);
    get_arg_int(argc, argv, "--record", P.record_every);

    get_arg_int(argc, argv, "--Mmul", P.M_mul);
    get_arg_int(argc, argv, "--Mcap", P.M_cap);

    get_arg_double(argc, argv, "--lr0", P.lr0);
    get_arg_double(argc, argv, "--lr1", P.lr_end);

    get_arg_double(argc, argv, "--sig0", P.sigma0);
    get_arg_double(argc, argv, "--sig1", P.sigma_end);

    get_arg_double(argc, argv, "--lam0", P.lambda0);
    get_arg_double(argc, argv, "--lam1", P.lambda_end);

    {
        int s; if (get_arg_int(argc, argv, "--seed", s)) P.seed = (unsigned)s;
    }

    if (argc >= 2 && string(argv[1]) == "--batch") {
        string datasetRoot, outDir;
        if (argc >= 3) datasetRoot = argv[2];
        else { cout << "Enter dataset root dir (e.g., dataset): "; cin >> datasetRoot; }

        if (argc >= 4) outDir = argv[3];
        else outDir = "output_bonus";

#ifdef _WIN32
        string cmd = "mkdir " + outDir + " >nul 2>nul";
#else
        string cmd = "mkdir -p " + outDir + " >/dev/null 2>&1";
#endif
        system(cmd.c_str());

        vector<string> dts = {"dt1", "dt2", "dt3", "dt4"}; // dt4 running longer but included in bonus

        for (auto& dt : dts) {
            string pointPath = join_path(join_path(datasetRoot, dt), "point.txt");
            string outPath   = join_path(outDir, "ans_" + dt + ".txt");
            string framesDir = join_path(outDir, dt + "_frames");

            vector<City> cities;
            if (!read_points(pointPath, cities)) {
                cerr << "[BONUS][SKIP] Cannot read " << pointPath << "\n";
                continue;
            }

            ENParams P2 = P;
            if (P2.eval_max <= 0) P2.eval_max = 10000LL * (long long)cities.size();

            double mean_best, best_all;
            vector<int> best_order;
            bool ok = solve_elastic_net_bonus(cities, P2, mean_best, best_all, best_order, framesDir);
            if (!ok) {
                cerr << "[BONUS][FAIL] " << dt << ": Elastic Net failed.\n";
                continue;
            }

            if (!write_summary_ans(outPath, mean_best, best_all, cities, best_order)) {
                cerr << "[BONUS][FAIL] Cannot write " << outPath << "\n";
                continue;
            }

            cout << "[BONUS][OK] " << dt
                 << " n=" << cities.size()
                 << " meanBest=" << fixed << setprecision(3) << mean_best
                 << " best=" << best_all
                 << " -> " << outPath
                 << " | frames: " << framesDir << "\n";
        }

        cout << "[BONUS] Batch done. (dt1~dt4)\n";
        return 0;
    }

    // Single-file mode
    string inputFile, outputFile;
    if (argc >= 2) inputFile = argv[1];
    else { cout << "Enter input point file path: "; cin >> inputFile; }

    if (argc >= 3) outputFile = argv[2];
    else { cout << "Enter output file path: "; cin >> outputFile; }

    vector<City> cities;
    if (!read_points(inputFile, cities)) {
        cerr << "ERROR: Cannot read input file: " << inputFile << "\n";
        return 1;
    }

    if (P.eval_max <= 0) P.eval_max = 10000LL * (long long)cities.size();

    // frames in single mode: outputFile + "_frames"
    string framesDir = outputFile + "_frames";

    double mean_best, best_all;
    vector<int> best_order;
    if (!solve_elastic_net_bonus(cities, P, mean_best, best_all, best_order, framesDir)) {
        cerr << "ERROR: Elastic Net failed.\n";
        return 1;
    }

    if (!write_summary_ans(outputFile, mean_best, best_all, cities, best_order)) {
        cerr << "ERROR: Cannot write output file: " << outputFile << "\n";
        return 1;
    }

    cout << "Solved with Elastic Net (Bonus).\n";
    cout << "Cities: " << cities.size() << "\n";
    cout << "Mean best (30 runs): " << fixed << setprecision(3) << mean_best << "\n";
    cout << "Best overall: " << best_all << "\n";
    cout << "Output: " << outputFile << "\n";
    cout << "Frames: " << framesDir << " (path_******.txt every " << P.record_every << " eval)\n";
    return 0;
}
