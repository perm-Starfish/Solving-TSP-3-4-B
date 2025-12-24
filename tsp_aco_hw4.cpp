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

// Compute tour length (order holds indices, does NOT repeat start at end)
static double tour_length(const vector<int>& order, const vector<vector<double>>& d) {
    int n = (int)order.size();
    double L = 0.0;
    for (int i = 0; i < n - 1; i++) L += d[order[i]][order[i + 1]];
    L += d[order.back()][order[0]]; // close tour
    return L;
}

// Optional 2-opt local search (kept OFF by default; can enable via flag)
// 2-opt switches edges to reduce tour length
static void two_opt(vector<int>& order, const vector<vector<double>>& d) {
    int n = (int)order.size();
    bool improved = true;
    while (improved) {
        improved = false;
        for (int i = 1; i < n - 2; i++) {
            for (int k = i + 1; k < n - 1; k++) {
                int a = order[i - 1], b = order[i];
                int c = order[k], d2 = order[k + 1];
                double before = d[a][b] + d[c][d2];
                double after  = d[a][c] + d[b][d2];
                if (after + 1e-12 < before) {
                    reverse(order.begin() + i, order.begin() + k + 1);
                    improved = true;
                }
            }
        }
    }
}

struct ACOParams {
    int run_times = 30;
    int iterations = 200;           // can set
    int population_size = 50;       // can set, the number of ants (paths) per iteration
    double alpha = 1.0;             // pheromone importance
    double beta  = 5.0;             // heuristic importance
    double rho   = 0.5;             // evaporation rate, avoid early convergence
    double Q     = 100.0;           // constant
    long long eval_max = -1;        // default: 10000 * n
    bool use_2opt = false;          // optional local search
    unsigned seed = 0;              // 0 => random_device
};

// Construct one ant solution (tour)
static vector<int> construct_ant_solution(
    int n,
    const vector<vector<double>>& tau,
    const vector<vector<double>>& eta,
    mt19937& rng
) {
    vector<int> tour;
    tour.reserve(n);

    // start from node 0
    int current = 0;
    tour.push_back(current);

    // visited set
    vector<char> visited(n, false);
    visited[current] = true;

    for (int step = 1; step < n; step++) {
        // build probability for choosing next city j
        vector<double> prob(n, 0.0);
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (visited[j]) continue;
            // prob ∝ (tau)^alpha * (eta)^beta
            // eta = 1 / distance
            double p = tau[current][j] * eta[current][j];
            // tau already stores tau^alpha.
            // will compute exactly as formula.
            p = pow(tau[current][j], 1.0) * pow(eta[current][j], 1.0); // placeholder; will overwrite below
        }
        // will re-loop once
        // NOTE: This function doesn't know alpha/beta; will compute in caller by passing pre-powered tables.
        // So in this helper, tau and eta are already pre-powered.

        // HERE: pass tauPow and etaPow in. So use directly:
        sum = 0.0;
        for (int j = 0; j < n; j++) {
            if (visited[j]) continue;
            double p = tau[current][j] * eta[current][j];
            prob[j] = p; // store p of unvisited j
            sum += p;
        }

        int next = -1;
        if (sum <= 0.0) {
            // fallback: choose any unvisited uniformly
            vector<int> cand;
            for (int j = 0; j < n; j++) if (!visited[j]) cand.push_back(j);
            uniform_int_distribution<int> uni(0, (int)cand.size() - 1);
            next = cand[uni(rng)]; // pick random
        } else {
            // roulette wheel
            uniform_real_distribution<double> uni(0.0, sum);
            double r = uni(rng); // random threshold
            double acc = 0.0;
            for (int j = 0; j < n; j++) {
                if (visited[j]) continue;
                acc += prob[j];
                if (r <= acc) { next = j; break; } // greater than threshold -> select as next
            }
            if (next == -1) {
                // fallback: numerical edge-case
                for (int j = n - 1; j >= 0; j--) if (!visited[j]) { next = j; break; }
            }
        }

        tour.push_back(next);
        visited[next] = true;
        current = next;
    }

    return tour;
}

static bool solve_tsp_aco(
    const vector<City>& cities,
    const ACOParams& P,
    double& mean_best_30,
    double& best_overall,
    vector<int>& best_order_idx,      // indices in cities
    vector<double>* best_record_100   // optional: record best every 100 eval for GIF
) {
    int n = (int)cities.size();
    if (n < 2) return false;

    // distances + eta
    vector<vector<double>> d(n, vector<double>(n, 0.0));
    vector<vector<double>> eta_raw(n, vector<double>(n, 0.0)); // 1/distance
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            d[i][j] = dist2d(cities[i], cities[j]);
            if (i == j) eta_raw[i][j] = 0.0;
            else eta_raw[i][j] = 1.0 / max(d[i][j], 1e-12); // eta = 1/distance, be selected more if closer
        }
    }

    long long eval_max = P.eval_max;
    if (eval_max <= 0) eval_max = 10000LL * n;

    // Random Number Generator
    mt19937 rng;
    if (P.seed == 0) rng.seed(random_device{}());
    else rng.seed(P.seed);

    vector<double> best_each_run;
    best_each_run.reserve(P.run_times);

    best_overall = 1e100;
    best_order_idx.clear();

    if (best_record_100) best_record_100->clear();

    for (int run = 0; run < P.run_times; run++) {
        // initialize pheromone trails tau (raw, not powered)
        // common choice: small constant
        vector<vector<double>> tau(n, vector<double>(n, 1.0)); // initialize tau

        double best_run = 1e100;
        vector<int> best_run_order;

        long long eval_count = 0;

        // allow stopping by eval_count even before iterations end
        for (int it = 0; it < P.iterations && eval_count < eval_max; it++) {
            // Precompute tau^alpha and eta^beta tables for speed & correctness
            vector<vector<double>> tauPow(n, vector<double>(n, 0.0));
            vector<vector<double>> etaPow(n, vector<double>(n, 0.0));
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    tauPow[i][j] = pow(max(tau[i][j], 1e-12), P.alpha);
                    etaPow[i][j] = pow(max(eta_raw[i][j], 1e-12), P.beta);
                }
            }
            // outputs tauPow & etaPow

            // Each iteration: each ant constructs one tour => population_size evaluations
            vector<vector<int>> ant_tours;
            ant_tours.reserve(P.population_size);

            vector<double> ant_lengths;
            ant_lengths.reserve(P.population_size);

            // ant tours and lengths
            for (int a = 0; a < P.population_size && eval_count < eval_max; a++) {
                vector<int> tour = construct_ant_solution(n, tauPow, etaPow, rng);

                if (P.use_2opt) two_opt(tour, d);

                double L = tour_length(tour, d);

                ant_tours.push_back(move(tour));
                ant_lengths.push_back(L); // update pheromone deposit = Q / L

                eval_count++; // an ant create one tour = one eval

                if (L < best_run) {
                    best_run = L; // best tour length in this run
                    best_run_order = ant_tours.back();
                }
                if (L < best_overall) {
                    best_overall = L; // best tour length overall (all runs)
                    best_order_idx = ant_tours.back();
                }

                // record best every 100 eval, no need
                if (best_record_100 && (eval_count % 100 == 0)) {
                    best_record_100->push_back(best_run);
                }
            }

            // Update pheromones:
            // tau_ij = (1-rho)*tau_ij + sum_k delta_tau_ij^k
            // delta_tau_ij^k = Q / L_k if ant k used edge (i,j) else 0
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    tau[i][j] *= (1.0 - P.rho);

            for (int a = 0; a < (int)ant_tours.size(); a++) {
                double Lk = ant_lengths[a];
                double deposit = P.Q / max(Lk, 1e-12);
                const auto& tour = ant_tours[a];
                for (int i = 0; i < n - 1; i++) {
                    int u = tour[i], v = tour[i + 1];
                    tau[u][v] += deposit;
                    tau[v][u] += deposit;
                }
                // closing edge
                int u = tour.back(), v = tour[0];
                tau[u][v] += deposit;
                tau[v][u] += deposit; // ant that walk shorter path deposit more pheromone
            }

            // (optional) can cap tau to avoid explosion
            const double TAU_MAX = 1e6;
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    tau[i][j] = min(tau[i][j], TAU_MAX);
        }

        best_each_run.push_back(best_run);
    }

    // mean best over 30 runs
    double sum = 0.0;
    for (double v : best_each_run) sum += v;
    mean_best_30 = sum / (double)best_each_run.size();

    return true;
}

static bool write_answer_hw4(const string& outPath,
                             double mean_best_30,
                             double best_overall,
                             const vector<City>& cities,
                             const vector<int>& best_order_idx) {
    ofstream fout(outPath);
    if (!fout) return false;
    fout << fixed << setprecision(3);
    fout << "mean distance: " << mean_best_30 << "\n";
    fout << "distance: " << best_overall << "\n";
    for (int idx : best_order_idx) fout << cities[idx].id << "\n";
    return true;
}

// Simple CLI parsing helper
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

    ACOParams P;

    // optional parameter overrides from CLI
    get_arg_int(argc, argv, "--runs", P.run_times);
    get_arg_int(argc, argv, "--iter", P.iterations);
    get_arg_int(argc, argv, "--pop",  P.population_size);
    get_arg_double(argc, argv, "--alpha", P.alpha);
    get_arg_double(argc, argv, "--beta",  P.beta);
    get_arg_double(argc, argv, "--rho",   P.rho);
    get_arg_double(argc, argv, "--Q",     P.Q);
    get_arg_ll(argc, argv, "--eval",      P.eval_max);
    {
        int s; if (get_arg_int(argc, argv, "--seed", s)) P.seed = (unsigned)s;
    }
    {
        string flag;
        if (get_arg(argc, argv, "--2opt", flag)) {
            // any value means enable, e.g. --2opt 1
            P.use_2opt = true;
        }
    }

    // enforce run_times=30 (but allow override)
    // P.run_times = 30;

    if (argc >= 2 && string(argv[1]) == "--batch") {
        string datasetRoot, outDir;
        if (argc >= 3) datasetRoot = argv[2];
        else { cout << "Enter dataset root dir (e.g., dataset): "; cin >> datasetRoot; }
        if (argc >= 4) outDir = argv[3];
        else outDir = "output_hw4";

#ifdef _WIN32
        string cmd = "mkdir " + outDir + " >nul 2>nul";
#else
        string cmd = "mkdir -p " + outDir + " >/dev/null 2>&1";
#endif
        system(cmd.c_str());

        vector<string> dts = {"dt1", "dt2", "dt3", "dt4"};
        for (auto& dt : dts) {
            string pointPath = join_path(join_path(datasetRoot, dt), "point.txt");
            string outPath   = join_path(outDir, "ans_" + dt + ".txt");

            vector<City> cities;
            if (!read_points(pointPath, cities)) {
                cerr << "[HW4][SKIP] Cannot read " << pointPath << "\n";
                continue;
            }

            // If eval_max not specified, use 10000 * n per run
            ACOParams P2 = P;
            if (P2.eval_max <= 0) P2.eval_max = 10000LL * (long long)cities.size();

            double mean_best, best_all;
            vector<int> best_order;
            bool ok = solve_tsp_aco(cities, P2, mean_best, best_all, best_order, nullptr);
            if (!ok) {
                cerr << "[HW4][FAIL] " << dt << ": ACO failed.\n";
                continue;
            }

            if (!write_answer_hw4(outPath, mean_best, best_all, cities, best_order)) {
                cerr << "[HW4][FAIL] Cannot write " << outPath << "\n";
                continue;
            }

            cout << "[HW4][OK] " << dt
                 << "  n=" << cities.size()
                 << "  meanBest=" << fixed << setprecision(3) << mean_best
                 << "  best=" << best_all
                 << "  -> " << outPath << "\n";
        }

        cout << "[HW4] Batch done. (dt1~dt4)\n";
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

    double mean_best, best_all;
    vector<int> best_order;
    if (!solve_tsp_aco(cities, P, mean_best, best_all, best_order, nullptr)) {
        cerr << "ERROR: ACO failed.\n";
        return 1;
    }

    if (!write_answer_hw4(outputFile, mean_best, best_all, cities, best_order)) {
        cerr << "ERROR: Cannot write output file: " << outputFile << "\n";
        return 1;
    }

    cout << "Solved with ACO.\n";
    cout << "Cities: " << cities.size() << "\n";
    cout << "Mean best (30 runs): " << fixed << setprecision(3) << mean_best << "\n";
    cout << "Best overall: " << best_all << "\n";
    cout << "Output: " << outputFile << "\n";
    return 0;
}
