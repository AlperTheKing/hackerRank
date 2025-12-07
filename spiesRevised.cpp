#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <numeric>
#include <random>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

struct LineKey {
    int a{0};
    int b{0};
    int c{0};
    bool operator==(const LineKey &other) const noexcept {
        return a == other.a && b == other.b && c == other.c;
    }
};

struct LineHash {
    std::size_t operator()(const LineKey &key) const noexcept {
        std::size_t h1 = std::hash<int>{}(key.a);
        std::size_t h2 = std::hash<int>{}(key.b);
        std::size_t h3 = std::hash<int>{}(key.c);
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};

static LineKey makeLineKey(int x1, int y1, int x2, int y2) {
    long long A = static_cast<long long>(y2) - static_cast<long long>(y1);
    long long B = static_cast<long long>(x1) - static_cast<long long>(x2);
    long long C = static_cast<long long>(x2) * y1 - static_cast<long long>(x1) * y2;
    long long g = std::gcd(std::gcd(std::llabs(A), std::llabs(B)), std::llabs(C));
    if (g == 0) g = 1;
    A /= g;
    B /= g;
    C /= g;
    if (A < 0 || (A == 0 && B < 0) || (A == 0 && B == 0 && C < 0)) {
        A = -A;
        B = -B;
        C = -C;
    }
    return {static_cast<int>(A), static_cast<int>(B), static_cast<int>(C)};
}

struct ValidationResult {
    bool ok;
    std::string message;
};

static ValidationResult validatePlacement(const std::vector<int> &rowsPerColumn) {
    const int n = static_cast<int>(rowsPerColumn.size());
    std::vector<char> rowUsed(n, 0), colUsed(n, 0);
    std::vector<char> diagMain(2 * n + 3, 0), diagAnti(2 * n + 3, 0);

    std::vector<std::pair<int, int>> processed;
    processed.reserve(n);
    std::unordered_map<LineKey, int, LineHash> linePairs;
    linePairs.reserve(static_cast<std::size_t>(n) * (n - 1) / 2);

    for (int col = 0; col < n; ++col) {
        int row = rowsPerColumn[col];
        if (row < 0 || row >= n) {
            return {false, "Row out of bounds at column " + std::to_string(col + 1)};
        }
        if (rowUsed[row]) return {false, "Row conflict at row " + std::to_string(row + 1)};
        if (colUsed[col]) return {false, "Column conflict at column " + std::to_string(col + 1)};
        rowUsed[row] = colUsed[col] = 1;

        int mainIdx = row + col;
        int antiIdx = row - col + n;
        if (diagMain[mainIdx]) return {false, "Main diagonal conflict"};
        if (diagAnti[antiIdx]) return {false, "Anti-diagonal conflict"};
        diagMain[mainIdx] = diagAnti[antiIdx] = 1;

        for (const auto &prev : processed) {
            LineKey key = makeLineKey(row, col, prev.first, prev.second);
            auto it = linePairs.find(key);
            if (it != linePairs.end() && it->second > 0) {
                return {false, "Collinear trio detected"};
            }
        }
        for (const auto &prev : processed) {
            LineKey key = makeLineKey(row, col, prev.first, prev.second);
            linePairs[key] += 1;
        }
        processed.emplace_back(row, col);
    }
    return {true, "ok"};
}

class Solver {
public:
    Solver(int n, std::uint64_t seed) : n(n), rng(seed) {
        initializeStructures();
        seedInitialPlacement();
        initialized = true;
    }

    bool solve(int maxIter = -1, std::atomic<bool> *stop = nullptr) {
        int iter = 0;
        while (!solved()) {
            if (stop && stop->load(std::memory_order_relaxed)) {
                return false;
            }
            if (maxIter >= 0 && iter == maxIter) {
                return false;
            }

            calculatePlacementOptions();
            if (changableColumns.empty()) {
                return false;
            }

            int col = getChangableCol();
            revalidate(col);
            int row = getPlacement(col);
            rows[col] = row;
            invalidate(col);
            ++iter;
        }
        return true;
    }

    std::vector<int> getSolution() const { return rows; }

private:
    int n;
    std::mt19937 rng;

    std::vector<std::vector<int>> numConflicts;
    std::vector<std::vector<std::pair<int, int>>> simpleInvalidations;
    std::vector<std::vector<std::vector<std::pair<int, int>>>> colinearInvalidations;

    std::vector<int> rows;
    std::vector<std::vector<int>> placementOptions;
    std::vector<int> changableColumns;
    bool initialized{false};

    void initializeStructures() {
        numConflicts.assign(n, std::vector<int>(n, 0));
        simpleInvalidations.assign(n, std::vector<std::pair<int, int>>());
        colinearInvalidations.assign(n,
                                     std::vector<std::vector<std::pair<int, int>>>(n));
        placementOptions.assign(n, std::vector<int>());
        rows.assign(n, n * 2 - 1);
    }

    void seedInitialPlacement() {
        std::uniform_int_distribution<int> dist(0, n - 1);
        for (int col = 0; col < n; col++) {
            rows[col] = dist(rng);
            revalidate(col);
            invalidate(col);
        }
    }

    void invalidateSpot(int row, int col, int originalColumn) {
        numConflicts[row][col] += 1;
        simpleInvalidations[originalColumn].push_back(std::make_pair(row, col));
    }

    void invalidateSimple(int col) {
        int row = rows[col];
        for (int i = 0; i < col; i++) {
            invalidateSpot(row, i, col);
        }
        for (int i = col + 1; i < n; i++) {
            invalidateSpot(row, i, col);
        }
        for (int i = 1; row + i < n && col + i < n; i++) {
            invalidateSpot(row + i, col + i, col);
        }
        for (int i = 1; 0 <= row - i && col + i < n; i++) {
            invalidateSpot(row - i, col + i, col);
        }
        for (int i = 1; row + i < n && 0 <= col - i; i++) {
            invalidateSpot(row + i, col - i, col);
        }
        for (int i = 1; 0 <= row - i && 0 <= col - i; i++) {
            invalidateSpot(row - i, col - i, col);
        }
    }

    void invalidateColinearPair(int col1, int col2) {
        int row1 = rows[col1];
        int row2 = rows[col2];
        int rowStep = row2 - row1;
        int colStep = col2 - col1;
        int g = std::abs(std::gcd(rowStep, colStep));
        rowStep /= g;
        colStep /= g;

        if (row1 == n * 2 - 1) {
            if (row2 != n * 2 - 1) {
                numConflicts[row2][col2] += 1;
                colinearInvalidations[col1][col2].push_back(std::make_pair(row2, col2));
            }
            return;
        }

        int row;
        int col;
        if (rowStep > 0) {
            row = row1;
            col = col1;
            while (col < n && row < n) {
                numConflicts[row][col] += 1;
                colinearInvalidations[col1][col2].push_back(std::make_pair(row, col));
                row += rowStep;
                col += colStep;
            }
            row = row1 - rowStep;
            col = col1 - colStep;
            while (0 <= col && 0 <= row) {
                numConflicts[row][col] += 1;
                colinearInvalidations[col1][col2].push_back(std::make_pair(row, col));
                row -= rowStep;
                col -= colStep;
            }
        } else {
            row = row1;
            col = col1;
            while (col < n && 0 <= row) {
                numConflicts[row][col] += 1;
                colinearInvalidations[col1][col2].push_back(std::make_pair(row, col));
                row += rowStep;
                col += colStep;
            }
            row = row1 - rowStep;
            col = col1 - colStep;
            while (0 <= col && row < n) {
                numConflicts[row][col] += 1;
                colinearInvalidations[col1][col2].push_back(std::make_pair(row, col));
                row -= rowStep;
                col -= colStep;
            }
        }
    }

    void invalidateColinear(int col) {
        for (int i = 0; i < col; i++) {
            invalidateColinearPair(i, col);
        }
        for (int i = col + 1; i < n; i++) {
            invalidateColinearPair(col, i);
        }
        numConflicts[rows[col]][col] -= n - 1;
    }

    void invalidate(int col) {
        invalidateSimple(col);
        invalidateColinear(col);
    }

    void revalidateSimple(int col) {
        for (const auto &spot : simpleInvalidations[col]) {
            numConflicts[spot.first][spot.second] -= 1;
        }
        simpleInvalidations[col].clear();
    }

    void revalidateColinear(int col) {
        for (int i = 0; i < col; i++) {
            for (auto spot : colinearInvalidations[i][col]) {
                numConflicts[spot.first][spot.second] -= 1;
            }
            colinearInvalidations[i][col].clear();
        }

        for (int i = col + 1; i < n; i++) {
            for (auto spot : colinearInvalidations[col][i]) {
                numConflicts[spot.first][spot.second] -= 1;
            }
            colinearInvalidations[col][i].clear();
        }

        if (initialized) {
            numConflicts[rows[col]][col] += n - 1;
        }
    }

    void revalidate(int col) {
        revalidateSimple(col);
        revalidateColinear(col);
    }

    void calculatePlacementOptions(int col) {
        placementOptions[col].clear();
        placementOptions[col].push_back(0);
        int minConflicts = numConflicts[0][col];
        for (int row = 1; row < n; row++) {
            int conflicts = numConflicts[row][col];
            if (conflicts < minConflicts) {
                placementOptions[col].clear();
                minConflicts = conflicts;
                placementOptions[col].push_back(row);
            } else if (conflicts == minConflicts) {
                placementOptions[col].push_back(row);
            }
        }
    }

    void calculatePlacementOptions() {
        changableColumns.clear();
        for (int col = 0; col < n; col++) {
            calculatePlacementOptions(col);
            if (placementOptions[col].size() > 1 || placementOptions[col][0] != rows[col]) {
                changableColumns.push_back(col);
            }
        }
    }

    int getChangableCol() {
        std::uniform_int_distribution<int> dist(0,
                                                static_cast<int>(changableColumns.size()) - 1);
        return changableColumns[dist(rng)];
    }

    int getPlacement(int col) {
        std::uniform_int_distribution<int> dist(0,
                                                static_cast<int>(placementOptions[col].size()) - 1);
        return placementOptions[col][dist(rng)];
    }

    bool solved() {
        for (int col = 0; col < n; col++) {
            if (numConflicts[rows[col]][col] != 0) {
                return false;
            }
        }
        return true;
    }
};

struct RunConfig {
    int n = 999;
    int maxIter = -1;
    int threads = 0;
};

static RunConfig parseArgs(int argc, char *argv[]) {
    RunConfig cfg;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--threads" && i + 1 < argc) {
            cfg.threads = std::max(1, std::stoi(argv[++i]));
        } else if (arg == "--max-iter" && i + 1 < argc) {
            cfg.maxIter = std::stoi(argv[++i]);
        } else {
            try {
                cfg.n = std::stoi(arg);
            } catch (...) {
                std::cerr << "Ignoring invalid argument: " << arg << '\n';
            }
        }
    }
    if (cfg.n < 4) {
        cfg.n = 4;
    }
    if (cfg.threads <= 0) {
        unsigned int hw = std::thread::hardware_concurrency();
        cfg.threads = hw ? static_cast<int>(hw) : 2;
    }
    return cfg;
}

struct SharedState {
    std::atomic<bool> done{false};
    std::vector<int> solution;
    std::mutex mtx;
};

static void worker(int id, const RunConfig &cfg, SharedState &state) {
    std::mt19937_64 seeder(std::random_device{}() ^ (static_cast<std::uint64_t>(id) << 16));
    while (!state.done.load(std::memory_order_relaxed)) {
        Solver solver(cfg.n, seeder());
        bool ok = solver.solve(cfg.maxIter, &state.done);
        if (!ok) continue;

        auto sol = solver.getSolution();
        auto res = validatePlacement(sol);
        if (!res.ok) continue;

        bool expected = false;
        if (state.done.compare_exchange_strong(expected, true, std::memory_order_acq_rel)) {
            std::lock_guard<std::mutex> lock(state.mtx);
            state.solution = std::move(sol);
        }
        return;
    }
}

int main(int argc, char *argv[]) {
    RunConfig cfg = parseArgs(argc, argv);
    SharedState state;
    state.solution.assign(cfg.n, -1);

    std::vector<std::thread> threads;
    threads.reserve(static_cast<std::size_t>(cfg.threads));

    auto start = std::chrono::steady_clock::now();
    for (int i = 0; i < cfg.threads; ++i) {
        threads.emplace_back(worker, i, std::cref(cfg), std::ref(state));
    }
    for (auto &t : threads) t.join();
    auto stop = std::chrono::steady_clock::now();
    auto elapsedMs = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

    if (!state.done.load()) {
        std::cerr << "Failed to find a solution for N=" << cfg.n << '\n';
        return 1;
    }

    std::vector<int> colsPerRow(cfg.n, -1);
    for (int col = 0; col < cfg.n; ++col) {
        int row = state.solution[col];
        if (row >= 0 && row < cfg.n) {
            colsPerRow[row] = col;
        }
    }

    std::cerr << "Solved N=" << cfg.n << " in " << elapsedMs << " ms using "
              << cfg.threads << " threads\n";

    std::cout << cfg.n << '\n';
    for (int row = 0; row < cfg.n; ++row) {
        if (row) std::cout << ' ';
        std::cout << colsPerRow[row] + 1;
    }
    std::cout << '\n';
    return 0;
}
