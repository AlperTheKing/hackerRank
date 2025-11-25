#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <mutex>
#include <numeric>
#include <random>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

struct Position {
    int row;  // 1-indexed
    int col;  // 1-indexed
    std::string notation;
};

std::string columnToLetters(int col) {
    std::string letters;
    while (col > 0) {
        col--;
        letters.push_back(static_cast<char>('a' + (col % 26)));
        col /= 26;
    }
    std::reverse(letters.begin(), letters.end());
    return letters;
}

std::string toChess(int row, int col) {
    return columnToLetters(col) + std::to_string(row);
}

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

LineKey makeLineKey(int x1, int y1, int x2, int y2) {
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

ValidationResult validatePermutation(const std::vector<int> &perm) {
    const int n = static_cast<int>(perm.size());
    std::vector<Position> positions;
    positions.reserve(n);
    for (int row = 0; row < n; ++row) {
        int col = perm[row];
        positions.push_back({row + 1, col + 1, toChess(row + 1, col + 1)});
    }

    std::sort(positions.begin(), positions.end(), [](const Position &lhs, const Position &rhs) {
        return lhs.notation < rhs.notation;
    });

    std::vector<char> rowUsed(n + 1, 0);
    std::vector<char> colUsed(n + 1, 0);
    std::vector<char> diagMain(2 * n + 3, 0);
    std::vector<char> diagAnti(2 * n + 3, 0);

    std::vector<Position> processed;
    processed.reserve(n);
    std::unordered_map<LineKey, int, LineHash> linePairs;
    linePairs.reserve(static_cast<std::size_t>(n) * (n - 1) / 2);

    for (const auto &pos : positions) {
        if (rowUsed[pos.row]) {
            return {false, "Row conflict detected while scanning " + pos.notation};
        }
        rowUsed[pos.row] = 1;

        if (colUsed[pos.col]) {
            return {false, "Column conflict detected while scanning " + pos.notation};
        }
        colUsed[pos.col] = 1;

        int mainIdx = pos.row + pos.col;
        if (diagMain[mainIdx]) {
            return {false, "Main diagonal conflict triggered by " + pos.notation};
        }
        diagMain[mainIdx] = 1;

        int antiIdx = pos.row - pos.col + n;
        if (diagAnti[antiIdx]) {
            return {false, "Anti-diagonal conflict triggered by " + pos.notation};
        }
        diagAnti[antiIdx] = 1;

        for (const auto &prev : processed) {
            LineKey key = makeLineKey(pos.row, pos.col, prev.row, prev.col);
            auto it = linePairs.find(key);
            if (it != linePairs.end() && it->second > 0) {
                return {false, "Collinear trio detected via " + prev.notation + " and " + pos.notation};
            }
        }

        for (const auto &prev : processed) {
            LineKey key = makeLineKey(pos.row, pos.col, prev.row, prev.col);
            linePairs[key] += 1;
        }
        processed.push_back(pos);
    }
    return {true, "Validated"};
}

long long pairPenalty(int pairs) {
    return (pairs >= 3) ? (pairs - 2) : 0;
}

void incrementLine(std::unordered_map<LineKey, int, LineHash> &linePairs,
                   const LineKey &key, long long &cost) {
    auto it = linePairs.find(key);
    int before = 0;
    if (it == linePairs.end()) {
        linePairs.emplace(key, 1);
    } else {
        before = it->second;
        it->second = before + 1;
    }
    int after = before + 1;
    cost += pairPenalty(after) - pairPenalty(before);
}

void decrementLine(std::unordered_map<LineKey, int, LineHash> &linePairs,
                   const LineKey &key, long long &cost) {
    auto it = linePairs.find(key);
    if (it == linePairs.end()) return;
    int before = it->second;
    int after = before - 1;
    cost += pairPenalty(after) - pairPenalty(before);
    if (after == 0) {
        linePairs.erase(it);
    } else {
        it->second = after;
    }
}

void initializeLineStructures(const std::vector<int> &perm,
                              std::vector<LineKey> &pairKeys,
                              std::unordered_map<LineKey, int, LineHash> &linePairs,
                              long long &cost) {
    const int n = static_cast<int>(perm.size());
    pairKeys.assign(static_cast<std::size_t>(n) * n, LineKey{});
    linePairs.clear();
    linePairs.reserve(static_cast<std::size_t>(n) * (n - 1) / 2);
    cost = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            LineKey key = makeLineKey(i + 1, perm[i] + 1, j + 1, perm[j] + 1);
            std::size_t idx = static_cast<std::size_t>(i) * n + j;
            pairKeys[idx] = key;
            incrementLine(linePairs, key, cost);
        }
    }
}

struct PairDelta {
    std::size_t index;
    LineKey oldKey;
    LineKey newKey;
};

struct SearchConfig {
    int maxRestarts;
    int maxIterations;
};

int swapScore(int n, int r1, int c1, int r2, int c2,
              const std::vector<int> &diagMain, const std::vector<int> &diagAnti) {
    const int offset = n - 1;
    auto diagIdx = [](int row, int col) { return row + col; };
    auto antiIdx = [offset](int row, int col) { return row - col + offset; };

    auto conflictsAfter = [&](int row, int newCol) {
        int dMain = diagIdx(row, newCol);
        int dAnti = antiIdx(row, newCol);
        int conflicts = diagMain[dMain] + diagAnti[dAnti];
        if (dMain == diagIdx(r1, c1)) conflicts--;
        if (dMain == diagIdx(r2, c2)) conflicts--;
        if (dAnti == antiIdx(r1, c1)) conflicts--;
        if (dAnti == antiIdx(r2, c2)) conflicts--;
        return conflicts;
    };

    return conflictsAfter(r1, c2) + conflictsAfter(r2, c1);
}

bool attemptSolve(std::vector<int> perm, std::mt19937_64 &rng, int maxIterations, std::vector<int> &out) {
    const int n = static_cast<int>(perm.size());
    std::vector<int> diagMain(2 * n, 0), diagAnti(2 * n, 0);
    auto diagIdx = [](int row, int col) { return row + col; };
    auto antiIdx = [offset = n - 1](int row, int col) { return row - col + offset; };

    auto rebuildDiagonals = [&]() {
        std::fill(diagMain.begin(), diagMain.end(), 0);
        std::fill(diagAnti.begin(), diagAnti.end(), 0);
        for (int row = 0; row < n; ++row) {
            diagMain[diagIdx(row, perm[row])]++;
            diagAnti[antiIdx(row, perm[row])]++;
        }
    };

    rebuildDiagonals();

    std::vector<int> rowConflicts(n, 0);
    auto recomputeConflicts = [&]() {
        for (int row = 0; row < n; ++row) {
            int mainConf = diagMain[diagIdx(row, perm[row])] - 1;
            int antiConf = diagAnti[antiIdx(row, perm[row])] - 1;
            rowConflicts[row] = mainConf + antiConf;
        }
    };

    recomputeConflicts();
    int stagnation = 0;

    for (int iter = 0; iter < maxIterations; ++iter) {
        std::vector<int> conflictRows;
        conflictRows.reserve(n);
        for (int row = 0; row < n; ++row) {
            if (rowConflicts[row] > 0) conflictRows.push_back(row);
        }
        if (conflictRows.empty()) {
            out = perm;
            return true;
        }

        int r1 = conflictRows[rng() % conflictRows.size()];
        int c1 = perm[r1];
        int bestRow = -1;
        int bestScore = std::numeric_limits<int>::max();
        for (int r2 = 0; r2 < n; ++r2) {
            if (r2 == r1) continue;
            int candidateScore = swapScore(n, r1, c1, r2, perm[r2], diagMain, diagAnti);
            if (candidateScore < bestScore) {
                bestScore = candidateScore;
                bestRow = r2;
            }
        }

        if (bestRow == -1) {
            std::shuffle(perm.begin(), perm.end(), rng);
            rebuildDiagonals();
            recomputeConflicts();
            stagnation = 0;
            continue;
        }

        int previousCombined = rowConflicts[r1] + rowConflicts[bestRow];
        int r2 = bestRow;
        int c2 = perm[r2];

        diagMain[diagIdx(r1, c1)]--;
        diagAnti[antiIdx(r1, c1)]--;
        diagMain[diagIdx(r2, c2)]--;
        diagAnti[antiIdx(r2, c2)]--;

        std::swap(perm[r1], perm[r2]);
        c1 = perm[r1];
        c2 = perm[r2];

        diagMain[diagIdx(r1, c1)]++;
        diagAnti[antiIdx(r1, c1)]++;
        diagMain[diagIdx(r2, c2)]++;
        diagAnti[antiIdx(r2, c2)]++;

        recomputeConflicts();

        if (bestScore >= previousCombined) {
            stagnation++;
        } else {
            stagnation = 0;
        }
        if (stagnation > n) {
            std::shuffle(perm.begin(), perm.end(), rng);
            rebuildDiagonals();
            recomputeConflicts();
            stagnation = 0;
        }
    }
    return false;
}

bool solveQueens(int n, std::mt19937_64 &rng, const SearchConfig &config, std::vector<int> &out) {
    std::vector<int> base(n);
    std::iota(base.begin(), base.end(), 0);
    for (int attempt = 0; attempt < config.maxRestarts; ++attempt) {
        std::shuffle(base.begin(), base.end(), rng);
        std::vector<int> candidate = base;
        if (attemptSolve(candidate, rng, config.maxIterations, out)) {
            return true;
        }
    }
    return false;
}

bool enforceCollinearity(std::vector<int> &perm, std::mt19937_64 &rng, int multiplier = 40) {
    const int n = static_cast<int>(perm.size());
    std::vector<int> diagMain(2 * n, -1);
    std::vector<int> diagAnti(2 * n, -1);
    auto mainIdx = [](int row, int col) { return row + col; };
    auto antiIdx = [offset = n - 1](int row, int col) { return row - col + offset; };
    for (int row = 0; row < n; ++row) {
        diagMain[mainIdx(row, perm[row])] = row;
        diagAnti[antiIdx(row, perm[row])] = row;
    }

    std::vector<int> colOwner(n, -1);
    for (int row = 0; row < n; ++row) {
        colOwner[perm[row]] = row;
    }

    std::vector<LineKey> pairKeys;
    std::unordered_map<LineKey, int, LineHash> linePairs;
    long long cost = 0;
    initializeLineStructures(perm, pairKeys, linePairs, cost);
    if (cost == 0) {
        return true;
    }

    const long long maxSteps = static_cast<long long>(n) * n * multiplier;
    std::uniform_real_distribution<double> prob(0.0, 1.0);

    auto diagFreeAfterSwap = [&](int row, int newCol, int allowedRow) {
        int idxMain = mainIdx(row, newCol);
        int idxAnti = antiIdx(row, newCol);
        if (diagMain[idxMain] != -1 && diagMain[idxMain] != allowedRow) return false;
        if (diagAnti[idxAnti] != -1 && diagAnti[idxAnti] != allowedRow) return false;
        return true;
    };

    auto trySwap = [&](int rowA, int rowB, double temperature) {
        if (rowA == rowB) return false;
        int colA = perm[rowA];
        int colB = perm[rowB];
        if (!diagFreeAfterSwap(rowA, colB, rowB)) return false;
        if (!diagFreeAfterSwap(rowB, colA, rowA)) return false;

        long long beforeCost = cost;

        // Remove old diagonals
        diagMain[mainIdx(rowA, colA)] = -1;
        diagAnti[antiIdx(rowA, colA)] = -1;
        diagMain[mainIdx(rowB, colB)] = -1;
        diagAnti[antiIdx(rowB, colB)] = -1;

        std::vector<PairDelta> deltas;
        deltas.reserve(static_cast<std::size_t>(n) * 2);

        auto updatePair = [&](int x, int y) {
            if (x > y) std::swap(x, y);
            std::size_t idx = static_cast<std::size_t>(x) * n + y;
            LineKey oldKey = pairKeys[idx];
            decrementLine(linePairs, oldKey, cost);
            LineKey newKey = makeLineKey(x + 1, perm[x] + 1, y + 1, perm[y] + 1);
            incrementLine(linePairs, newKey, cost);
            pairKeys[idx] = newKey;
            deltas.push_back({idx, oldKey, newKey});
        };

        std::swap(perm[rowA], perm[rowB]);
        colOwner[perm[rowA]] = rowA;
        colOwner[perm[rowB]] = rowB;
        diagMain[mainIdx(rowA, perm[rowA])] = rowA;
        diagAnti[antiIdx(rowA, perm[rowA])] = rowA;
        diagMain[mainIdx(rowB, perm[rowB])] = rowB;
        diagAnti[antiIdx(rowB, perm[rowB])] = rowB;

        for (int k = 0; k < n; ++k) {
            if (k == rowA || k == rowB) continue;
            updatePair(rowA, k);
            updatePair(rowB, k);
        }
        updatePair(rowA, rowB);

        long long afterCost = cost;
        bool accept = (afterCost <= beforeCost);
        if (!accept && temperature > 0.0) {
            double p = std::exp(static_cast<double>(beforeCost - afterCost) /
                                std::max(temperature, 1e-9));
            if (prob(rng) < p) accept = true;
        }

        if (!accept) {
            for (auto it = deltas.rbegin(); it != deltas.rend(); ++it) {
                decrementLine(linePairs, it->newKey, cost);
                incrementLine(linePairs, it->oldKey, cost);
                pairKeys[it->index] = it->oldKey;
            }
            diagMain[mainIdx(rowA, perm[rowA])] = -1;
            diagAnti[antiIdx(rowA, perm[rowA])] = -1;
            diagMain[mainIdx(rowB, perm[rowB])] = -1;
            diagAnti[antiIdx(rowB, perm[rowB])] = -1;
            std::swap(perm[rowA], perm[rowB]);
            colOwner[colA] = rowA;
            colOwner[colB] = rowB;
            diagMain[mainIdx(rowA, colA)] = rowA;
            diagAnti[antiIdx(rowA, colA)] = rowA;
            diagMain[mainIdx(rowB, colB)] = rowB;
            diagAnti[antiIdx(rowB, colB)] = rowB;
            cost = beforeCost;
            return false;
        }
        return true;
    };

    std::vector<int> conflictRows;
    std::vector<char> conflictMarks(n, 0);
    std::vector<char> forbiddenCols(n, 0);
    auto refreshConflictRows = [&]() {
        std::fill(conflictMarks.begin(), conflictMarks.end(), 0);
        conflictRows.clear();
        for (const auto &entry : linePairs) {
            if (entry.second < 3) continue;
            const auto &key = entry.first;
            for (int row = 0; row < n; ++row) {
                if (conflictMarks[row]) continue;
                long long lhs = 1LL * key.a * (row + 1) + 1LL * key.b * (perm[row] + 1) + key.c;
                if (lhs == 0) {
                    conflictMarks[row] = 1;
                    conflictRows.push_back(row);
                }
            }
        }
    };

    refreshConflictRows();

    auto markForbidden = [&](int row) {
        std::fill(forbiddenCols.begin(), forbiddenCols.end(), 0);
        for (const auto &entry : linePairs) {
            if (entry.second < 3) continue;
            const auto &key = entry.first;
            long long numerator = -1LL * key.a * (row + 1) - key.c;
            long long denom = key.b;
            if (numerator % denom != 0) continue;
            long long col = numerator / denom - 1;
            if (col >= 0 && col < n) {
                forbiddenCols[static_cast<std::size_t>(col)] = 1;
            }
        }
        forbiddenCols[perm[row]] = 1;
    };

    auto targetedSwap = [&](int rowA, double temperature) {
        markForbidden(rowA);
        for (int attempt = 0; attempt < n / 2 + 1; ++attempt) {
            int col = rng() % n;
            if (forbiddenCols[col]) continue;
            int rowB = colOwner[col];
            if (rowA == rowB) continue;
            if (trySwap(rowA, rowB, temperature)) {
                return true;
            }
        }
        return false;
    };

    for (long long step = 0; step < maxSteps; ++step) {
        if (cost == 0) return true;
        if (conflictRows.empty()) {
            refreshConflictRows();
            if (conflictRows.empty()) {
                return cost == 0;
            }
        }
        double temperature = 1.0 - static_cast<double>(step) / maxSteps;
        int rowA = conflictRows[rng() % conflictRows.size()];
        bool improved = targetedSwap(rowA, temperature);
        if (!improved) {
            int rowB = rng() % n;
            if (rowA != rowB) {
                improved = trySwap(rowA, rowB, temperature);
            }
        }
        if (improved) {
            refreshConflictRows();
        }
    }
    refreshConflictRows();
    return cost == 0 && conflictRows.empty();
}

void searchWorker(int workerId, int n, const SearchConfig &config,
                  bool skipCollinearity,
                  std::atomic<bool> &solved, std::vector<int> &solution,
                  std::mutex &solutionMutex, std::atomic<std::uint64_t> &attempts) {
    std::mt19937_64 rng(std::random_device{}() ^ (static_cast<std::uint64_t>(workerId) << 32));
    while (!solved.load(std::memory_order_relaxed)) {
        std::vector<int> perm;
        if (!solveQueens(n, rng, config, perm)) {
            continue;
        }
        if (!skipCollinearity) {
            if (!enforceCollinearity(perm, rng)) {
                continue;
            }
        }
        attempts.fetch_add(1, std::memory_order_relaxed);
        if (!skipCollinearity) {
            auto result = validatePermutation(perm);
            if (!result.ok) {
                continue;
            }
        }
        bool expected = false;
        if (solved.compare_exchange_strong(expected, true)) {
            std::lock_guard<std::mutex> guard(solutionMutex);
            solution = std::move(perm);
        }
        return;
    }
}

int main(int argc, char *argv[]) {
    int N = 999;
    bool skipCollinearity = false;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--skip-collinear" || arg == "--skip-collinearity") {
            skipCollinearity = true;
        } else {
            try {
                N = std::stoi(arg);
            } catch (...) {
                std::cerr << "Ignoring invalid argument: " << arg << '\n';
            }
        }
    }
    if (skipCollinearity) {
        std::cerr << "[warn] Skipping collinearity validation; use only for diagnostics." << std::endl;
    }
    if (N < 4) {
        std::cerr << "N must be at least 4." << std::endl;
        return 1;
    }

    SearchConfig config;
    config.maxRestarts = 64;
    config.maxIterations = N * 12;

    std::atomic<bool> solved{false};
    std::atomic<std::uint64_t> attempts{0};
    std::vector<int> solution;
    std::mutex solutionMutex;

    const unsigned int hw = std::thread::hardware_concurrency();
    const unsigned int numThreads = std::max(2u, hw == 0 ? 2u : hw);

    std::vector<std::thread> workers;
    workers.reserve(numThreads);

    auto start = std::chrono::steady_clock::now();
    for (unsigned int i = 0; i < numThreads; ++i) {
        workers.emplace_back(searchWorker, static_cast<int>(i), N, std::cref(config),
                             skipCollinearity,
                             std::ref(solved), std::ref(solution), std::ref(solutionMutex), std::ref(attempts));
    }
    for (auto &t : workers) {
        t.join();
    }
    auto stop = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

    if (!solved.load()) {
        std::cerr << "Failed to find a valid configuration for N=" << N << std::endl;
        return 2;
    }

    std::cerr << "Solved N=" << N << " in " << elapsed << " ms after "
              << attempts.load() << " validated permutations." << std::endl;

    std::cout << N << '\n';
    for (int row = 0; row < N; ++row) {
        if (row) std::cout << ' ';
        std::cout << solution[row] + 1;
    }
    std::cout << '\n';
    return 0;
}
