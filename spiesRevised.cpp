#include <iostream>
#include <vector>
#include <thread>
#include <atomic>
#include <mutex>
#include <random>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <limits>

using namespace std;

struct SpinBarrier {
    int total;
    atomic<int> arrived{0};
    atomic<bool> go{false};

    explicit SpinBarrier(int t) : total(t) {}

    void wait() {
        int a = arrived.fetch_add(1, memory_order_acq_rel) + 1;
        if (a == total) go.store(true, memory_order_release);
        while (!go.load(memory_order_acquire)) {
            this_thread::yield();
        }
    }
};

class MinConflictSolver {
public:
    MinConflictSolver(int n, uint64_t seed)
        : n_(n),
          UNASSIGNED_(2 * n - 1),
          initialized_(false),
          conflicts_(n, vector<int>(n, 0)),          // conflicts_[col][row]
          simpleInv_(n),
          colInv_(n),
          rowOfCol_(n, UNASSIGNED_),
          rng_(seed32(seed)) {

        for (int c = 0; c < n_; ++c) simpleInv_[c].reserve(3 * n_);
        for (int i = 0; i < n_; ++i) {
            int len = n_ - i - 1;
            if (len > 0) colInv_[i].resize(len);
        }
        greedyInitialize();
    }

    void greedyInitialize() {
        for (int c = 0; c < n_; ++c) {
            fill(conflicts_[c].begin(), conflicts_[c].end(), 0);
            simpleInv_[c].clear();
        }
        for (int i = 0; i < n_; ++i) {
            for (auto &v : colInv_[i]) v.clear();
        }
        fill(rowOfCol_.begin(), rowOfCol_.end(), UNASSIGNED_);
        initialized_ = false;

        vector<int> bestRows;
        bestRows.reserve(n_);
        for (int col = 0; col < n_; ++col) {
            revalidate(col);

            int best = numeric_limits<int>::max();
            bestRows.clear();
            for (int row = 0; row < n_; ++row) {
                int c = conflicts_[col][row];
                if (c < best) {
                    best = c;
                    bestRows.clear();
                    bestRows.push_back(row);
                } else if (c == best) {
                    bestRows.push_back(row);
                }
            }
            int chosen = bestRows[randIndex((int)bestRows.size())];
            rowOfCol_[col] = chosen;
            invalidate(col);
        }
        initialized_ = true;
    }

    bool solved() const {
        for (int col = 0; col < n_; ++col) {
            int row = rowOfCol_[col];
            if (row < 0 || row >= n_) return false;
            if (conflicts_[col][row] != 0) return false;
        }
        return true;
    }

    bool solve(int maxSteps, atomic<bool>* stopFlag, uint64_t* stepsDone) {
        vector<int> worstCols;
        worstCols.reserve(n_);
        vector<int> candidates;
        candidates.reserve(n_);
        vector<int> secondCandidates;
        secondCandidates.reserve(n_);

        for (int step = 0; step < maxSteps; ++step) {
            if (stepsDone) ++(*stepsDone);
            if (stopFlag && stopFlag->load(memory_order_relaxed)) return false;

            int maxConf = 0;
            worstCols.clear();
            for (int col = 0; col < n_; ++col) {
                int c = conflicts_[col][rowOfCol_[col]];
                if (c > maxConf) {
                    maxConf = c;
                    worstCols.clear();
                    if (c > 0) worstCols.push_back(col);
                } else if (c == maxConf && c > 0) {
                    worstCols.push_back(col);
                }
            }
            if (maxConf == 0) return true;

            int col = worstCols[randIndex((int)worstCols.size())];
            int curRow = rowOfCol_[col];

            int minC = numeric_limits<int>::max();
            for (int row = 0; row < n_; ++row) minC = min(minC, conflicts_[col][row]);

            candidates.clear();
            for (int row = 0; row < n_; ++row) {
                if (conflicts_[col][row] == minC) candidates.push_back(row);
            }

            int newRow = curRow;
            if ((int)candidates.size() == 1 && candidates[0] == curRow && minC > 0) {
                int second = numeric_limits<int>::max();
                for (int row = 0; row < n_; ++row) {
                    if (row == curRow) continue;
                    second = min(second, conflicts_[col][row]);
                }
                secondCandidates.clear();
                for (int row = 0; row < n_; ++row) {
                    if (row == curRow) continue;
                    if (conflicts_[col][row] == second) secondCandidates.push_back(row);
                }
                if (!secondCandidates.empty()) {
                    newRow = secondCandidates[randIndex((int)secondCandidates.size())];
                } else {
                    newRow = (curRow + 1 + randIndex(n_ - 1)) % n_;
                }
            } else {
                newRow = candidates[randIndex((int)candidates.size())];
            }

            if (newRow != curRow) {
                revalidate(col);
                rowOfCol_[col] = newRow;
                invalidate(col);
            }
        }
        return solved();
    }

    const vector<int>& rowOfCol() const { return rowOfCol_; }

private:
    int n_;
    int UNASSIGNED_;
    bool initialized_;

    vector<vector<int>> conflicts_;
    vector<vector<uint32_t>> simpleInv_;
    vector<vector<vector<uint32_t>>> colInv_;
    vector<int> rowOfCol_;
    mt19937 rng_;

    static uint32_t seed32(uint64_t s) {
        s ^= (s >> 33);
        s *= 0xff51afd7ed558ccdULL;
        s ^= (s >> 33);
        s *= 0xc4ceb9fe1a85ec53ULL;
        s ^= (s >> 33);
        return (uint32_t)(s & 0xFFFFFFFFu);
    }

    inline int randIndex(int m) { return (int)(rng_() % (uint32_t)m); }

    static inline uint32_t enc(int row, int col) { return (uint32_t(row) << 16) | uint32_t(col); }
    static inline int decRow(uint32_t v) { return int(v >> 16); }
    static inline int decCol(uint32_t v) { return int(v & 0xFFFFu); }

    inline void inc(int row, int col) { ++conflicts_[col][row]; }
    inline void dec(int row, int col) { --conflicts_[col][row]; }

    void invalidateSpotSimple(int row, int col, int ownerCol) {
        inc(row, col);
        simpleInv_[ownerCol].push_back(enc(row, col));
    }

    void invalidateSimple(int col) {
        int row = rowOfCol_[col];
        for (int c = 0; c < n_; ++c) {
            if (c == col) continue;
            invalidateSpotSimple(row, c, col);
        }
        for (int d = 1; row + d < n_ && col + d < n_; ++d) invalidateSpotSimple(row + d, col + d, col);
        for (int d = 1; row - d >= 0 && col + d < n_; ++d) invalidateSpotSimple(row - d, col + d, col);
        for (int d = 1; row + d < n_ && col - d >= 0; ++d) invalidateSpotSimple(row + d, col - d, col);
        for (int d = 1; row - d >= 0 && col - d >= 0; ++d) invalidateSpotSimple(row - d, col - d, col);
    }

    void revalidateSimple(int col) {
        for (uint32_t idx : simpleInv_[col]) {
            int r = decRow(idx);
            int c = decCol(idx);
            dec(r, c);
        }
        simpleInv_[col].clear();
    }

    void invalidateColinearPair(int col1, int col2) {
        int row1 = rowOfCol_[col1];
        int row2 = rowOfCol_[col2];
        auto &vec = colInv_[col1][col2 - col1 - 1];

        if (row1 == UNASSIGNED_ && row2 == UNASSIGNED_) return;
        if (row1 == UNASSIGNED_ && row2 != UNASSIGNED_) {
            inc(row2, col2);
            vec.push_back(enc(row2, col2));
            return;
        }
        if (row2 == UNASSIGNED_ && row1 != UNASSIGNED_) {
            inc(row1, col1);
            vec.push_back(enc(row1, col1));
            return;
        }

        int dr = row2 - row1;
        int dc = col2 - col1;
        int g = std::gcd(std::abs(dr), dc);
        dr /= g;
        dc /= g;

        int r = row1, c = col1;
        while (0 <= r && r < n_ && 0 <= c && c < n_) {
            inc(r, c);
            vec.push_back(enc(r, c));
            r += dr;
            c += dc;
        }
        r = row1 - dr;
        c = col1 - dc;
        while (0 <= r && r < n_ && 0 <= c && c < n_) {
            inc(r, c);
            vec.push_back(enc(r, c));
            r -= dr;
            c -= dc;
        }
    }

    void invalidateColinear(int col) {
        for (int i = 0; i < col; ++i) invalidateColinearPair(i, col);
        for (int j = col + 1; j < n_; ++j) invalidateColinearPair(col, j);
        conflicts_[col][rowOfCol_[col]] -= (n_ - 1);
    }

    void revalidateColinear(int col) {
        for (int i = 0; i < col; ++i) {
            auto &v = colInv_[i][col - i - 1];
            for (uint32_t idx : v) {
                int r = decRow(idx);
                int c = decCol(idx);
                dec(r, c);
            }
            v.clear();
        }
        for (int j = col + 1; j < n_; ++j) {
            auto &v = colInv_[col][j - col - 1];
            for (uint32_t idx : v) {
                int r = decRow(idx);
                int c = decCol(idx);
                dec(r, c);
            }
            v.clear();
        }
        if (initialized_) {
            conflicts_[col][rowOfCol_[col]] += (n_ - 1);
        }
    }

    void invalidate(int col) {
        invalidateSimple(col);
        invalidateColinear(col);
    }

    void revalidate(int col) {
        revalidateSimple(col);
        revalidateColinear(col);
    }
};

int main(int argc, char** argv) {
    int N = 999;
    if (argc >= 2) {
        try {
            N = max(4, stoi(argv[1]));
        } catch (...) {
            N = 999;
        }
    }

    unsigned hw = std::thread::hardware_concurrency();
    int threads = (hw ? (int)hw : 4);
    threads = std::min(threads, 8);
    if (argc >= 3) {
        try {
            threads = max(1, stoi(argv[2]));
        } catch (...) {
        }
    }

    int maxSteps = 200000;
    if (argc >= 4) {
        try {
            maxSteps = max(1000, stoi(argv[3]));
        } catch (...) {
        }
    }
    bool printStats = (argc >= 5);

    atomic<bool> found(false);
    mutex resultMutex;
    vector<int> bestRowOfCol;

    SpinBarrier barrier(threads);
    vector<atomic<uint64_t>> steps(threads);
    vector<atomic<uint64_t>> restarts(threads);
    for (int i = 0; i < threads; ++i) {
        steps[i].store(0);
        restarts[i].store(0);
    }

    auto worker = [&](int id) {
        uint64_t seed =
            (uint64_t)chrono::high_resolution_clock::now().time_since_epoch().count()
            ^ (0x9e3779b97f4a7c15ULL * (uint64_t)(id + 1))
            ^ (((uint64_t)random_device{}() << 32) | (uint64_t)random_device{}());

        MinConflictSolver solver(N, seed);
        barrier.wait();

        uint64_t localSteps = 0;
        uint64_t localRestarts = 0;

        while (!found.load(memory_order_relaxed)) {
            if (solver.solve(maxSteps, &found, &localSteps)) {
                bool expected = false;
                if (found.compare_exchange_strong(expected, true, memory_order_acq_rel)) {
                    lock_guard<mutex> lk(resultMutex);
                    bestRowOfCol = solver.rowOfCol();
                }
                break;
            }
            ++localRestarts;
            solver.greedyInitialize();
        }

        steps[id].store(localSteps, memory_order_relaxed);
        restarts[id].store(localRestarts, memory_order_relaxed);
    };

    vector<thread> pool;
    pool.reserve((size_t)threads);
    for (int t = 0; t < threads; ++t) pool.emplace_back(worker, t);
    for (auto &th : pool) th.join();

    if (bestRowOfCol.empty()) {
        cerr << "Failed to find a solution for N=" << N << "\n";
        return 1;
    }

    if (printStats) {
        for (int i = 0; i < threads; ++i) {
            cerr << "thread " << i
                 << " steps=" << steps[i].load()
                 << " restarts=" << restarts[i].load() << "\n";
        }
    }

    cout << N << "\n";
    vector<int> rowToCol(N, -1);
    for (int col = 0; col < N; ++col) {
        int row = bestRowOfCol[col];
        if (row >= 0 && row < N) rowToCol[row] = col;
    }
    for (int row = 0; row < N; ++row) {
        if (row) cout << ' ';
        cout << rowToCol[row] + 1;
    }
    cout << "\n";
    return 0;
}
