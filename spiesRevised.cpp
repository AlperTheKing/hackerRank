#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <utility>
#include <cstdlib>
#include <chrono>

class Solver {
public:
    explicit Solver(int n) : n(n) {
        std::srand(static_cast<unsigned int>(std::time(nullptr)));
        numConflicts.assign(n, std::vector<int>(n, 0));
        simpleInvalidations.assign(n, {});
        colinearInvalidations.assign(n, std::vector<std::vector<std::pair<int, int>>>(n));
        placementOptions.assign(n, {});
        rows.assign(n, n * 2 - 1);

        for (int col = 0; col < n; col++) {
            rows[col] = std::rand() % n;
            revalidate(col);
            invalidate(col);
        }
        initialized = true;
    }

    bool solve(int maxIter = -1) {
        int i = 0;
        while (!solved()) {
            if (maxIter >= 0 && i == maxIter) return false;

            calculatePlacementOptions();
            if (changableColumns.empty()) return false;

            int col = getChangableCol();
            revalidate(col);
            int row = getPlacement(col);
            rows[col] = row;
            invalidate(col);
            i++;
        }
        return true;
    }

    std::vector<int> getSolution() const { return rows; }

private:
    int n;
    std::vector<std::vector<int>> numConflicts;
    std::vector<std::vector<std::pair<int, int>>> simpleInvalidations;
    std::vector<std::vector<std::vector<std::pair<int, int>>>> colinearInvalidations;
    std::vector<int> rows;
    std::vector<std::vector<int>> placementOptions;
    std::vector<int> changableColumns;
    bool initialized = false;

    void invalidateSpot(int row, int col, int originalColumn) {
        numConflicts[row][col] += 1;
        simpleInvalidations[originalColumn].push_back(std::make_pair(row, col));
    }

    void invalidateSimple(int col) {
        int row = rows[col];
        for (int i = 0; i < col; i++) invalidateSpot(row, i, col);
        for (int i = col + 1; i < n; i++) invalidateSpot(row, i, col);
        for (int i = 1; row + i < n && col + i < n; i++) invalidateSpot(row + i, col + i, col);
        for (int i = 1; 0 <= row - i && col + i < n; i++) invalidateSpot(row - i, col + i, col);
        for (int i = 1; row + i < n && 0 <= col - i; i++) invalidateSpot(row + i, col - i, col);
        for (int i = 1; 0 <= row - i && 0 <= col - i; i++) invalidateSpot(row - i, col - i, col);
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

        int row = row1;
        int col = col1;
        if (rowStep > 0) {
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
        for (int i = 0; i < col; i++) invalidateColinearPair(i, col);
        for (int i = col + 1; i < n; i++) invalidateColinearPair(col, i);
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
        if (initialized) numConflicts[rows[col]][col] += n - 1;
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
        return changableColumns[std::rand() % changableColumns.size()];
    }

    int getPlacement(int col) {
        return placementOptions[col][std::rand() % placementOptions[col].size()];
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

int main(int argc, char *argv[]) {
    unsigned seed = static_cast<unsigned>(
        std::chrono::steady_clock::now().time_since_epoch().count());
    std::srand(seed);

    int N = 999;
    if (argc > 1) {
        try {
            N = std::stoi(argv[1]);
        } catch (...) {
            N = 999;
        }
    }
    if (N < 1) N = 999;

    Solver solver(N);
    solver.solve();
    auto sol = solver.getSolution();

    std::cout << N << '\n';
    for (int i = 0; i < N; i++) {
        if (i) std::cout << ' ';
        std::cout << sol[i] + 1;
    }
    std::cout << '\n';
    return 0;
}
