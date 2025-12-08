#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <unistd.h>
#include <vector>

using namespace std;

// --- Config ---
// Time limit lowered slightly for safety

struct City {
  int x, y;
  int p;
  double distHQ;
};

// --- Globals ---
int N;
double C_cost, D_decay;
vector<City> cities;
vector<char> visited;
int sold_count = 0;

// --- Dynamic Grid ---
int MIN_X, MAX_X, MIN_Y, MAX_Y;
int CELLS_X, CELLS_Y;
const int CELL_SIZE = 20;

// Flat grid logic
vector<vector<int>> grid_vec;

// --- State ---
bool query_printed = false;

// --- Helper Functions ---
static inline double euclid(double x1, double y1, double x2, double y2) {
  double dx = x1 - x2;
  double dy = y1 - y2;
  return sqrt(dx * dx + dy * dy);
}

struct Trip {
  double profit = -1e100;
  vector<int> order;
  Trip() {}
  Trip(double p, vector<int> o) : profit(p), order(o) {}
};

// --- Price Calc ---
vector<double> stageMul(10, 1.0);
void init_stages() {
  for (int i = 1; i < 10; ++i)
    stageMul[i] = stageMul[i - 1] * D_decay;
}
double get_sale_factor(int sale_idx) {
  int B = max(1, N / 10);
  int stage = sale_idx / B;
  if (stage > 9)
    stage = 9;
  return stageMul[stage];
}

// Logic
vector<int> buildLcand(int remaining) {
  // Aggressive pruning for speed
  int maxL;
  if (C_cost <= 0.5)
    maxL = 6;
  else if (C_cost <= 1.5)
    maxL = 5;
  else
    maxL = 4;

  maxL = min(maxL, remaining);

  vector<int> base = {1, 2, 3, 4, 5, 6};
  vector<int> res;
  for (int L : base)
    if (L <= maxL)
      res.push_back(L);
  if (res.empty() && remaining > 0)
    res.push_back(1);
  return res;
}

double tripProfit(const vector<int> &order, bool returnToHQ, int startSold) {
  int L = (int)order.size();
  double profit = 0.0;
  double cx = 0, cy = 0;
  int remaining = L;
  for (int i = 0; i < L; ++i) {
    int id = order[i];
    double d = euclid(cx, cy, cities[id].x, cities[id].y);
    profit -= d * (1.0 + remaining * C_cost);

    double mult = get_sale_factor(startSold + i);
    profit += cities[id].p * mult;

    cx = cities[id].x;
    cy = cities[id].y;
    remaining--;
  }
  if (returnToHQ && L > 0) {
    profit -= euclid(cx, cy, 0, 0);
  }
  return profit;
}

Trip optimizeTrip(vector<int> order, bool returnToHQ, int startSold) {
  double bestP = tripProfit(order, returnToHQ, startSold);
  bool improved = true;
  int L = (int)order.size();

  for (int pass = 0; pass < 2 && improved; ++pass) {
    improved = false;
    for (int i = 0; i < L; ++i) {
      for (int j = i + 1; j < L; ++j) {
        swap(order[i], order[j]);
        double p = tripProfit(order, returnToHQ, startSold);
        if (p > bestP + 1e-9) {
          bestP = p;
          improved = true;
        } else {
          swap(order[i], order[j]);
        }
      }
    }
  }
  return {bestP, order};
}

vector<int> candMark;
int candToken = 1;
vector<int> usedMark;
int usedToken = 1;

void bumpToken(int &token, vector<int> &arr) {
  token++;
  if (token >= 2000000000) {
    fill(arr.begin(), arr.end(), 0);
    token = 1;
  }
}

// Global Sort Caches
vector<int> byPrice;
vector<int> byGain;
int dirty_count = 0;

void compact_vectors() {
  vector<int> newPrice, newGain;
  newPrice.reserve(N - sold_count);
  newGain.reserve(N - sold_count);
  for (int id : byPrice)
    if (!visited[id])
      newPrice.push_back(id);
  for (int id : byGain)
    if (!visited[id])
      newGain.push_back(id);
  byPrice = move(newPrice);
  byGain = move(newGain);
  dirty_count = 0;
}

clock_t start_time_global;

Trip planTrip(bool returnToHQ) {
  Trip best;
  int remaining = N - sold_count;
  if (remaining <= 0)
    return best;

  if (dirty_count > 1000 && dirty_count > remaining * 0.5) {
    compact_vectors();
  }

  vector<int> Lcand = buildLcand(remaining);

  int topPriceP = 25;
  int topGainP = 25;
  int nearLimit = 50;

  for (int L : Lcand) {

    vector<int> order;
    order.reserve(L);
    double cx = 0, cy = 0;
    double curDistHQ = 0;

    bumpToken(usedToken, usedMark);

    bool ok = true;
    for (int step = 0; step < L; ++step) {
      int k_carry = L - step;
      double mult = get_sale_factor(sold_count + step);

      bumpToken(candToken, candMark);
      vector<int> cands;
      cands.reserve(topPriceP + topGainP + nearLimit);

      // Add Top Price
      int added = 0;
      for (int id : byPrice) {
        if (visited[id] || candMark[id] == candToken)
          continue;
        candMark[id] = candToken;
        cands.push_back(id);
        if (++added >= topPriceP)
          break;
      }

      // Add Top Gain
      added = 0;
      for (int id : byGain) {
        if (visited[id] || candMark[id] == candToken)
          continue;
        candMark[id] = candToken;
        cands.push_back(id);
        if (++added >= topGainP)
          break;
      }

      // Add Near
      if (CELLS_X > 0 && CELLS_Y > 0) {
        int gx = (int)((cx - MIN_X) / CELL_SIZE);
        int gy = (int)((cy - MIN_Y) / CELL_SIZE);
        gx = max(0, min(CELLS_X - 1, gx));
        gy = max(0, min(CELLS_Y - 1, gy));

        int rad = 0;
        int found_near = 0;
        while (found_near < nearLimit && rad < 40) {
          int r_min = max(0, gy - rad);
          int r_max = min(CELLS_Y - 1, gy + rad);
          int c_min = max(0, gx - rad);
          int c_max = min(CELLS_X - 1, gx + rad);

          for (int r = r_min; r <= r_max; ++r) {
            for (int c = c_min; c <= c_max; ++c) {
              if (rad > 0 && r > r_min && r < r_max && c > c_min && c < c_max)
                continue;
              int g_idx = c * CELLS_Y + r;
              if (g_idx < 0 || g_idx >= (int)grid_vec.size())
                continue;

              for (int id : grid_vec[g_idx]) {
                if (visited[id] || candMark[id] == candToken)
                  continue;
                candMark[id] = candToken;
                cands.push_back(id);
                found_near++;
                if (found_near >= nearLimit)
                  goto end_near;
              }
            }
          }
          rad++;
        }
      end_near:;
      }

      int bestId = -1;
      double bestVal = -1e100;

      for (int id : cands) {
        if (usedMark[id] == usedToken)
          continue;

        double d = euclid(cx, cy, cities[id].x, cities[id].y);
        double gain = mult * cities[id].p - d * (1.0 + k_carry * C_cost);
        if (returnToHQ)
          gain -= (cities[id].distHQ - curDistHQ);

        if (gain > bestVal) {
          bestVal = gain;
          bestId = id;
        }
      }

      if (bestId == -1) {
        ok = false;
        break;
      }
      order.push_back(bestId);
      usedMark[bestId] = usedToken;
      cx = cities[bestId].x;
      cy = cities[bestId].y;
      curDistHQ = cities[bestId].distHQ;
    }

    if (!ok)
      continue;
    Trip optimized = optimizeTrip(order, returnToHQ, sold_count);
    if (optimized.profit > best.profit)
      best = optimized;
  }
  return best;
}

void solve() {
  start_time_global = clock();
  cities.clear();
  cities.resize(N);
  visited.assign(N, 0);
  sold_count = 0;
  dirty_count = 0;
  query_printed = false;

  // 1. Read Cities & Find Bounds
  MIN_X = 1e9;
  MAX_X = -1e9;
  MIN_Y = 1e9;
  MAX_Y = -1e9;

  for (int i = 0; i < N; ++i) {
    if (scanf("%d %d %d", &cities[i].x, &cities[i].y, &cities[i].p) != 3) {
      break;
    }
    cities[i].distHQ = euclid(cities[i].x, cities[i].y, 0, 0);
    if (cities[i].x < MIN_X)
      MIN_X = cities[i].x;
    if (cities[i].x > MAX_X)
      MAX_X = cities[i].x;
    if (cities[i].y < MIN_Y)
      MIN_Y = cities[i].y;
    if (cities[i].y > MAX_Y)
      MAX_Y = cities[i].y;
  }

  // 2. Setup Dynamic Grid
  MIN_X -= 100;
  MIN_Y -= 100;

  CELLS_X = (MAX_X - MIN_X) / CELL_SIZE + 5;
  CELLS_Y = (MAX_Y - MIN_Y) / CELL_SIZE + 5;

  // Safety cap
  if ((long long)CELLS_X * CELLS_Y > 500000) {
    CELLS_X = 0;
    CELLS_Y = 0;
    grid_vec.clear();
  } else {
    grid_vec.assign(CELLS_X * CELLS_Y, vector<int>());
    for (int i = 0; i < N; ++i) {
      int gx = (cities[i].x - MIN_X) / CELL_SIZE;
      int gy = (cities[i].y - MIN_Y) / CELL_SIZE;
      if (gx >= 0 && gx < CELLS_X && gy >= 0 && gy < CELLS_Y) {
        grid_vec[gx * CELLS_Y + gy].push_back(i);
      }
    }
  }

  init_stages();
  candMark.assign(N, 0);
  candToken = 1;
  usedMark.assign(N, 0);
  usedToken = 1;
  byPrice.resize(N);
  iota(byPrice.begin(), byPrice.end(), 0);
  byGain.resize(N);
  iota(byGain.begin(), byGain.end(), 0);

  sort(byPrice.begin(), byPrice.end(), [&](int a, int b) {
    if (cities[a].p != cities[b].p)
      return cities[a].p > cities[b].p;
    return cities[a].distHQ < cities[b].distHQ;
  });
  sort(byGain.begin(), byGain.end(), [&](int a, int b) {
    double ga = cities[a].p - (2.0 + C_cost) * cities[a].distHQ;
    double gb = cities[b].p - (2.0 + C_cost) * cities[b].distHQ;
    return ga > gb;
  });

  while (true) {

    Trip t = planTrip(true);
    if (t.order.empty() || t.profit <= 1e-9)
      break;

    for (size_t i = 0; i < t.order.size(); ++i) {
      int id = t.order[i];
      visited[id] = 1;
      sold_count++;
      if (i == 0)
        printf("%d %d %d\n", cities[id].x, cities[id].y, (int)t.order.size());
      else
        printf("%d %d\n", cities[id].x, cities[id].y);
    }
    printf("0 0\n");
    fflush(stdout);
    query_printed = true;
    dirty_count += t.order.size();
  }

  Trip last = planTrip(false);
  if (!last.order.empty() && last.profit > 1e-9) {
    for (size_t i = 0; i < last.order.size(); ++i) {
      int id = last.order[i];
      visited[id] = 1;
      sold_count++;
      if (i == 0)
        printf("%d %d %d\n", cities[id].x, cities[id].y,
               (int)last.order.size());
      else
        printf("%d %d\n", cities[id].x, cities[id].y);
    }
    fflush(stdout);
    query_printed = true;
  }

  if (!query_printed && sold_count == 0 && N > 0) {
    int bestId = 0;
    double bestScore = -1e100;
    for (int i = 0; i < N; ++i) {
      double d = cities[i].distHQ;
      double score = (double)cities[i].p - d * (1.0 + C_cost);
      if (score > bestScore) {
        bestScore = score;
        bestId = i;
      }
    }
    printf("%d %d 1\n", cities[bestId].x, cities[bestId].y);
    fflush(stdout);
  }
}

int main() {
  setvbuf(stdout, NULL, _IONBF, 0);
  while (scanf("%d %lf %lf", &N, &C_cost, &D_decay) == 3) {
    solve();
  }
  return 0;
}
