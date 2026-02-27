#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <numeric>
#include <thread>
#include <vector>

using namespace std;

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
vector<double> saleMul;
void init_stages() {
  int B = max(1, N / 10);
  vector<double> stageMul(10, 1.0);
  for (int i = 1; i < 10; ++i) stageMul[i] = stageMul[i - 1] * D_decay;
  saleMul.assign(N + 8, stageMul[9]);
  for (int i = 0; i < N + 8; ++i) {
    int stage = i / B;
    if (stage > 9) stage = 9;
    saleMul[i] = stageMul[stage];
  }
}

// Logic
vector<int> buildLcand(int remaining) {
  int maxL;
  vector<int> base;
  if (N <= 500) {
    if (C_cost <= 0.5)
      maxL = 6;
    else if (C_cost <= 1.5)
      maxL = 5;
    else
      maxL = 4;
    maxL = min(maxL, remaining);
    base = {1, 2, 3, 4, 5, 6};
  } else {
    if (C_cost <= 0.6)
      maxL = 10;
    else if (C_cost <= 1.0)
      maxL = 9;
    else if (C_cost <= 1.8)
      maxL = 8;
    else if (C_cost <= 2.8)
      maxL = 7;
    else
      maxL = 6;
    maxL = min(maxL, remaining);
    for (int L = maxL; L >= 1; --L) base.push_back(L);
  }

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

    double mult = saleMul[startSold + i];
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

  int passLimit = (N <= 20 ? 2 : 3);
  for (int pass = 0; pass < passLimit && improved; ++pass) {
    improved = false;

    // Swap neighborhoods.
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

    if (N > 20) {
      // Insertion neighborhoods.
      for (int i = 0; i < L; ++i) {
        int keep = order[i];
        for (int j = 0; j < L; ++j) {
          if (i == j) continue;
          vector<int> cand = order;
          cand.erase(cand.begin() + i);
          cand.insert(cand.begin() + j, keep);
          double p = tripProfit(cand, returnToHQ, startSold);
          if (p > bestP + 1e-9) {
            bestP = p;
            order = std::move(cand);
            improved = true;
          }
        }
      }
    }
  }
  return {bestP, order};
}

vector<int> candMark;
int candToken = 1;

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
  byPrice = std::move(newPrice);
  byGain = std::move(newGain);
  dirty_count = 0;
}

static constexpr int MAX_ROUTE = 12;
struct BeamState {
  double profit = 0.0;
  double cx = 0.0, cy = 0.0;
  double distHQ = 0.0;
  int len = 0;
  std::array<int, MAX_ROUTE> path{};
};

Trip planTripGreedy(bool returnToHQ) {
  Trip best;
  int remaining = N - sold_count;
  if (remaining <= 0) return best;

  vector<int> Lcand = buildLcand(remaining);

  int topPriceP = 35;
  int topGainP = 35;
  int nearLimit = 70;
  int topPool = 180;
  if (N <= 500) {
    topPriceP = 25;
    topGainP = 25;
    nearLimit = 50;
    topPool = remaining;
  }

  vector<int> topPriceList;
  vector<int> topGainList;
  topPriceList.reserve(topPool);
  topGainList.reserve(topPool);
  for (int id : byPrice) {
    if (!visited[id]) topPriceList.push_back(id);
    if ((int)topPriceList.size() >= topPool) break;
  }
  for (int id : byGain) {
    if (!visited[id]) topGainList.push_back(id);
    if ((int)topGainList.size() >= topPool) break;
  }

  vector<int> cands;
  cands.reserve(topPriceP + topGainP + nearLimit + 16);

  for (int L : Lcand) {
    vector<int> order;
    order.reserve(L);
    double cx = 0, cy = 0;
    double curDistHQ = 0;
    bool ok = true;

    for (int step = 0; step < L; ++step) {
      int k_carry = L - step;
      double mult = saleMul[sold_count + step];
      bumpToken(candToken, candMark);
      cands.clear();

      int added = 0;
      for (int id : topPriceList) {
        if (visited[id] || candMark[id] == candToken) continue;
        candMark[id] = candToken;
        cands.push_back(id);
        if (++added >= topPriceP) break;
      }
      added = 0;
      for (int id : topGainList) {
        if (visited[id] || candMark[id] == candToken) continue;
        candMark[id] = candToken;
        cands.push_back(id);
        if (++added >= topGainP) break;
      }

      if (CELLS_X > 0 && CELLS_Y > 0) {
        int gx = (int)((cx - MIN_X) / CELL_SIZE);
        int gy = (int)((cy - MIN_Y) / CELL_SIZE);
        gx = max(0, min(CELLS_X - 1, gx));
        gy = max(0, min(CELLS_Y - 1, gy));
        int rad = 0;
        int found_near = 0;
        while (found_near < nearLimit && rad < 45) {
          int r_min = max(0, gy - rad);
          int r_max = min(CELLS_Y - 1, gy + rad);
          int c_min = max(0, gx - rad);
          int c_max = min(CELLS_X - 1, gx + rad);
          for (int r = r_min; r <= r_max; ++r) {
            for (int c = c_min; c <= c_max; ++c) {
              if (rad > 0 && r > r_min && r < r_max && c > c_min && c < c_max) continue;
              int g_idx = c * CELLS_Y + r;
              if (g_idx < 0 || g_idx >= (int)grid_vec.size()) continue;
              for (int id : grid_vec[g_idx]) {
                if (visited[id] || candMark[id] == candToken) continue;
                candMark[id] = candToken;
                cands.push_back(id);
                if (++found_near >= nearLimit) goto end_near_g;
              }
            }
          }
          rad++;
        }
      end_near_g:;
      }

      int bestId = -1;
      double bestVal = -1e100;
      for (int id : cands) {
        bool inOrder = false;
        for (int x : order) {
          if (x == id) {
            inOrder = true;
            break;
          }
        }
        if (inOrder) continue;

        double d = euclid(cx, cy, cities[id].x, cities[id].y);
        double gain = mult * cities[id].p - d * (1.0 + k_carry * C_cost);
        if (returnToHQ) gain -= (cities[id].distHQ - curDistHQ);
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
      cx = cities[bestId].x;
      cy = cities[bestId].y;
      curDistHQ = cities[bestId].distHQ;
    }

    if (!ok) continue;
    Trip opt = optimizeTrip(order, returnToHQ, sold_count);
    if (opt.profit > best.profit) best = std::move(opt);
  }
  return best;
}

Trip planTrip(bool returnToHQ) {
  Trip best;
  int remaining = N - sold_count;
  if (remaining <= 0)
    return best;
  if (N <= 20) {
    return planTripGreedy(returnToHQ);
  }
  if (N > 3000) {
    int beamSoldLimit = N / 9;
    if (N >= 30000) beamSoldLimit = N / 12;
    if (sold_count >= beamSoldLimit) {
      return planTripGreedy(returnToHQ);
    }
  }

  if (dirty_count > 1000 && dirty_count > remaining * 0.5) {
    compact_vectors();
  }

  vector<int> Lcand = buildLcand(remaining);

  int topPriceP = 60;
  int topGainP = 60;
  int nearLimit = 120;
  int topPool = 300;
  int beamWidth = (N >= 20000 ? 18 : 28);
  int branchPerState = 12;
  if (N >= 30000) {
    topPriceP = 35;
    topGainP = 35;
    nearLimit = 70;
    topPool = 180;
    beamWidth = 8;
    branchPerState = 6;
  } else if (N >= 10000) {
    topPriceP = 45;
    topGainP = 45;
    nearLimit = 90;
    topPool = 220;
    beamWidth = 12;
    branchPerState = 8;
  }
  if (N <= 2000) {
    topPool = min(1000, remaining);
    topPriceP = min(160, topPool);
    topGainP = min(160, topPool);
    nearLimit = 220;
    beamWidth = 64;
    branchPerState = 24;
  } else if (N <= 10000) {
    topPool = min(700, remaining);
    topPriceP = min(120, topPool);
    topGainP = min(120, topPool);
    nearLimit = 160;
    beamWidth = 40;
    branchPerState = 16;
  }

  vector<int> topPriceList;
  vector<int> topGainList;
  topPriceList.reserve(topPool);
  topGainList.reserve(topPool);

  for (int id : byPrice) {
    if (!visited[id]) topPriceList.push_back(id);
    if ((int)topPriceList.size() >= topPool) break;
  }
  for (int id : byGain) {
    if (!visited[id]) topGainList.push_back(id);
    if ((int)topGainList.size() >= topPool) break;
  }

  vector<int> cands;
  cands.reserve(topPriceP + topGainP + nearLimit + 24);
  vector<pair<double, int>> scored;
  scored.reserve(topPriceP + topGainP + nearLimit + 24);
  vector<BeamState> beam, nextBeam;
  beam.reserve(beamWidth * 2);
  nextBeam.reserve(beamWidth * 4);

  for (int L : Lcand) {
    if (L > MAX_ROUTE) continue;
    beam.clear();
    BeamState root;
    beam.push_back(root);
    bool ok = true;

    for (int step = 0; step < L; ++step) {
      int k_carry = L - step;
      double mult = saleMul[sold_count + step];
      nextBeam.clear();

      for (const BeamState &st : beam) {
        bumpToken(candToken, candMark);
        cands.clear();

        int added = 0;
        for (int id : topPriceList) {
          if (visited[id] || candMark[id] == candToken) continue;
          candMark[id] = candToken;
          cands.push_back(id);
          if (++added >= topPriceP) break;
        }

        added = 0;
        for (int id : topGainList) {
          if (visited[id] || candMark[id] == candToken) continue;
          candMark[id] = candToken;
          cands.push_back(id);
          if (++added >= topGainP) break;
        }

        if (CELLS_X > 0 && CELLS_Y > 0) {
          int gx = (int)((st.cx - MIN_X) / CELL_SIZE);
          int gy = (int)((st.cy - MIN_Y) / CELL_SIZE);
          gx = max(0, min(CELLS_X - 1, gx));
          gy = max(0, min(CELLS_Y - 1, gy));

          int rad = 0;
          int found_near = 0;
          while (found_near < nearLimit && rad < 50) {
            int r_min = max(0, gy - rad);
            int r_max = min(CELLS_Y - 1, gy + rad);
            int c_min = max(0, gx - rad);
            int c_max = min(CELLS_X - 1, gx + rad);

            for (int r = r_min; r <= r_max; ++r) {
              for (int c = c_min; c <= c_max; ++c) {
                if (rad > 0 && r > r_min && r < r_max && c > c_min && c < c_max) continue;
                int g_idx = c * CELLS_Y + r;
                if (g_idx < 0 || g_idx >= (int)grid_vec.size()) continue;
                for (int id : grid_vec[g_idx]) {
                  if (visited[id] || candMark[id] == candToken) continue;
                  candMark[id] = candToken;
                  cands.push_back(id);
                  if (++found_near >= nearLimit) goto end_near2;
                }
              }
            }
            rad++;
          }
        end_near2:;
        }

        scored.clear();
        for (int id : cands) {
          bool usedInPath = false;
          for (int i = 0; i < st.len; ++i) {
            if (st.path[i] == id) {
              usedInPath = true;
              break;
            }
          }
          if (usedInPath) continue;

          double d = euclid(st.cx, st.cy, cities[id].x, cities[id].y);
          double gain = mult * cities[id].p - d * (1.0 + k_carry * C_cost);
          if (returnToHQ) gain -= (cities[id].distHQ - st.distHQ);
          scored.push_back({st.profit + gain, id});
        }

        if (scored.empty()) continue;
        if ((int)scored.size() > branchPerState) {
          nth_element(scored.begin(), scored.begin() + branchPerState, scored.end(),
                      [](const pair<double, int> &a, const pair<double, int> &b) {
                        return a.first > b.first;
                      });
          scored.resize(branchPerState);
        }

        for (const auto &pr : scored) {
          int id = pr.second;
          BeamState ns = st;
          ns.profit = pr.first;
          ns.cx = cities[id].x;
          ns.cy = cities[id].y;
          ns.distHQ = cities[id].distHQ;
          ns.path[ns.len++] = id;
          nextBeam.push_back(ns);
        }
      }

      if (nextBeam.empty()) {
        ok = false;
        break;
      }
      if ((int)nextBeam.size() > beamWidth) {
        nth_element(nextBeam.begin(), nextBeam.begin() + beamWidth, nextBeam.end(),
                    [](const BeamState &a, const BeamState &b) { return a.profit > b.profit; });
        nextBeam.resize(beamWidth);
      }
      beam.swap(nextBeam);
    }

    if (!ok) continue;
    for (const BeamState &st : beam) {
      if (st.len != L) continue;
      vector<int> order(st.path.begin(), st.path.begin() + L);
      bool uniqueOk = true;
      for (int i = 0; i < L && uniqueOk; ++i) {
        for (int j = i + 1; j < L; ++j) {
          if (order[i] == order[j]) {
            uniqueOk = false;
            break;
          }
        }
      }
      if (!uniqueOk) continue;
      Trip optimized = optimizeTrip(order, returnToHQ, sold_count);
      if (optimized.profit > best.profit) best = std::move(optimized);
    }
  }
  if (N <= 1500) {
    Trip greedy = planTripGreedy(returnToHQ);
    if (greedy.profit > best.profit) best = std::move(greedy);
  }
  return best;
}

void solve() {
  cities.clear();
  cities.resize(N);
  visited.assign(N, 0);
  sold_count = 0;
  dirty_count = 0;

  // 1. Read Cities & Find Bounds
  MIN_X = 1e9;
  MAX_X = -1e9;
  MIN_Y = 1e9;
  MAX_Y = -1e9;

  for (int i = 0; i < N; ++i) {
    if (scanf("%d %d %d", &cities[i].x, &cities[i].y, &cities[i].p) != 3) {
      break;
    }
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
  int hw = (int)thread::hardware_concurrency();
  int threads = (hw > 0 ? hw : 4);
  threads = max(1, min(threads, 8));

  if (threads > 1 && N >= 20000) {
    vector<thread> workers;
    workers.reserve(threads);
    int chunk = (N + threads - 1) / threads;
    for (int t = 0; t < threads; ++t) {
      int lo = t * chunk;
      int hi = min(N, lo + chunk);
      if (lo >= hi) break;
      workers.emplace_back([&, lo, hi]() {
        for (int i = lo; i < hi; ++i) {
          cities[i].distHQ = euclid(cities[i].x, cities[i].y, 0, 0);
        }
      });
    }
    for (auto &th : workers) th.join();
  } else {
    for (int i = 0; i < N; ++i) {
      cities[i].distHQ = euclid(cities[i].x, cities[i].y, 0, 0);
    }
  }

  candMark.assign(N, 0);
  candToken = 1;
  byPrice.resize(N);
  iota(byPrice.begin(), byPrice.end(), 0);
  byGain.resize(N);
  iota(byGain.begin(), byGain.end(), 0);

  sort(byPrice.begin(), byPrice.end(), [&](int a, int b) {
    if (cities[a].p != cities[b].p)
      return cities[a].p > cities[b].p;
    return cities[a].distHQ < cities[b].distHQ;
  });
  vector<double> gainScore(N, 0.0);
  if (threads > 1 && N >= 20000) {
    vector<thread> workers;
    workers.reserve(threads);
    int chunk = (N + threads - 1) / threads;
    for (int t = 0; t < threads; ++t) {
      int lo = t * chunk;
      int hi = min(N, lo + chunk);
      if (lo >= hi) break;
      workers.emplace_back([&, lo, hi]() {
        for (int i = lo; i < hi; ++i) {
          gainScore[i] = cities[i].p - (2.0 + C_cost) * cities[i].distHQ;
        }
      });
    }
    for (auto &th : workers) th.join();
  } else {
    for (int i = 0; i < N; ++i) {
      gainScore[i] = cities[i].p - (2.0 + C_cost) * cities[i].distHQ;
    }
  }
  sort(byGain.begin(), byGain.end(), [&](int a, int b) { return gainScore[a] > gainScore[b]; });

  while (true) {

    Trip t = planTrip(true);
    if (t.order.empty() || t.profit <= 1e-9)
      break;

    vector<int> clean;
    clean.reserve(t.order.size());
    for (int id : t.order) {
      if (visited[id]) continue;
      bool inClean = false;
      for (int x : clean) {
        if (x == id) {
          inClean = true;
          break;
        }
      }
      if (!inClean) clean.push_back(id);
    }
    if (clean.empty()) break;

    for (size_t i = 0; i < clean.size(); ++i) {
      int id = clean[i];
      visited[id] = 1;
      sold_count++;
      if (i == 0)
        printf("%d %d %d\n", cities[id].x, cities[id].y, (int)clean.size());
      else
        printf("%d %d\n", cities[id].x, cities[id].y);
    }
    printf("0 0\n");
    dirty_count += clean.size();
  }

  Trip last = planTrip(false);
  if (!last.order.empty() && last.profit > 1e-9) {
    vector<int> clean;
    clean.reserve(last.order.size());
    for (int id : last.order) {
      if (visited[id]) continue;
      bool inClean = false;
      for (int x : clean) {
        if (x == id) {
          inClean = true;
          break;
        }
      }
      if (!inClean) clean.push_back(id);
    }

    for (size_t i = 0; i < clean.size(); ++i) {
      int id = clean[i];
      visited[id] = 1;
      sold_count++;
      if (i == 0)
        printf("%d %d %d\n", cities[id].x, cities[id].y, (int)clean.size());
      else
        printf("%d %d\n", cities[id].x, cities[id].y);
    }
  }

  if (sold_count == 0 && N > 0) {
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
  }
}

int main() {
  while (scanf("%d %lf %lf", &N, &C_cost, &D_decay) == 3) {
    solve();
  }
  return 0;
}
