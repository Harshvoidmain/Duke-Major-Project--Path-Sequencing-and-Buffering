#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif
#include <map>
#include <set>
#include <algorithm>
#include <iomanip>

using namespace std;

static const double EPS        = 1e-6;
static const double CLIP_BOUND = 20000.0;
static const int    N_WEDGES   = 128;

struct Vec2 { double x, y; };
using Poly = vector<Vec2>;

struct Hit {
    long long hit_id;
    double    x, y, z, r, phi;
    int       volume_id, layer_id, module_id, wedge_index;
};

struct LayerData {
    int         layer_id;
    double      radius;
    vector<Hit> hits;   
};

struct Superpoint {
    double      z_min = 0, z_max = 0;
    int         start_idx = -1;
    vector<Hit> hits;
};

struct SeedPatch {
    int                global_index;
    int                wedge_index;
    int                column_index;
    int                patch_in_column;
    bool               is_rectangular;
    vector<Superpoint> superpoints;
};

// Polygon Clipping & Geometry
Poly clip_hp(const Poly& p, double a, double b, double c) {
    Poly res;
    int n = (int)p.size();
    for (int i = 0; i < n; i++) {
        const Vec2& pi = p[i];
        const Vec2& pj = p[(i + 1) % n];
        double di = a * pi.x + b * pi.y - c;
        double dj = a * pj.x + b * pj.y - c;
        if (di <= EPS) res.push_back(pi);
        if ((di < -EPS && dj > EPS) || (di > EPS && dj < -EPS)) {
            double t = di / (di - dj);
            res.push_back({pi.x + t * (pj.x - pi.x), pi.y + t * (pj.y - pi.y)});
        }
    }
    return res;
}

Poly compute_polygon(const vector<pair<double,double>>& sps, const vector<double>& alphas) {
    Poly p = {{-CLIP_BOUND, -CLIP_BOUND}, {CLIP_BOUND, -CLIP_BOUND},
              { CLIP_BOUND,  CLIP_BOUND}, {-CLIP_BOUND, CLIP_BOUND}};
    for (size_t l = 0; l < sps.size() && !p.empty(); l++) {
        double a = 1.0 - alphas[l], b = alphas[l];
        p = clip_hp(p,  a,  b, sps[l].second);
        p = clip_hp(p, -a, -b, -sps[l].first);
    }
    return p;
}

bool is_axis_aligned_rect(const Poly& p) {
    if (p.size() != 4) return false;
    double xmin = p[0].x, xmax = p[0].x, ymin = p[0].y, ymax = p[0].y;
    for (auto& v : p) {
        xmin = min(xmin, v.x); xmax = max(xmax, v.x);
        ymin = min(ymin, v.y); ymax = max(ymax, v.y);
    }
    for (auto& v : p) {
        if ((fabs(v.x - xmin) >= 1e-3 && fabs(v.x - xmax) >= 1e-3) ||
            (fabs(v.y - ymin) >= 1e-3 && fabs(v.y - ymax) >= 1e-3)) return false;
    }
    return true;
}

// Field Boundaries & Superpoint Search

struct FieldBounds { double z1_b_mm, zL_b_mm; };

FieldBounds compute_field_bounds(const vector<double>& z_in, const vector<double>& z_out, double r1, double rL) {
    auto get_pct = [](vector<double> zv) {
        sort(zv.begin(), zv.end());
        int n = (int)zv.size();
        return pair<double,double>{zv[max(0, (int)(0.01 * n))] / 10.0, zv[min(n - 1, (int)(0.99 * n))] / 10.0};
    };
    auto [z0_min, z0_max] = get_pct(z_in);
    auto [zL_min, zL_max] = get_pct(z_out);
    double ratio = r1 / rL;
    return {(z0_max * 10.0) * (1.0 - ratio) + (zL_max * 10.0) * ratio, zL_max * 10.0};
}

Superpoint find_rj_sp(const vector<Hit>& hits, double z_target, int N) {
    if ((int)hits.size() < N) return {};
    auto it = upper_bound(hits.begin(), hits.end(), z_target + EPS, [](double val, const Hit& h){ return val < h.z; });
    int ri = (int)(it - hits.begin()) - 1;
    if (ri < N - 1) return {};
    int si = ri - N + 1;
    return {hits[si].z, hits[ri].z, si, vector<Hit>(hits.begin() + si, hits.begin() + ri + 1)};
}

Superpoint find_spanning_sp(const vector<Hit>& hits, double z_need_min, double z_need_max, int N) {
    if ((int)hits.size() < N) return {};
    auto it = lower_bound(hits.begin(), hits.end(), z_need_max - EPS, [](const Hit& h, double val){ return h.z < val; });
    if (it == hits.end()) return {};
    int ri = (int)(it - hits.begin());
    int si = ri - N + 1;
    if (si < 0 || hits[si].z > z_need_min + EPS) return {};
    return {hits[si].z, hits[ri].z, si, vector<Hit>(hits.begin() + si, hits.begin() + ri + 1)};
}


// Patch Formation
vector<SeedPatch> form_seed_patches(int wedge_idx, const vector<LayerData>& layers, int N, const FieldBounds& fb, int& global_ctr) {
    vector<SeedPatch> patches;
    int L = (int)layers.size();
    if (L < 2) return patches;

    double r0 = layers[0].radius, rL = layers[L-1].radius;
    vector<double> alpha(L);
    for (int l = 0; l < L; l++) alpha[l] = (layers[l].radius - r0) / (rL - r0);

    double z1_target = fb.z1_b_mm;
    int col = 0, patch_in_col = 0;
    bool flagged = false;

    while (!flagged) {
        Superpoint sp1 = find_rj_sp(layers[0].hits, z1_target, N);
        if (sp1.start_idx < 0) break;

        double zL_target = fb.zL_b_mm;
        patch_in_col = 0;

        while (!flagged) {
            Superpoint sp_L = find_rj_sp(layers[L-1].hits, zL_target, N);
            if (sp_L.start_idx < 0) break;

            vector<Superpoint> sps(L);
            sps[0] = sp1; sps[L-1] = sp_L;
            bool analytic_rect = true;

            for (int l = 1; l < L - 1; l++) {
                double z_min = (1.0 - alpha[l]) * sp1.z_min + alpha[l] * sp_L.z_min;
                double z_max = (1.0 - alpha[l]) * sp1.z_max + alpha[l] * sp_L.z_max;
                Superpoint sp_l = find_spanning_sp(layers[l].hits, z_min, z_max, N);
                if (sp_l.start_idx < 0) {
                    analytic_rect = false;
                    sp_l = find_rj_sp(layers[l].hits, z_max, N);
                }
                sps[l] = sp_l;
            }

            vector<pair<double,double>> sp_ranges(L);
            bool valid = true;
            for (int l = 0; l < L; l++) {
                if (sps[l].start_idx < 0) { valid = false; break; }
                sp_ranges[l] = {sps[l].z_min, sps[l].z_max};
            }
            if (!valid) { flagged = true; break; }

            Poly poly = compute_polygon(sp_ranges, alpha);
            if (poly.empty()) { flagged = true; break; }

            bool is_rect = analytic_rect && is_axis_aligned_rect(poly);
            patches.push_back({global_ctr++, wedge_idx, col, patch_in_col++, is_rect, sps});

            if (is_rect) zL_target = sp_L.z_min - EPS;
            else        { flagged = true; break; }
        }
        z1_target = sp1.z_min - EPS;
        col++;
    }
    return patches;
}

int main(int argc, char* argv[]) {
    int    N           = (argc > 1) ? atoi(argv[1]) : 16;
    string input_csv   = (argc > 2) ? argv[2] : "volume8_wedge_assignments.csv";
    string output_hits = (argc > 3) ? argv[3] : "seedpatch_hits.csv";

    ifstream fin(input_csv);
    if (!fin.is_open()) return 1;
    string line; getline(fin, line);

    map<int, map<int, vector<Hit>>> raw;
    while (getline(fin, line)) {
        if (line.empty()) continue;
        istringstream ss(line);
        string tok;
        try {
            Hit h;
            getline(ss, tok, ','); h.hit_id      = stoll(tok);
            getline(ss, tok, ','); h.x         = stod(tok);
            getline(ss, tok, ','); h.y         = stod(tok);
            getline(ss, tok, ','); h.z         = stod(tok);
            getline(ss, tok, ','); h.volume_id = stoi(tok);
            getline(ss, tok, ','); h.layer_id  = stoi(tok);
            getline(ss, tok, ','); h.module_id = stoi(tok);
            getline(ss, tok, ','); h.phi       = stod(tok);
            getline(ss, tok, ','); h.wedge_index = stoi(tok);
            h.r = sqrt(h.x * h.x + h.y * h.y);
            raw[h.wedge_index][h.layer_id].push_back(h);
        } catch (...) {}
    }
    fin.close();

    set<int> lids;
    for (auto& [w, lm] : raw) for (auto& [l, _] : lm) lids.insert(l);
    vector<int> layer_order(lids.begin(), lids.end());

    map<int, pair<double, int>> r_sum;
    vector<double> z_in, z_out;
    for (auto& [w, lm] : raw) {
        for (auto& [l, hits] : lm) {
            for (auto& h : hits) {
                r_sum[l].first += h.r; r_sum[l].second++;
                if (l == layer_order.front()) z_in.push_back(h.z);
                if (l == layer_order.back())  z_out.push_back(h.z);
            }
        }
    }

    FieldBounds fb = compute_field_bounds(z_in, z_out, r_sum[layer_order.front()].first / r_sum[layer_order.front()].second,
                                                       r_sum[layer_order.back()].first / r_sum[layer_order.back()].second);

    ofstream fout(output_hits);
    if (!fout.is_open()) return 1;

    fout << "global_patch_index,wedge_index,column_index,patch_in_column,"
         << "is_rectangular,layer_id,hit_index_in_sp,hit_id,"
         << "x_mm,y_mm,z_mm,r_mm,phi_rad,module_id\n" << fixed << setprecision(3);

    int global_ctr = 0;
    for (int w = 0; w < N_WEDGES; w++) {
        if (!raw.count(w)) continue;
        vector<LayerData> layers;
        for (int lid : layer_order) {
            if (!raw[w].count(lid)) continue;
            auto hits = raw[w][lid];
            sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) { return a.z < b.z; });
            double r_mean = 0; for (auto& h : hits) r_mean += h.r;
            layers.push_back({lid, r_mean / hits.size(), hits});
        }

        for (auto& sp : form_seed_patches(w, layers, N, fb, global_ctr)) {
            for (auto& superpoint : sp.superpoints) {
                for (size_t hi = 0; hi < superpoint.hits.size(); hi++) {
                    const auto& h = superpoint.hits[hi];
                    fout << sp.global_index << "," << sp.wedge_index << "," << sp.column_index << ","
                         << sp.patch_in_column << "," << sp.is_rectangular << "," << h.layer_id << ","
                         << hi << "," << h.hit_id << "," << h.x << "," << h.y << "," << h.z << ","
                         << h.r << "," << h.phi << "," << h.module_id << "\n";
                }
            }
        }
    }
    fout.close();
    return 0;
}
