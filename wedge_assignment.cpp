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
#include <algorithm>
#include <iomanip>

using namespace std;

static const int    N_WEDGES    = 128;
static const int    TARGET_VOL  = 8;
static const double TWO_PI      = 2.0 * M_PI;
static const double WEDGE_WIDTH = TWO_PI / N_WEDGES;

struct Hit {
    long long hit_id;
    double    x, y, z, phi;
    int       volume_id, layer_id, module_id, wedge_index;
};

int compute_wedge(double x, double y) {
    double phi = atan2(y, x);
    double phi_shifted = phi + M_PI;
    if (phi_shifted >= TWO_PI) 
         phi_shifted -= TWO_PI;
    int idx = static_cast<int>(phi_shifted / WEDGE_WIDTH);
    return max(0, min(idx, N_WEDGES - 1));
}

int main(int argc, char* argv[]) {
    string in_path  = (argc > 1) ? argv[1] : "hits_truth1000.csv";
    string out_path = (argc > 2) ? argv[2] : "volume8_wedge_assignments.csv";

    ifstream fin(in_path);
    if (!fin.is_open()) return 1;

    string line;
    getline(fin, line); // Skip header

    vector<Hit> hits;
    while (getline(fin, line)) {
        if (line.empty()) continue;
        istringstream ss(line);
        string tok;
        try {
            Hit h;
            getline(ss, tok, ','); h.hit_id    = stoll(tok);
            getline(ss, tok, ',');
            getline(ss, tok, ',');
            getline(ss, tok, ',');
            getline(ss, tok, ',');
            getline(ss, tok, ',');
            getline(ss, tok, ',');
            getline(ss, tok, ',');
            getline(ss, tok, ',');
            getline(ss, tok, ','); h.x         = stod(tok);
            getline(ss, tok, ','); h.y         = stod(tok);
            getline(ss, tok, ','); h.z         = stod(tok);
            getline(ss, tok, ','); h.volume_id = stoi(tok);
            getline(ss, tok, ','); h.layer_id  = stoi(tok);
            getline(ss, tok, ','); h.module_id = stoi(tok);

            if (h.volume_id != TARGET_VOL) continue;

            h.phi         = atan2(h.y, h.x);
            h.wedge_index = compute_wedge(h.x, h.y);
            hits.push_back(h);
        } catch (...) {}
    }
    fin.close();

    sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) {
        if (a.wedge_index != b.wedge_index) return a.wedge_index < b.wedge_index;
        if (a.layer_id    != b.layer_id)    return a.layer_id    < b.layer_id;
        return a.z < b.z;
    });

    ofstream fout(out_path);
    if (!fout.is_open()) return 1;

    fout << "hit_id,x,y,z,volume_id,layer_id,module_id,phi_rad,wedge_index\n" << fixed << setprecision(6);
    for (const auto& h : hits) {
        fout << h.hit_id << "," << h.x << "," << h.y << "," << h.z << ","
             << h.volume_id << "," << h.layer_id << "," << h.module_id << ","
             << h.phi << "," << h.wedge_index << "\n";
    }
    fout.close();

    return 0;
}
