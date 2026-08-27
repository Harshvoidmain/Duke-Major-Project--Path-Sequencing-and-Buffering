#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>

using namespace std;

static const int    L = 5;
static const double R[L] = {5.0, 10.0, 15.0, 20.0, 25.0};   


double SEED_SUPERPOINTS[L][2] = {
    { 2.0,  22.0},   
    { 8.0,  31.0}, 
    { 1.9,  40.0},   
    { 1.8,  44.0},   
    {-4.6,  50.4},   
};

struct Vec2 { double x, y; };
using Poly = vector<Vec2>;


struct Corners { double z_a, z_b, z_c, z_d; };

Poly clip_hp(Poly p, double a, double b, double c) {
    Poly res;
    int n = (int)p.size();
    for (int i = 0; i < n; i++) {
        double xi = p[i].x, yi = p[i].y;
        double xj = p[(i+1)%n].x, yj = p[(i+1)%n].y;
        double di = a*xi + b*yi - c;
        double dj = a*xj + b*yj - c;
        if (di <= 1e-9) res.push_back({xi, yi});
        if ((di < -1e-9 && dj > 1e-9) || (di > 1e-9 && dj < -1e-9)) {
            double t = di / (di - dj);
            res.push_back({xi + t*(xj-xi), yi + t*(yj-yi)});
        }
    }
    return res;
} 

inline double alpha_l(int l) { return (R[l] - R[0]) / (R[L-1] - R[0]); }

Poly compute_poly(const double sp[L][2]) {
    Poly p = {{-300,-300},{300,-300},{300,300},{-300,300}};
    for (int l = 0; l < L; l++) {
        double zmax = sp[l][1];
        double zmin = sp[l][0];
        p = clip_hp(p,  (1.0 - alpha_l(l)),  alpha_l(l),  zmax);  if (p.empty()) return p;
        p = clip_hp(p, -(1.0 - alpha_l(l)), -alpha_l(l), -zmin);  if (p.empty()) return p;
    }
    return p;
}

Corners get_corners(const Poly& p, double z1_min, double z1_max) {
    Corners c = {-1e9, -1e9, 1e9, 1e9}; 
    for (auto& v : p) {
        if (abs(v.x - z1_min) < 1e-5) {
            c.z_a = max(c.z_a, v.y); // top-left
            c.z_c = min(c.z_c, v.y); // bottom-left
        }
        if (abs(v.x - z1_max) < 1e-5) {
            c.z_b = max(c.z_b, v.y); // top-right
            c.z_d = min(c.z_d, v.y); // bottom-right
        }
    }
    return c;
}

int main() {
    cout << fixed << setprecision(6);

    Poly seed_patch = compute_poly(SEED_SUPERPOINTS);
    if (seed_patch.empty()) {
        cout << "Error: Seed patch is empty.\n";
        return 1;
    }
    Corners seed_c = get_corners(seed_patch, SEED_SUPERPOINTS[0][0], SEED_SUPERPOINTS[0][1]);

    double W[L];
    for (int l = 0; l < L; l++) {
        W[l] = SEED_SUPERPOINTS[l][1] - SEED_SUPERPOINTS[l][0];
    }
    double best_zLc = 0;
    Poly comp_patch;
    double comp_sp[L][2];
    double min_overlap = 1e9;
    double z1_c = SEED_SUPERPOINTS[0][0]; 
    
    for (double trial_zLc = -60.0; trial_zLc <= 60.0; trial_zLc += 0.1) {
        double trial_sp[L][2];
        for (int l = 0; l < L; l++) {
            double zl_min = (1.0 - alpha_l(l)) * z1_c + alpha_l(l) * trial_zLc;
            trial_sp[l][0] = zl_min;
            trial_sp[l][1] = zl_min + W[l];
        }
        Poly p = compute_poly(trial_sp);
        if (p.empty()) continue;
        Corners p_c = get_corners(p, trial_sp[0][0], trial_sp[0][1]);
        double delta_ac = p_c.z_a - seed_c.z_c;
        double delta_bd = p_c.z_b - seed_c.z_d;
        if (delta_ac > 0 && delta_bd > 0) {
            double overlap = min(delta_ac, delta_bd);
            if (overlap < min_overlap) {
                min_overlap = overlap;
                best_zLc = trial_zLc;
                comp_patch = p;
                for(int l=0; l<L; l++) {
                    comp_sp[l][0] = trial_sp[l][0];
                    comp_sp[l][1] = trial_sp[l][1];
                }
            }
        }
    }

    if (comp_patch.empty()) {
        cout << "Error: Could not form a valid overlapping complementary patch.\n";
        return 1;
    }

    cout << "=== COMPLEMENTARY PATCH SUPERPOINTS ===\n";
    cout << left << setw(6) << "layer" << setw(12) << "z_min" << setw(12) << "z_max" << "\n";
    for (int l = 0; l < L; l++) {
        cout << left << setw(6) << l << setw(12) << comp_sp[l][0] << setw(12) << comp_sp[l][1] << "\n";
    }


    double z1_min = -30, z1_max = 35, zL_min = -60, zL_max = 60;
    double PX0 = 60, PY0 = 30, PW = 560, PH = 460;
    auto X = [&](double z1){ return PX0 + (z1 - z1_min)/(z1_max-z1_min)*PW; };
    auto Y = [&](double zL){ return PY0 + PH - (zL - zL_min)/(zL_max-zL_min)*PH; };

    ofstream svg("superpatch.svg");
    svg << "<?xml version=\"1.0\"?>\n<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"680\" height=\"540\" viewBox=\"0 0 680 540\">\n";
    svg << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    svg << "<rect x=\"" << PX0 << "\" y=\"" << PY0 << "\" width=\"" << PW << "\" height=\"" << PH << "\" fill=\"none\" stroke=\"#888\"/>\n";
    
    // Axes
    svg << "<line x1=\"" << PX0 << "\" y1=\"" << Y(0) << "\" x2=\"" << PX0+PW << "\" y2=\"" << Y(0) << "\" stroke=\"#aaa\"/>\n";
    svg << "<line x1=\"" << X(0) << "\" y1=\"" << PY0 << "\" x2=\"" << X(0) << "\" y2=\"" << PY0+PH << "\" stroke=\"#aaa\"/>\n";
    for (int v = -30; v <= 30; v += 10)
        svg << "<text x=\"" << X(v) << "\" y=\"" << PY0+PH+16 << "\" font-size=\"10\" text-anchor=\"middle\">" << v << "</text>\n";
    for (int v = -60; v <= 60; v += 20)
        svg << "<text x=\"" << PX0-8 << "\" y=\"" << Y(v)+3 << "\" font-size=\"10\" text-anchor=\"end\">" << v << "</text>\n";
    
    svg << "<text x=\"" << PX0+PW/2 << "\" y=\"" << PY0+PH+32 << "\" font-size=\"11\" text-anchor=\"middle\">z1 (cm)</text>\n";
    svg << "<text x=\"15\" y=\"" << PY0+PH/2 << "\" font-size=\"11\" text-anchor=\"middle\" transform=\"rotate(-90,15," << PY0+PH/2 << ")\">zL (cm)</text>\n";
    
    // Draw Seed Patch (Orange)
    svg << "<polygon points=\"";
    for (auto& v : seed_patch) svg << X(v.x) << "," << Y(v.y) << " ";
    svg << "\" fill=\"#D85A30\" fill-opacity=\"0.5\" stroke=\"#993C1D\" stroke-width=\"1.2\"/>\n";

    // Draw Complementary Patch (Blue)
    svg << "<polygon points=\"";
    for (auto& v : comp_patch) svg << X(v.x) << "," << Y(v.y) << " ";
    svg << "\" fill=\"#305AD8\" fill-opacity=\"0.5\" stroke=\"#1D3C99\" stroke-width=\"1.2\"/>\n";
    
    svg << "</svg>\n";
    svg.close();
    
    cout << "\nGenerated 'superpatch.svg' mapping both the Seed Patch (Orange) and Complementary Patch (Blue).\n";
    return 0;
}