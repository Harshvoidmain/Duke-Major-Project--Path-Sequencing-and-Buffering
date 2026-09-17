// https://fcrit-my.sharepoint.com/personal/megha_kolhekar_fcrit_ac_in/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fmegha%5Fkolhekar%5Ffcrit%5Fac%5Fin%2FDocuments%2FDuke%20Projects%20Sh%2D2026%2F1%2E%20Clustering%20TY%20CSE%2B%20EXTC&viewid=d558c038%2Df94d%2D455c%2D8927%2D168886314f6a&ga=1
// > g++ -std=c++17 integrated.cpp -o integrated.exe
// > ./integrated.exe

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <string>
#include <map>
#include <utility>

using namespace std;


// Create a structure to store the boundary information of one wedge in one layer
struct LayerBounds {
    int    layer_index;
    double wedge_start_rad;      // Store the starting angle of the wedge in radians
    double wedge_end_rad;        // Store the ending   angle of the wedge in radians
    double wedge_start_Degrees;  // Store the starting angle of the wedge in degrees
    double wedge_end_Degrees;    // Store the ending   angle of the wedge in degrees
};


//Common Variable Declaration begins here***
const  int    nLayers = 5;  //Layer 0 (beam axis) plus 4 physical layers
const  int    number_of_wedges = 128;
const  double M_PI = 3.14159265358979323846;

double Radius[nLayers]  = {0.0,   3.231,  7.215, 11.609, 17.206};
double zLength[nLayers] = {15.0, 49.078, 49.078, 49.078, 49.078};

LayerBounds layer_bounds[128][nLayers];

double B_Tesla = 2.0; // Store the magnetic field strength in Tesla
double pT_Gev = 10.0; // Store the transverse momentum in GeV/c 

//Declare variables for radius of curvature, Rad to Degree conversion
double R_Curv, rad2deg;

//Deflection Angle
double Deflection_Phi[nLayers];

//Common Variable Declaration ends here***



// Calculate Luminous Region
map<int, pair<double, double>> Compute_Field_Bounds() {
    map<int, pair<double, double>> bounds;

    // Calculate the slope using the outer radius and z coordinates
    double tan_a = (zLength[4] - zLength[0])/Radius[4];

    // Calculate and store the maximum and minimum allowed z-coordinate
    for (int iLy=1; iLy < nLayers; ++iLy){
        double z_max = zLength[0] + Radius[iLy] * tan_a; 
        double z_min = -z_max;
        bounds[iLy] = {z_min, z_max};
    }
    return bounds;
}



// Calculate deflection angle for each layer w.r.t. Layer 1
void Calculate_dPhi(){
    Deflection_Phi[0] = 0.;
    const double layer_1_phi = asin(Radius[1] / (2.0 * R_Curv));
    for (int iLy = 1; iLy < nLayers; iLy++){
        // Deflection angle of the current layer with respect to Layer 1.
        Deflection_Phi[iLy] =
            asin(Radius[iLy] / (2.0 * R_Curv)) - layer_1_phi;
    }
}



//Wedge creation: Calculate the angular width of one wedge in radians
void Compute_Wedges() {
    double dphi = 2.0 * M_PI/number_of_wedges;

    // Loop through all wedges from 0 to 127 and Layers 1 to 4
    for (int i=0; i<number_of_wedges; i++){
        for (int iLy=1; iLy<5; iLy++){
            layer_bounds[i][iLy].layer_index = iLy;  // Current Layer Index

            // Calculate the starting and ending wedge boundary in radians
            double phi_start = i * dphi     - Deflection_Phi[iLy];
            double phi_end   = (i+1) * dphi + Deflection_Phi[iLy];

            if(phi_start<0){phi_start  = 2*M_PI  + phi_start;}
            if(phi_end>2*M_PI){phi_end = phi_end - 2*M_PI;}

            layer_bounds[i][iLy].wedge_start_rad     = phi_start;           //Radians
            layer_bounds[i][iLy].wedge_end_rad       = phi_end;             //Radians
            layer_bounds[i][iLy].wedge_start_Degrees = phi_start * rad2deg; //Degrees
            layer_bounds[i][iLy].wedge_end_Degrees   = phi_end * rad2deg;   //Degrees
        }  // for (int iLy=1; iLy<5; iLy++
    } // for (int i=0; i<number_of_wedges


    //Stroe the Wedge boundary data
    ofstream MyFile("wedge_bounds.csv");
    MyFile << "Layer_Index,Wedge_Index,Wedge_Start_Rad,Wedge_End_Rad,Wedge_Start_Degrees,Wedge_End_Degrees\n";
    MyFile << fixed << setprecision(6);

    //Loop through all 128 wedges and 4 detector layers
    for (int i = 0; i < number_of_wedges; i++){
        for (int iLy = 1; iLy < 5; iLy++){
            int il = layer_bounds[i][iLy].layer_index;
            double phi_min_rad = layer_bounds[i][iLy].wedge_start_rad;
            double phi_max_rad = layer_bounds[i][iLy].wedge_end_rad;
            double phi_min_deg = layer_bounds[i][iLy].wedge_start_Degrees;
            double phi_max_deg = layer_bounds[i][iLy].wedge_end_Degrees;

            MyFile << il << "," << i << "," << phi_min_rad << ","
                   << phi_max_rad << "," << phi_min_deg << "," << phi_max_deg << endl;
        }
    }
    MyFile.close();
}



double Norm_Phi(double p) {
    // If the angle is negative, add 2*PI to make it positive
    if (p < 0) p += 2 * M_PI;

    // If the angle is greater than or equal to 2*PI, subtract 2*PI
    if (p >= 2 * M_PI) p -= 2 * M_PI; 
    return p;
}



// To check whether an angle lies inside a wedge
bool In_Wedge(double phi, double start, double end) {
    if (start <= end) return phi >= start && phi <= end;
    return phi >= start || phi <= end;
}



// To filter hits from the input CSV file
//Input File Data Format
//0:event_id,  1:particle_id,  2:hit_id,  3:volume_id,     4:layer_id,  5:module_id,
//6:x,  7:y,   8:z,  9:tx,    10:ty,     11:tz,  12:tpx,  13:tpy,      14:tpz,  15:weight
void Filter_Hits(const map<int, pair<double, double>>& field_bounds){
    ifstream inputfile("hits_truth1000 (1).csv");
    ofstream outputfile("filtered_wedge_hits.csv");

    if (!inputfile.is_open()) {
        cerr << "Cannot open input file" << "\n";
        return;
    }

    outputfile << "event,particle_id,hit_id,layer,wedge_id,module_id,x,y,z,tpx,tpy,tpz,phi\n";
    // Read the input file
    string line;    getline(inputfile, line);
    while (getline(inputfile, line)){
        string field;            // Temporarily stores each CSV value
        vector<string> f;        // Stores all the values from the current CSV line        
        stringstream ss(line);   // Create a string stream from the current CSV line

        // Read each value from the line using comma as the separator
        while (getline(ss, field, ',')) {f.push_back(field);}
        if (stoi(f[3]) != 8) continue; //Select hits with volume ID == 8

        int iEvent = stod(f[0]), iLayer = stoi(f[4])/2,   module_id = stoi(f[5]);
        double x   = stod(f[6]) * 0.1,   y   = stod(f[7]) * 0.1,   z   = stod(f[8]) * 0.1; // mm to cm
        double tpx = stod(f[12]),  tpy = stod(f[13]),  tpz = stod(f[14]);
        long long hit_id = stoll(f[2]),   part_id = stoll(f[1]);

        auto it = field_bounds.find(iLayer);
        // Get the minimum and maximum allowed z value
        double z_min = it->second.first,   z_max = it->second.second; 

        // Check whether the z-coordinate is inside the luminous-region bounds
        if (!(z > z_min && z < z_max)) continue; 

        // Calculate the azimuthal angle of the hit using x and y
        double phi = Norm_Phi(atan2(y, x));

        // Loop through all Wedges
        for (int i = 0; i < number_of_wedges; i++){
            double phi_start = layer_bounds[i][iLayer].wedge_start_rad;
            double phi_end   = layer_bounds[i][iLayer].wedge_end_rad;
            bool checkPhi    = In_Wedge(phi, phi_start, phi_end);
            if (checkPhi){
                outputfile << iEvent << "," << part_id << "," << hit_id << "," 
                << iLayer << "," << i << "," << module_id << "," << x << "," 
                << y << "," << z  << "," << tpx << "," << tpy << "," << tpz << ","
                << phi*rad2deg << "\n";
            }
        } // for (int i = 0; i < number_of_wedge
    }    // while (getline(inputfile
}



#include <algorithm>
#include <array>

const int Hits_Per_Superpoint = 16;

struct SeedHit {
    long long event, particle_id, hit_id;
    int layer, wedge_id;
    long long module_id;
    double x, y, z, tpx, tpy, tpz, phi; // cm; phi in degrees
};
struct Superpoint {
    vector<SeedHit> hits;
    double z_min = 0, z_max = 0;
};
struct PatchCorner { double z1, z4; };

// The supplied list is sorted by increasing z.
bool Select_Superpoint(const vector<SeedHit>& hits, double target, Superpoint& sp) {
    int end = 0;
    while (end < (int)hits.size() && hits[end].z <= target) ++end;
    if (end < Hits_Per_Superpoint) return false;
    sp.hits.assign(hits.begin() + end - Hits_Per_Superpoint, hits.begin() + end);
    sp.z_min = sp.hits.front().z;
    sp.z_max = sp.hits.back().z;
    return true;
}

// Intersection with a*z1 + b*z4 <= c.
vector<PatchCorner> Clip_Seed_Patch(const vector<PatchCorner>& polygon,
                                   double a, double b, double c) {
    vector<PatchCorner> result;
    for (size_t i = 0; i < polygon.size(); ++i) {
        auto p = polygon[i], q = polygon[(i+1) % polygon.size()];
        double dp = a*p.z1 + b*p.z4 - c, dq = a*q.z1 + b*q.z4 - c;
        if (dp <= 0) result.push_back(p);
        if ((dp <= 0) != (dq <= 0)) {
            double t = dp / (dp-dq);
            result.push_back({p.z1+t*(q.z1-p.z1), p.z4+t*(q.z4-p.z4)});
        }
    }
    return result;
}

// Stage 1: ONE first seed patch per event/wedge, starting from high z.
// Input coordinates must already be cm; no filtering or unit conversion here.
void Form_Seed_Patches() {
    ifstream input("filtered_wedge_hits.csv");
    if (!input.is_open()) {
        cerr << "Cannot open filtered_wedge_hits.csv\n";
        return;
    }

    // Group each event and wedge into separate layer lists.
    map<pair<long long,int>, array<vector<SeedHit>,nLayers>> groups;
    string line;
    getline(input, line); // Skip the CSV header.
    while (getline(input, line)) {
        stringstream ss(line);
        string field;
        vector<string> f;
        while (getline(ss, field, ',')) f.push_back(field);
        if (f.size() != 13) continue;

        SeedHit h;
        h.event = stoll(f[0]); h.particle_id = stoll(f[1]); h.hit_id = stoll(f[2]);
        h.layer = stoi(f[3]); h.wedge_id = stoi(f[4]); h.module_id = stoll(f[5]);
        h.x = stod(f[6]); h.y = stod(f[7]); h.z = stod(f[8]);
        h.tpx = stod(f[9]); h.tpy = stod(f[10]); h.tpz = stod(f[11]);
        h.phi = stod(f[12]);
        if (h.layer < 1 || h.layer >= nLayers) continue;
        groups[{h.event,h.wedge_id}][h.layer].push_back(h);
    }
    input.close();

    ofstream hits_out("seedpatch_hits.csv"), corners_out("seedpatch_corners.csv");
    if (!hits_out || !corners_out) {
        cerr << "Cannot create patch output files\n";
        return;
    }
    hits_out << setprecision(17)
        << "patch_id,event,wedge_id,layer,hit_index,particle_id,hit_id,module_id,"
           "x_cm,y_cm,z_cm,tpx,tpy,tpz,phi_deg,sp_z_min_cm,sp_z_max_cm\n";
    corners_out << setprecision(17)
        << "patch_id,event,wedge_id,corner_index,z1_cm,z4_cm\n";
    int patch_id = 0, insufficient = 0, empty = 0;
    for (auto& entry : groups) {
        auto& layers = entry.second;
        bool enough = true;
        for (int layer = 1; layer < nLayers; ++layer) {
            auto& hits = layers[layer];
            sort(hits.begin(), hits.end(), [](const SeedHit& a, const SeedHit& b) {
                return a.z < b.z || (a.z == b.z && a.hit_id < b.hit_id);
            });
            if (hits.size() < Hits_Per_Superpoint) enough = false;
        }
        if (!enough) { ++insufficient; continue; }
        Superpoint sp[nLayers];
        Select_Superpoint(layers[1], layers[1].back().z, sp[1]);
        Select_Superpoint(layers[4], layers[4].back().z, sp[4]);
        bool valid = true;
        for (int layer = 2; layer <= 3; ++layer) {
            double alpha = (Radius[layer]-Radius[1])/(Radius[4]-Radius[1]);
            double target = (1-alpha)*sp[1].z_max + alpha*sp[4].z_max;
            // Right-justify: take the last 16 hits at or below the max-to-max line.
            if (!Select_Superpoint(layers[layer], target, sp[layer])) valid = false;
        }
        if (!valid) { ++insufficient; continue; }
        vector<PatchCorner> polygon = {
            {sp[1].z_min,sp[4].z_min}, {sp[1].z_max,sp[4].z_min},
            {sp[1].z_max,sp[4].z_max}, {sp[1].z_min,sp[4].z_max}
        };
        for (int layer = 2; layer <= 3; ++layer) {
            double alpha = (Radius[layer]-Radius[1])/(Radius[4]-Radius[1]);
            polygon = Clip_Seed_Patch(polygon, 1-alpha, alpha, sp[layer].z_max);
            polygon = Clip_Seed_Patch(polygon, alpha-1, -alpha, -sp[layer].z_min);
        }
        double area2 = 0;
        for (size_t i = 0; i < polygon.size(); ++i) {
            auto p = polygon[i], q = polygon[(i+1)%polygon.size()];
            area2 += p.z1*q.z4-q.z1*p.z4;
        }
        if (polygon.size() < 3 || fabs(area2) <= 1e-12) { ++empty; continue; }
        for (int layer = 1; layer < nLayers; ++layer) {
            for (size_t i = 0; i < sp[layer].hits.size(); ++i) {
                const auto& h = sp[layer].hits[i];
                hits_out << patch_id << ',' << h.event << ',' << h.wedge_id << ','
                    << layer << ',' << i << ',' << h.particle_id << ',' << h.hit_id
                    << ',' << h.module_id << ',' << h.x << ',' << h.y << ',' << h.z
                    << ',' << h.tpx << ',' << h.tpy << ',' << h.tpz << ',' << h.phi
                    << ',' << sp[layer].z_min << ',' << sp[layer].z_max << '\n';
            }
        }
        for (size_t i = 0; i < polygon.size(); ++i)
            corners_out << patch_id << ',' << entry.first.first << ','
                << entry.first.second << ',' << i << ',' << polygon[i].z1
                << ',' << polygon[i].z4 << '\n';
        ++patch_id;
    }
    hits_out.close(); corners_out.close();
    cout << "First seed patches: " << patch_id
         << "\nInsufficient selectable hits: " << insufficient
         << "\nEmpty/zero-area patches: " << empty << '\n';
}

int main(){
    rad2deg = 180.0 / M_PI;
    // Calculate the radius of curvature R
    R_Curv = (pT_Gev / (0.3 * B_Tesla)) * 100.0;

    // Calculate LUMINOUS REGION by calculating z_min and z_max for every layer
    auto bounds = Compute_Field_Bounds();

    // Calculate the deflection angle w.r.t. Layer 1
    Calculate_dPhi();

    // Calculate Wedges boundaries for all layers and store data i a csv file
    Compute_Wedges();

    // Data Filtering
    Filter_Hits(bounds);

    //patch formation 
    Form_Seed_Patches();
    return 0;
}
