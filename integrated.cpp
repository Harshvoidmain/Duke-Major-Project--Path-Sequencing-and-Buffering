// g++ -std=c++17 integrated.cpp -o integrated.exe
// .\integrated.exe


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

double B_Tesla = 2.0; // Store the magnetic field strength in Tesla
double pT_gevc = 10.0; // Store the transverse momentum in GeV/c 
double R, r_c, r_p; // Declare variables for radius of curvature, current radius, and previous radius
double Theta[3];
double Radius[4] = {3.231, 7.215, 11.609, 17.206};
int number_of_wedges = 128;
const double M_PI = 3.14159265358979323846; // Define the mathematical constant PI

// Create a structure to store the boundary information of one wedge in one layer
struct LayerBounds {
    int layer_index;
    double wedge_start_rad; // Store the starting angle of the wedge in radians
    double wedge_end_rad; // Store the ending angle of the wedge in radians
    double wedge_start_Degrees;
    double wedge_end_Degrees;
};

// WEDGE CREATION CODE

void compute_wedges(LayerBounds layer_bounds[128][4]) {
    double phi = 2.0 * M_PI / number_of_wedges; // Calculate the angular width of one wedge

    for (int i = 0; i < number_of_wedges; i++) // Loop through all wedges from 0 to 127
    {
        layer_bounds[i][0].layer_index = 0;
        layer_bounds[i][0].wedge_start_rad = i * phi; // Calculate the starting angle of the wedge in the first layer
        layer_bounds[i][0].wedge_end_rad   = (i + 1) * phi; // Calculate the ending angle of the wedge in the first layer
        layer_bounds[i][0].wedge_start_Degrees = layer_bounds[i][0].wedge_start_rad * (180.0 / M_PI);
        layer_bounds[i][0].wedge_end_Degrees   = layer_bounds[i][0].wedge_end_rad * (180.0 / M_PI);

        for (int j = 1; j < 4; j++)  // Loop through the remaining detector layers from 0 to 3
         {
            layer_bounds[i][j].layer_index = j; // Store the current layer index
            layer_bounds[i][j].wedge_start_rad = layer_bounds[i][j-1].wedge_start_rad - Theta[j-1]; // Calculate the starting wedge boundary for the current layer by subtracting the corresponding Theta value
            layer_bounds[i][j].wedge_end_rad   = layer_bounds[i][j-1].wedge_end_rad + Theta[j-1]; // Calculate the ending wedge boundary for the current layer by adding the corresponding Theta value
            layer_bounds[i][j].wedge_start_Degrees = layer_bounds[i][j].wedge_start_rad * (180.0 / M_PI);
            layer_bounds[i][j].wedge_end_Degrees   = layer_bounds[i][j].wedge_end_rad * (180.0 / M_PI);
        }
    }
}

void calculate_theta() {
    for (int i = 1; i < 4; i++) // Loop through the layer radii starting from layer 1
     {
        r_c = Radius[i]; // Store the radius of the current layer
        r_p = Radius[i-1]; // Store the radius of the previous layer
        Theta[i-1] = std::asin((r_c) / (2.0 * R))-std::asin((r_p) / (2.0 * R));
    }
}

// LUMINOUS REGION CODE

map<int, pair<double, double>>
compute_field_bounds(const vector<double>& radii, double z0, double r_outer, double z_outer) // Receive the radius of the outermost layer (r_outer), Receive the z-coordinate corresponding to the outer layer (z_outer)
{
    double tan_a = (z_outer - z0) / r_outer; // Calculate the slope using the outer radius and z coordinates
    map<int, pair<double, double>> bounds;

    for (size_t i = 0; i < radii.size(); ++i) 
    {
        double z_max = z0 + radii[i] * tan_a; // Calculate the maximum allowed z-coordinate for the current layer
        double z_min = -z_max; // Calculate the minimum z-coordinate as the negative of z_max
        bounds[i + 1] = {z_min, z_max}; // Store the z_min and z_max values using the layer number as the key
    }

    return bounds;
}

// DATA FILTER CODE

double norm_phi(double p) {
    if (p < 0) p += 2 * M_PI; // If the angle is negative, add 2*PI to make it positive
    if (p >= 2 * M_PI) p -= 2 * M_PI;  // If the angle is greater than or equal to 2*PI, subtract 2*PI
    return p;
}

bool in_wedge(double phi, double start, double end) // To check whether an angle lies inside a wedge
 {
    if (start <= end) return phi >= start && phi <= end;
    return phi >= start || phi <= end;
}

void filter_hits(const string& input_csv, const string& output_csv,
                 LayerBounds layer_bounds[128][4],
                 const map<int, pair<double, double>>& field_bounds) // To filter hits from the input CSV file
                  {

    ifstream inputfile(input_csv);
    ofstream outputfile(output_csv);

    if (!inputfile.is_open()) {
        cerr << "Cannot open " << input_csv << "\n";
        return;
    }

    string line;
    getline(inputfile, line);

    outputfile << "hit_id,layer_index,wedge_index,x,y,z,r,phi\n";

    while (getline(inputfile, line)) // Continue reading the input file until there are no more lines
    {
        stringstream ss(line); // Create a string stream from the current CSV line
        string field; // Create a string variable to temporarily store each CSV value
        vector<string> f; // Create a vector to store all values from the current CSV row

        while (getline(ss, field, ',')) // Read each value from the line using comma as the separator
            f.push_back(field);

        // event_id,particle_id,hit_id,volume_id,layer_id,module_id,x,y,z

        if (stoi(f[3]) != 8) continue;

        int layer_id = stoi(f[4]);
        int layer_index = layer_id / 2 - 1;

        auto it = field_bounds.find(layer_index + 1);

        double z_min = it->second.first; // Get the minimum allowed z value for the current layer
        double z_max = it->second.second; // Get the maximum allowed z value for the current layer

        double x = stod(f[6]);
        double y = stod(f[7]);
        double z = stod(f[8]);
        long long hit_id = stoll(f[2]);

        if (!(z > z_min && z < z_max)) continue; // Check whether the z-coordinate is inside the luminous-region bounds

        double r = sqrt(x * x + y * y); // Calculate the radial distance of the hit from the origin
        double phi = norm_phi(atan2(y, x)); // Calculate the azimuthal angle of the hit using x and y

        for (int i = 0; i < number_of_wedges; i++) // Loop through all wedges
         {
            if (in_wedge(phi,
                         layer_bounds[i][layer_index].wedge_start_rad,
                         layer_bounds[i][layer_index].wedge_end_rad)) {

                outputfile << hit_id << "," << layer_index << "," << i << ","
                           << x << "," << y << "," << z << ","
                           << r << "," << phi << "\n";
            }
        }
    }
}

int main() {

    // 1. WEDGE CREATION

    LayerBounds layer_bounds[128][4]; // Create a 2D array to store wedge boundary information for 128 wedges and 4 detector layers

    R = (pT_gevc / (0.3 * B_Tesla)) * 100.0; // Calculate the radius of curvature R

    calculate_theta(); // Calculate the Theta values between the detector layers
    compute_wedges(layer_bounds); // Calculate all wedge boundaries and store them in layer_bounds

    ofstream MyFile("wedge_bounds.csv");

    MyFile << "Layer_Index,Wedge_Index,Wedge_Start_Rad,Wedge_End_Rad,Wedge_Start_Degrees,Wedge_End_Degrees\n";

    MyFile << fixed << setprecision(6);

    for (int i = 0; i < number_of_wedges; i++) // Loop through all 128 wedges
     {
        for (int j = 0; j < 4; j++) // Loop through all 4 detector layers
         {
            MyFile << layer_bounds[i][j].layer_index << ","
                   << i << ","
                   << layer_bounds[i][j].wedge_start_rad << ","
                   << layer_bounds[i][j].wedge_end_rad << ","
                   << layer_bounds[i][j].wedge_start_Degrees << ","
                   << layer_bounds[i][j].wedge_end_Degrees << "\n";
        }
    }

    MyFile.close();


    // 2. LUMINOUS REGION

    vector<double> radii = {3.231, 7.215, 11.609, 17.206}; // Store the radii of the four detector layers
    double z0 = 15.0;  // Define the reference z-coordinate
    double r_outer = 17.206;  // Define the radius of the outermost layer
    double z_outer = 49.078; // Define the radius of the z coordinate of the outermost layer

    auto bounds = compute_field_bounds(radii, z0, r_outer, z_outer); // Calculate z_min and z_max for every detector layer and store the results in bounds


    // 3. DATA FILTERING

    filter_hits("hits_truth1000_.csv", "filtered_wedge_hits.csv",
                layer_bounds, bounds); // Call the filtering function

    return 0;
}