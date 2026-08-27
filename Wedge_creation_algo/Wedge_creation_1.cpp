#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

double B_Tesla = 2.0;    //taken from paper 
double pT_gevc = 10.0;   //taken from paper 
double R, r_c, r_p;
double Theta[3];
double Radius[4] = {3.231, 7.215, 11.609, 17.206};
int number_of_wedges = 128;     //taken from paper 
const double M_PI = 3.14159265358979323846;

struct LayerBounds {
    int layer_index;
    double wedge_start_rad;
    double wedge_end_rad;
    double wedge_start_Degrees;
    double wedge_end_Degrees;
};

void compute_wedges(LayerBounds layer_bounds[128][4]) {
    double phi = 2.0 * M_PI / number_of_wedges;       //phi is the angle of each wedge in layer 0
    for (int i = 0; i < number_of_wedges; i++) {      //wedge loop
        layer_bounds[i][0].layer_index = 0;
        layer_bounds[i][0].wedge_start_rad = i * phi;
        layer_bounds[i][0].wedge_end_rad   = (i + 1) * phi;
        layer_bounds[i][0].wedge_start_Degrees = layer_bounds[i][0].wedge_start_rad * (180.0 / M_PI);
        layer_bounds[i][0].wedge_end_Degrees   = layer_bounds[i][0].wedge_end_rad * (180.0 / M_PI);
        
        for (int j = 1; j < 4; j++) {     //layer loop
            layer_bounds[i][j].layer_index = j;
            layer_bounds[i][j].wedge_start_rad = layer_bounds[i][j-1].wedge_start_rad - Theta[j-1];
            layer_bounds[i][j].wedge_end_rad   = layer_bounds[i][j-1].wedge_end_rad + Theta[j-1];
            layer_bounds[i][j].wedge_start_Degrees = layer_bounds[i][j].wedge_start_rad * (180.0 / M_PI);
            layer_bounds[i][j].wedge_end_Degrees   = layer_bounds[i][j].wedge_end_rad * (180.0 / M_PI);
        }
    }
}

void calculate_theta() {                                                            //theta is the extra angle that the wedge needs to be extended to account for the curvature of the particle trajectory
    for (int i = 1; i < 4; i++) {
        r_c = Radius[i];
        r_p = Radius[i-1];
        Theta[i-1] = std::asin((r_c) / (2.0 * R))-std::asin((r_p) / (2.0 * R));
    }
}

int main() {
    LayerBounds layer_bounds[128][4];

    R = (pT_gevc / (0.3 * B_Tesla)) * 100.0;
    calculate_theta();
    compute_wedges(layer_bounds);

                //for (int i = 0; i < 128; i++) {
                //    for (int j = 0; j < 4; j++) {//        std::cout << "Layer " << j 
                //                 << "Start (deg): " << layer_bounds[i][j].wedge_start_rad
                //                 << "End (deg): " << layer_bounds[i][j].wedge_end_rad
                //                << "Width (deg): " << (layer_bounds[i][j].wedge_end_rad - layer_bounds[i][j].wedge_start_rad)
                //                << "\n";
                //    }
                //}

    std::ofstream MyFile("wedge_bounds.csv");
    MyFile << "Layer_Index,Wedge_Index,Wedge_Start_Rad,Wedge_End_Rad,Wedge_Start_Degrees,Wedge_End_Degrees\n";
    MyFile << std::fixed << std::setprecision(6);
    for (int i = 0; i < number_of_wedges; i++) {
    for (int j = 0; j < 4; j++) {
        MyFile << layer_bounds[i][j].layer_index << ","
             << i << ","                                           //wedge 
             << layer_bounds[i][j].wedge_start_rad << ","
             << layer_bounds[i][j].wedge_end_rad << ","
             << layer_bounds[i][j].wedge_start_Degrees << ","
             << layer_bounds[i][j].wedge_end_Degrees << "\n";
    }
    }
    MyFile.close();

    return 0;
}