// TODO: Add an active transport feature, which delivers cargo to a small range at a constant rate

#include <iostream>
#include <cmath>
#include <array>
#include <fstream>
#include "Parameters.h"

std::string outputResults(const std::array<double,2*X_ELEMENTS> &conc, const int time_step);
void diffusionSolver();

int main() {

    diffusionSolver();
    std::system(TIMELAPSE);
    return 0;
}

std::string outputResults(const std::array<double,2*X_ELEMENTS> &conc, const int time_step) {

    /* Set up the filename including the correct buffer of 0's */
    std::string filename = "data/diffusion_1D_";
    std::string extension = ".txt";
    std::string number;
    if(time_step < 10) {
        number = "000" + std::to_string(time_step);
    }
    else if(time_step < 100) {
        number = "00" + std::to_string(time_step);
    }
    else if(time_step < 1000) {
        number = "0" + std::to_string(time_step);
    }
    else {
        number = std::to_string(time_step);
    }
    filename = filename + number + extension;

    /* Open file with sanity checking */
    std::ofstream output_file(filename);
    if(!output_file) {
        std::cerr << "Cannot open the output file: " << filename << std::endl;
        exit(1);
    }

    for(int i = 0; i < 2*X_ELEMENTS; i++) {
        output_file << (i-X_ELEMENTS)*DX << "\t" << conc[i] << std::endl;
    }

    output_file.close();
    return filename;
}

void diffusionSolver() {

    /* Set up starting conditions, t = 0 */
    std::array<double,2*X_ELEMENTS> conc_old = {0}, conc = {0};
    for(int i = 0; i < 2*X_ELEMENTS; i++) {
        conc_old[i] = iniDist((i - X_ELEMENTS)*DX);
    }

    std::string command_line = PY_SCRIPT + outputResults(conc_old,0) + " 0";
    std::system(command_line.c_str());

    for(int t = 1; t <= T_ELEMENTS; t++) {

        /* Update the system a single time step */
        for(int x = 0; x < 2*X_ELEMENTS - 1; x++) {
            conc[x] = conc_old[x] + ALPHA*(conc_old[x-1] - 2*conc_old[x] + conc_old[x+1]);
        }
        conc[0] = conc[1], conc[2*X_ELEMENTS - 1] = conc[2*X_ELEMENTS - 2]; // Update boundaries

        /* Output the results to a file and call the plotting script */
        if(t % STEPS_BETWEEN_PLOTS == 0) {
            command_line = PY_SCRIPT + outputResults(conc_old,t/STEPS_BETWEEN_PLOTS) + " " + std::to_string(t*DT);
            std::system(command_line.c_str());
        }

        for(int x = 0; x < 2*X_ELEMENTS; x++) {
            conc_old[x] = conc[x];
        }
    }

}