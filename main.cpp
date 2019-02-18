// TODO: Adjust active transport so it only delivers what is needed to maintain a maximum (in params) at x=0


#include <iostream>
#include <cmath>
#include <array>
#include <fstream>
#include <algorithm>
#include "Parameters.h"


/* Global variables (yuck) are suffixed with 'G', don't mess with them */
double stream_seg_remainder_G = 0.0;


std::string outputResults(const std::array<double,2*X_ELEMENTS> &conc, const int time_step);
void runDiffusion(std::array<double,2*X_ELEMENTS> &rho, std::array<double,2*X_ELEMENTS> &rho_old);
void cytoplasmicStream(std::array<double,2*X_ELEMENTS> &rho, std::array<double,2*X_ELEMENTS> &stream);
std::array<double, 2*X_ELEMENTS> visibleConc(std::array<double,2*X_ELEMENTS> &rho, std::array<double,2*X_ELEMENTS> &stream);
void progressBar(const int current_step);
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

    if(time_step < 10)
        number = "000" + std::to_string(time_step);
    else if(time_step < 100)
        number = "00" + std::to_string(time_step);
    else if(time_step < 1000)
        number = "0" + std::to_string(time_step);
    else
        number = std::to_string(time_step);

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

void setupSystem(std::array<double,2*X_ELEMENTS> &rho_old, std::array<double,2*X_ELEMENTS> &stream) {
    for(int i = 0; i < 2*X_ELEMENTS; i++) {
        rho_old[i] = iniDist((i - X_ELEMENTS)*DX);
    }
    std::system((PY_SCRIPT + outputResults(rho_old,0) + " 0").c_str());     // Call the plotting script at t=0
    for(int i = 0; i < (int)STREAM_VEL*DT/DX; i++) {
        stream[i] = STREAM_CONC;
    }
}

/* Update the system a single time step */
void runDiffusion(std::array<double,2*X_ELEMENTS> &rho, std::array<double,2*X_ELEMENTS> &rho_old) {
    for(int x = 0; x < 2*X_ELEMENTS - 1; x++) {
        rho[x] = rho_old[x] + ALPHA*(rho_old[x-1] - 2*rho_old[x] + rho_old[x+1]);
    }
    rho[0] = rho[1], rho[2*X_ELEMENTS - 1] = rho[2*X_ELEMENTS - 2]; // Update boundaries
}

/* Model delivery of vesicles on cytoplasmic stream. They come from the left (-x) */
void cytoplasmicStream(std::array<double,2*X_ELEMENTS> &rho, std::array<double,2*X_ELEMENTS> &stream) {
    double dist_moved = STREAM_VEL*DT;
    int x_segs_travelled = (int)(dist_moved/DX);
    stream_seg_remainder_G += dist_moved/DX - x_segs_travelled;
    if(stream_seg_remainder_G >= 1) {
        x_segs_travelled++;
        stream_seg_remainder_G -= 1;
    }
    for(int i = 0; i < 2*X_ELEMENTS - x_segs_travelled; i++) {
        stream[i + x_segs_travelled] = stream[i];
    }
    for(int i = 0; i < x_segs_travelled; i++) {
        stream[i] = STREAM_CONC;
    }

    int delivery_zone = (int)(DELIVERY_RADIUS/DX);
    for(int i = X_ELEMENTS - delivery_zone; i < X_ELEMENTS + delivery_zone; i++) {
        rho[i] += stream[i];
        stream[i] = 0.0;

        double conc_above_max = rho[i] - MAXIMUM_CONC;
        if(conc_above_max > 0.0) {
            rho[i] -= conc_above_max;
            stream[i] += conc_above_max;
        }
    }
}

/* Despite the stream being bound (and not diffusing) it will still fluoresce) */
std::array<double, 2*X_ELEMENTS> visibleConc(std::array<double,2*X_ELEMENTS> &rho, std::array<double,2*X_ELEMENTS> &stream) {
    std::array<double,2*X_ELEMENTS> vis_conc = {0};
    for(int i = 0; i < 2*X_ELEMENTS; i++) {
        vis_conc[i] = rho[i] + stream[i];
    }
    return vis_conc;
}

void progressBar(const int current_step) {
        int bar_width = 70;
        double progress = (double)current_step/(double)T_ELEMENTS;

        std::cout << "[";
        int pos = (int)(bar_width * progress);
        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
        }
        std::cout << "] " << int(progress * 100.0) << " %\r";
        std::cout.flush();
}

void diffusionSolver() {

    /* Set up starting conditions, t = 0 */
    std::array<double,2*X_ELEMENTS> rho_old = {0}, rho = {0}, stream = {0};
    setupSystem(rho_old, stream);

    for(int t = 1; t <= T_ELEMENTS; t++) {
        runDiffusion(rho, rho_old);

        /* Output the results to a file and call the plotting script only after a certain number of steps */
        if(t % STEPS_BETWEEN_PLOTS == 0) {
            std::system((PY_SCRIPT + outputResults(visibleConc(rho, stream),t/STEPS_BETWEEN_PLOTS) + " " + std::to_string(t*DT)).c_str());
        }
        cytoplasmicStream(rho, stream);
        std::copy(rho.begin(), rho.end(), rho_old.begin());

        progressBar(t);
    }
    std::cout << std::endl;
}