// TODO: Split main.cpp into two separate files
// TODO: Put some physical units and values in

#include <iostream>
#include <cmath>
#include <array>
#include <fstream>
#include <algorithm>
#include <random>
#include "prms/Parameters.h"


/* Global variables (yuck) are suffixed with 'G', don't mess with them */
double stream_seg_remainder_G = 0.0;


void testStability();
std::string outputResults(const std::array<double,2*X_ELEMENTS> &conc, const int time_step);
void runDiffusion(std::array<double,2*X_ELEMENTS> &rho, std::array<double,2*X_ELEMENTS> &rho_old);
double gaussian(double x, double mean, double stddev, double amp);
double RNG(double min, double max);
bool streamDelivers(int current_pos_element);
void cytoplasmicStream(std::array<double,2*X_ELEMENTS> &rho, std::array<double,2*X_ELEMENTS> &stream);
std::array<double, 2*X_ELEMENTS> visibleConc(std::array<double,2*X_ELEMENTS> &rho, std::array<double,2*X_ELEMENTS> &stream);
void progressBar(const int current_step);
void diffusionSolver();
void plotResults();


int main() {
    testStability();
    std::system("rm data/*.txt\nrm img/*.png");     // Clear previous data so that ffmpeg doesn't get confused
    diffusionSolver();
    plotResults();
    return 0;
}

/* Test stability of alpha. Results are meaningless if it is unstable */
void testStability() {
    const bool alpha_is_unstable = (bool)(ALPHA >= 0.5);
    std::string input = "N";
    if(alpha_is_unstable) {
        std::cerr << "Error! Alpha is unstable! (" << ALPHA << ")" << std::endl;
        exit(2);
    }
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
    outputResults(rho_old,0);
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

/* returns the values of f(x) for the function 'f' being a gaussian */
double gaussian(double x, double mean, double stddev, double amp) {
    double prefactor = amp/sqrt(2.0*M_PI*pow(stddev,2));
    double exponential = exp(-pow(x - mean, 2)/(2*pow(stddev,2)));
    return prefactor*exponential;
}

/* Returns a random number in the range [min, max) */
double RNG(double min, double max) {
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 number(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> distribution(min, max);
    return distribution(number);
}

bool streamDelivers(int current_pos_element) {
    double position = DX*(current_pos_element - X_ELEMENTS);    // simulates around x=0 so total elements is 2*X_ELEMENTS
    double random = RNG(0.0, GAUSSIAN_AMP);
    double delivery_prob = gaussian(position, GAUSSIAN_MEAN, DELIVERY_RADIUS, GAUSSIAN_AMP);
    if(delivery_prob > random) {
        return true;
    }
    return false;
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

    /* Looping in the opposite direction of the stream so to ensure that no new data gets interpreted as old */
    for(int i = 2*X_ELEMENTS - (x_segs_travelled + 1); i >= 0; i--) {
        stream[i + x_segs_travelled] = stream[i];
    }
    for(int i = 0; i < x_segs_travelled; i++) {
        stream[i] = STREAM_CONC;
    }

    double stream_delivery_vol = STREAM_DELIVERY_RATE * DT;
    for(int i = 0; i < 2*X_ELEMENTS; i++) {
        if(streamDelivers(i)) {
            if (stream[i] < STREAM_DELIVERY_RATE) {
                rho[i] += stream[i];
                stream[i] = 0.0;
            } else {
                rho[i] += STREAM_DELIVERY_RATE;
                stream[i] -= STREAM_DELIVERY_RATE;
            }

            double conc_above_max = rho[i] - MAXIMUM_CONC;
            if (conc_above_max > 0.0) {
                rho[i] -= conc_above_max;
                stream[i] += conc_above_max;
            }
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
    std::cout << "Running simulation..." << std::endl;
    std::array<double,2*X_ELEMENTS> rho_old = {0}, rho = {0}, stream = {0};
    setupSystem(rho_old, stream);

    for(int t = 1; t <= T_ELEMENTS; t++) {
        runDiffusion(rho, rho_old);

        /* Output the results to a file and call the plotting script only after a certain number of steps */
        if(t % STEPS_BETWEEN_PLOTS == 0) {
            outputResults(visibleConc(rho, stream),t/STEPS_BETWEEN_PLOTS);
        }
        cytoplasmicStream(rho, stream);
        std::copy(rho.begin(), rho.end(), rho_old.begin());

        progressBar(t);
    }
}

void plotResults() {
    std::cout << std::endl << "Plotting results..." << std::endl;
    double time_between_plots = DT*STEPS_BETWEEN_PLOTS;
    int files_to_plot = (T_ELEMENTS/STEPS_BETWEEN_PLOTS) + 1;
    std::string plotting = PLOT_ALL + std::to_string(time_between_plots) + " " + std::to_string(files_to_plot);
    std::system(plotting.c_str());
    std::system(TIMELAPSE);
}