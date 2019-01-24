//
// Created by Jordan Hembrow on 2019-01-19.
//
#include <cmath>


#define DIFFUSION_CONST 0.003

#define X_SIZE 10
#define X_ELEMENTS 100
#define DX ((double)X_SIZE/(double)X_ELEMENTS)

#define T_SIZE 1200
#define T_ELEMENTS 1000
#define DT ((double)T_SIZE/(double)T_ELEMENTS)

/* Dimensionless constant that is used in finite element analysis (Forward Euler) */
#define ALPHA (DIFFUSION_CONST*DT/pow(DX,2))

/* Initial distribution of particles */
//#define STDEV 0.5
//#define iniDist(x) (exp(-(x*x)/(2*STDEV*STDEV)))

#define PY_SCRIPT "python3 ../DiffusionPlot.py "
#define TIMELAPSE "ffmpeg -f image2 -r 10 -i img/diffusion_1D_%04d.png -vcodec mpeg4 -y Diffusion.mp4"

#define STEPS_BETWEEN_PLOTS 20

//#define iniDist(x) (int)(x>=0)                        //Step function at x=0
//#define iniDist(x) (int)(abs(x)<=1)                     //Top-hat function at for -1 < x < 1
#define iniDist(x) ((int)(abs(x)))%2