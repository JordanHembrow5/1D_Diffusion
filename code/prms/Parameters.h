//
// Created by Jordan Hembrow on 2019-01-19.
//
#include <cmath>
#include <iostream>

#define DIFFUSION_CONST 0.003

#define STREAM_VEL 0.3
#define STREAM_CONC 0.02
#define DELIVERY_RADIUS 1.5
#define MAXIMUM_CONC (1.0 - STREAM_CONC)
#define STREAM_DELIVERY_RATE ((1.0/6.0)*STREAM_CONC)

#define X_SIZE 10
#define X_ELEMENTS 100
#define DX ((double)X_SIZE/(double)X_ELEMENTS)

#define T_SIZE 600
#define T_ELEMENTS 600
#define DT ((double)T_SIZE/(double)T_ELEMENTS)

/* Dimensionless constant that is used in finite element analysis (Forward Euler) */
#define ALPHA (DIFFUSION_CONST*DT/pow(DX,2))

/* Test stability of alpha. Results are meaningless if it is unstable */
const bool alpha_is_unstable = (ALPHA >= 0.5);
#if (alpha_is_unstable)
    #error Alpha unstable. Adjust parameters to ensure it is < 0.5
#endif

/* Probability of delivery from cytoplasmic stream */
#define GAUSSIAN_AMP 1.0
#define GAUSSIAN_MEAN 0.0

#define PY_SCRIPT "python3 ../code/DiffusionPlot.py "
#define TIMELAPSE "ffmpeg -f image2 -r 10 -i img/diffusion_1D_%04d.png -vcodec mpeg4 -b:v 25M -y Diffusion.mp4 -hide_banner -loglevel panic"

#define STEPS_BETWEEN_PLOTS 30

/* Initial distribution of particles */
#define iniDist(x) ((MAXIMUM_CONC/5.0)*(abs(x) > 8.0))  // Just bleached