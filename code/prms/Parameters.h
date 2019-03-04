//
// Created by Jordan Hembrow on 2019-01-19.
//
#include <cmath>
#include <iostream>

#define DIFFUSION_CONST 0.003

#define STREAM_VEL 0.3
#define STREAM_CONC 0.05
#define DELIVERY_RADIUS 1.5
#define MAXIMUM_CONC (1.0 - STREAM_CONC)
#define STREAM_DELIVERY_RATE 0.003      /* conc per second delivered */

#define X_SIZE 10.0
#define X_ELEMENTS 100
#define DX ((double)X_SIZE/(double)X_ELEMENTS)

#define T_SIZE 600.0
#define T_ELEMENTS 600
#define DT ((double)T_SIZE/(double)T_ELEMENTS)

/* Dimensionless constant that is used in finite element analysis (Forward Euler) */
#define ALPHA (DIFFUSION_CONST*DT/pow(DX,2))

/* Probability of delivery from cytoplasmic stream */
#define GAUSSIAN_AMP 1.0
#define GAUSSIAN_MEAN 0.0

#define PLOT_ALL "python3 ../code/DiffusionPlotAll.py data/ "
#define TIMELAPSE "ffmpeg -f image2 -r 10 -i img/diffusion_1D_%04d.png -vcodec mpeg4 -b:v 2M -y Diffusion.mp4 -hide_banner -loglevel panic"

#define STEPS_BETWEEN_PLOTS 3

/* Initial distribution of particles */
#define iniDist(x) ((MAXIMUM_CONC/5.0)*(abs(x) > 5.0))  // Just bleached