#include <vector>
#include <math.h>
#include <stdio.h>
#include <random>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

typedef struct particle_particle
{
	double x[3];
} particle;

typedef struct _translational_data
{
	int pick;
	double phi, gamma, delta;
} translational_data;

typedef struct _removal_data
{
	double phi, gamma, delta;
} removal_data;

const int box_side_length = 100;

typedef struct _GCMC_System
{
        FILE * positions;
        FILE * energies;
        //FILE * unweightedradial;
        FILE * weightedradial;
        FILE * particlecount;
	std::vector <particle> particles;
	translational_data move;
	removal_data destroy;
        //Lennard-Jones parameters
	double epsilon,
               particle_mass,
               sigma,
               sigma_squared,//powers of sigma for PE
               sigma_sixth,
               sigma_twelfth;
        char particle_type[25];
        //system variables
        double system_temp,
               cutoff;
        int    maxStep,
               volume;
        //for averaging
        double sumparticles,
               sumenergy;
        //next three lines are for radial distribution function
        static const int nBins = 400; 
        double BinSize = .25; 
        double boxes[nBins] = {0};
        clock_t start_time;
} GCMC_System;

enum MoveType { TRANSLATE, CREATE_PARTICLE, DESTROY_PARTICLE };

const double k = 1.0; //boltzmann constant
const double h = 6.626e-34;//planck constant
const double conv_factor = 0.0073389366;//converts ATM to K/A^3

void input(GCMC_System *sys);


double calculate_PE(GCMC_System *sys);
double distfinder(GCMC_System *sys, int id_a, int id_b);

double randomish();
MoveType make_move(GCMC_System *sys);

void create_particle(GCMC_System *sys);
void move_particle(GCMC_System *sys, int pick);
void destroy_particle(GCMC_System *sys, int pick);

bool move_accepted(double cpe, double npe, MoveType move_type,\
                   GCMC_System *sys);

void undo_move(GCMC_System *sys, MoveType move);
void undo_insertion(GCMC_System *sys);
void unmove_particle(GCMC_System *sys);

double sphere_volume(GCMC_System *sys,double diameter);
void radialDistribution( GCMC_System *sys,int step);

void output(GCMC_System *sys,double accepted_energy , int step);

#endif
