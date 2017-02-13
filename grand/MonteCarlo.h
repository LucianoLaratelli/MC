#include <vector>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <random>
#include <stdlib.h>
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

typedef struct _creation_data
{
	int pick;
	double phi, gamma, delta;
} creation_data;

typedef struct _removal_data
{
	int pick;
	double phi, gamma, delta;
} removal_data;
const int box_side_length = 22;

typedef struct _GCMC_System
{
        FILE * positions;
        FILE * energies;
        FILE * unweightedradial;
        FILE * weightedradial;
        FILE * chemicalpotential;
	std::vector <particle> particles;
	translational_data move;
	creation_data creator;
	removal_data destroy;
        //lennard-jones parameters
	double sigma,
	       epsilon,
               particle_mass;
        char particle_type[25];
        //system variables
        double system_temp,
               cutoff;
        int    maxStep;
        //next three lines are for radial distribution function
        static constexpr int nBins = (box_side_length/2 + 1)*100;
        double BinSize = box_side_length/(double)nBins;
        double boxes[nBins] = {0};
} GCMC_System;


enum MoveType { TRANSLATE, CREATE_PARTICLE, DESTROY_PARTICLE };
const double k = 1.0; //boltzmann constant
const double h = 6.626e-34;//planck constant

bool positionchecker(GCMC_System *sys, int particleID);
double randomish();
void unmove_particle(GCMC_System *sys);
void move_particle(GCMC_System *sys, int pick);
void create_particle(GCMC_System *sys);
void destroy_particle(GCMC_System *sys, int pick);
MoveType make_move(GCMC_System *sys);
void undo_move(GCMC_System *sys, MoveType move);
double distfinder(GCMC_System *sys, int id_a, int id_b);
double calculate_PE(GCMC_System *sys);
bool move_accepted(double cpe, double npe, MoveType move_type,\
                   int n, GCMC_System *sys);
double qst_calc(int N, double energy, int c, double system_temp);
double sphere_volume(double diameter);
void radialDistribution( GCMC_System *sys, int n,int step);
void input(GCMC_System *sys);
void output(GCMC_System *sys,double accepted_energy , int step);

#endif
