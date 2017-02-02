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

typedef struct _GCMC_System {
	std::vector <particle> particles;
	translational_data move;
	creation_data creator;
	removal_data destroy;
	double sigma,
	       epsilon,
               particle_mass;
        double system_temp,
               cutoff;
        int    maxStep;
} GCMC_System;


enum MoveType { TRANSLATE, CREATE_PARTICLE, DESTROY_PARTICLE };

const int box_side_length = 22;
const double k = 1.38065 * pow(10,-23); //boltzmann constant
const double h = 6.626 * pow(10, -34);//planck constant

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
bool move_accepted(double cpe, double npe, int c, MoveType move_type,\
                   int n, GCMC_System *sys);
double qst_calc(int N, double energy, int c, double system_temp);
double sphere_volume(double diameter);
void radialDistribution( GCMC_System *sys, int n);
void input(GCMC_System *sys, char *particle_type);
void output(GCMC_System *sys, char *particle_type);

#endif
