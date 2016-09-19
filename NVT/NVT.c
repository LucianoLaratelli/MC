#include <stdio.H>
#include <math.H>
#include <stdlib.H>

struct particle
{
    double x[3];
};

double sigma, epsilon, m;

#define N 5 //number of particles
struct particle particles[N];

/************************
distfinder finds the distance between any two distinct particles of
number "id_a" and "id_b" respectively. this distance is used to
calculate the potential energy in PEfinder using the lennard-jones
potential, so it's important
************************/
double distfinder(int id_a, int id_b)
{

    return dist
}

/**********************
PEfinder finds the total potential energy between every particle
by using the lennard-jones potential
values of one are used for sigma and epsilon due to laziness but those
can be changed if you need a better model
**********************/
double PEfinder()
{

    return pe;
}

/*********************
this next function is part of the monte carlo method; it moves a random 
particle a random distance away
*********************/
void rand_p_mover()
{

    return;
}

/********************
this is the second part of the monte carlo method; it takes the PE from PE finderfor the current configuration and subtracts it from the previous configuration. 
it substitutes that delta E into the e^-beta equation and compares the value of
THAT exponential to a random number between 0 and 1
if the result of the exponential is less than the random number, it returns 1
if the result of the exponential is greater than the random number, it returns 0
********************/
bool E_checker()
{

    return True;
    return False;
}

void position_receiver()
{

    return;
}

int main()
{
    return 0;
}
