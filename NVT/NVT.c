/****************************************************************************
BASIC NVT MONTE CARLO
*****************************************************************************  
HOW TO USE THIS PROGRAM:
************************
You need a .txt file with initial coordinates for your atoms. If you got this
from me or from github, you can use startgenerator to do this.
That program only works well with particle numbers with integer cube roots,
so this program uses similar values of N.
This program asks for one input from the user: the number of times
it is to "try" - that is, move a particle p some random displacement and 
see what happens to the PE.

-Luciano Laratelli
***************************************************************************/

/***************************************************************************
CURRENT ISSUES:

Energy is way too high.
A fundamental misunderstanding of the monte carlo method.
***************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

struct particle
{
    double x[3];
};

struct move_values
{
    int p;
    double delta, gamma, zeta;
};

const int N = 2; //number of particles
const int T = 300; //kelvin
const int L = 5; //length of one side of the cube, L = sigma

struct particle particles[N];
struct move_values move;

/************************
distfinder finds the distance between any two distinct particles of
number "id_a" and "id_b" respectively. this distance is used to
calculate the potential energy in PEfinder using the lennard-jones potential 
************************/
double distfinder(int id_a, int id_b)
{
    int i;
    double dist, dist2,delta;
    dist2 = 0; //clearing junk from memory just in case
    for(i=0;i<3;i++)
    {
        delta = particles[id_a].x[i] - particles[id_b].x[i];
        if(fabs(delta) > 0.5 * L) //this is for boundary conditions
            {
                if(delta > 0)
                {
                    delta -= L;
                }
                if(delta < 0)
                {
                    delta += L;
                }
            }
        dist2 += delta * delta;
        //finds distance coordinate by coordinate,
        //keeping the sum of the squares in memory
    }
    dist = sqrt(dist2); 
    //distance formula says dist = sqrt of the sum of the squares
    return dist;
}

/**********************
PEfinder finds the total potential energy between every particle
by using the lennard-jones potential
values of one are used for sigma and epsilon 
**********************/
double PEfinder()
{
    double sigma, epsilon, r,pe;
    double s2,s4,s6,s12;//exponents of sigma
    double rinv,r2,r4,r6,r12;//exponents of r
    double sor6,sor12;
    int b,c;
    epsilon = 1; //this can be changed to experimental values
    sigma = 1; // same as above
    s2 = sigma * sigma;
    s6 = s2*s2*s2;
    s12 = s6*s6;
    for(b=0;b<N;b++)
    {
        for(c=b+1;c<N;c++)
        {
            r = distfinder(b,c);
            rinv = 1/r;
            //formula has sigma/r; this is a bit easier to work with
            r2 = rinv * rinv;
            r6 = r2*r2*r2;
            r12 = r6* r6;
            sor6 = r6 * s6;
            sor12 = r12 * s12;
            //makes it harder to mess up the formula, h/t adam for the tip
            pe += (4 * epsilon * (sor12-sor6));
        }
    }
    return pe;
}

/**********************
Positionchecker is called by rand_p_mover.
it looks at every particle p in the system; if rand_p moves p out of the box,
positionchecker makes it so that the particle "pops in" from the other side
of the box. think pacman or snake.
**********************/
void positionchecker(int id_a)
{
    int i;
    for(i=0;i<3;i++)
    {
        if(particles[id_a].x[i] > L)
        {
            particles[id_a].x[i] -= L;
        }
        if(particles[id_a].x[i] < 0)
        {
            particles[id_a].x[i] += L;
        }
    }
    return;
}

/*********************
this next function is the first step of the monte carlo method; 
it moves a random particle a random distance away from its current position
and then calls positionchecker to enforce boundary conditions
*********************/
void rand_p_mover()
{
    int i,p;
    double delta,gamma,zeta;
    p = rand() / (RAND_MAX / N + 1); //p is a random int between 0 and N
    delta = (rand()/(RAND_MAX)) - 0.5; //a random float between -0.5 and 0.5
    gamma = (rand()/(RAND_MAX)) - 0.5; //same as above
    zeta = (rand()/(RAND_MAX)) - 0.5; //same as above
    move.p = p;
    move.delta = delta;
    move.gamma = gamma;
    move.zeta = zeta;
    particles[p].x[0] += delta;
    particles[p].x[1] += gamma;
    particles[p].x[2] += zeta;
    //we can't leave the box (nor do we want to be at the surface;)
    //if p manages to escape the box, positionchecker() bosses it around 
    positionchecker(p);
    return;
}

/****************
rand_p_unmover does something very important: if a move is rejected, it
makes the system go back to the state it was in before the move!
it uses the struct move_values to do this. 
****************/


void rand_p_unmover()
{
    int p;
    p = move.p;
    particles[p].x[0] -= move.delta;
    particles[p].x[1] -= move.gamma;
    particles[p].x[1] -= move.zeta;
    return;
}


/********************
this is the second part of the monte carlo method; it takes the PE from PE finder
for the current configuration and subtracts it from the previous configuration. 
it substitutes delta E into the e^-beta equation and compares the value of
THAT exponential to a random number between 0 and 1
if the result of the exponential is greater than the random number, the change
is accepted
********************/
bool E_checker(double cpe, double npe)
{
    double deltaE,guess,k,beta,prob; //no ints in this function, no sir!
    deltaE = cpe - npe; // between the "current" PE and the "new" PE
    guess=((double)rand()/(double)RAND_MAX);
    k = 1; //Boltzmann constant in the natural units for the system
    beta = 1/(k*T); //temperature as set above
    prob = exp(-beta*deltaE);
    if(prob > guess)
    {
        return true;
    }
    else
    {
        return false;
    }
}

/*******************
starting_positions takes a .txt file with coordinates 
for the starting positions of the particles in this simulation.
each line of the file is an (x,y,z) with format x y z\n
*******************/
/*
void starting_positions()
{
    int p;
    FILE * startingpositions; //any .txt file
    startingpositions = fopen("startingpositions.txt", "r");
    for(p=0;p<N;p++)
    {
        fscanf(startingpositions,"%lf %lf %lf\n", &particles[p].x[0],&particles[p].x[1],&particles[p].x[2]);
    }
    fclose(startingpositions);
    return;
}*/

/******************
output_to_file lives up to its name;
It ouputs the coordinates of the particles at every frame into an .xyz file
*******************/
void output_to_file()
{
    FILE * positions;
    positions = fopen("positions.xyz","a");
    int p;
    for(p=0;p<N;p++)
    {
        fprintf(positions,"H %lf %lf %lf\n", particles[p].x[0], particles[p].x[1], particles[p].x[2]);
        //hoping this works lmao
        //update: it works!
    }
    fclose(positions);
    return;
}

int main()
{
    FILE * positions; 
    positions = fopen("positions.xyz","w");
    fclose(positions);
    double cpe, npe;
    FILE * energies; 
    energies = fopen("energies.dat","w");
    bool guess; 
    int c; //count
    int m; //maximum number of tries
    printf("How many tries do you want to do?\n"); //user-directed!
    scanf("%d", &m);
    //starting_positions();
    particles[0].x[0] = 0;
    particles[0].x[1] = 0;
    particles[0].x[2] = 0;
    particles[1].x[0] = 2;
    particles[1].x[1] = 2;
    particles[1].x[2] = 2;
    fprintf(positions, "%d \n\n",N);
    cpe = PEfinder(); //"current potential energy"
    for(c=0;c<m;c++)
    {
        rand_p_mover(); //monte carlo step one
        npe = PEfinder(); //"new potential"
        guess = E_checker(cpe,npe); //sets guess to bool, step two of monte carlo
        if (guess == true) //if new energy is less than the random number, accept
        {
            output_to_file(); //for just the positions
            cpe = npe; //updates energy
            fprintf(energies,"%d %f\n", c, cpe); //the new one
            // if the new energy is accepted, it becomes the 
            // "current" energy for the next loop 
        }
        else //hopefully self-explanatory
        {
            //have to actually reject the move!
            rand_p_unmover(); //returns the system to the previous state
            continue;
        }
       
    }
    fclose(energies);
    
    printf("Done! Hope it worked out.\n"); //it never does
    return 0;
}
