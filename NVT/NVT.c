#include <stdio.H>
#include <math.H>
#include <stdlib.H>
#include <stdbool.H>

struct particle
{
    double x[3];
};

double sigma, epsilon, m;
/****************************************************************************
HOW TO USE THIS PROGRAM
You need a .txt file with initial coordinates for your atoms. If you got this
from me or from github, you can use startgenerator to do this.
That program works best with integer particle numbers with integer cube roots,
so this program uses similar values of N.
This program asks for one input from the user: the number of times
it is to "try" - that is, move a particle p some random displacement and 
see what happens to the PE.

-Luciano Laratelli
***************************************************************************/
#define N 27//number of particles
#define T 300//kelvin
#define L 20//length of one side of the cube, 20L = 20 atomic radii I guess

struct particle particles[N];

/************************
distfinder finds the distance between any two distinct particles of
number "id_a" and "id_b" respectively. this distance is used to
calculate the potential energy in PEfinder using the lennard-jones potential 
************************/
double distfinder(int id_a, int id_b)
{
    int i;
    double dist, dist2;
    dist2 = 0;//clearing junk from memory
    for(i=0;i<3;i++)
    {
        dist2 += pow(particles[id_a].x[i]-particles[id_b].x[i],2);
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
values of one are used for sigma and epsilon due to laziness but those
can be changed if you need a better model
**********************/
double PEfinder()
{
    double sigma, epsilon, r;
    double s2,s4,s6,s12;//exponents of sigma
    double rinv,r2,r4,r6,r12;//exponents of r
    double sor6,sor12;
    int b,c,pe;
    epsilon = 1; //this can be changed to experimental values
    sigma = 1; // same as above
    s2 = sigma * sigma;
    s6 = s2*s2*s2;
    s12 = s6*s6;
    for(b=0;b<N;b++)
    {
        for(c=0;c<N;c++)
        {
            if(b!=c)
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
                pe += 0.5 * (4 * epsilon * (sor12-sor6));
                // each sum does half the PE per pair
                // because the way this loop is written 
                // has each pair included twice
                // typing a "0.5 *" is easier than thinking about loops
            }
        }
    }           
    return pe;
}

/*********************
this next function is part of the monte carlo method; it moves a random 
particle a random distance away from its current position
*********************/
void rand_p_mover()
{
    int i,p;
    double delta
    p = rand() / (RAND_MAX / N + 1);//p is a random int between 0 and N
    i = rand() / (RAND_MAX / 4);//i is a random int between 0 and 3
    delta = rand()/(RAND_MAX+ 1.0);//delta is a random int between 0 and 1
    //we can't leave the box (nor do we want to be at the surface,)
    //so p tries to escape V
    //we multiply by a new random int between 0 and 1 to correct it
    //this is probably the wrong way to do it
    particles[p].x[i] += delta;
    if(particles[p].x[i] >= L) 
    { 
        delta = rand()/(RAND_MAX + 1.0);
        particles[p].x[i] =  particles[p].x[i]* delta;
    }
    return;
}

/********************
this is the second part of the monte carlo method; it takes the PE from PE finderfor the current configuration and subtracts it from the previous configuration. 
it substitutes that delta E into the e^-beta equation and compares the value of
THAT exponential to a random number between 0 and 1
if the result of the exponential is less than the random number, it returns 1
if the result of the exponential is greater than the random number, it returns 0
********************/
bool E_checker(double cpe, double npe)
{
    double deltaE,guess,k,beta,prob;
    deltaE = cpe - npe; //delta between the "current" PE and the "new" PE
    guess=((double)rand()/(double)RAND_MAX);
    k = 1.3806485279 * pow(10,-23);//Boltzmann constant in joules per kelvin
    beta = 1/(k*T);//temperature as defined above
    prob = exp(-1*beta*deltaE);
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
for the starting positions of the atoms in this simulation.
The first line of the file is the x coordinate of particle p
the second line of the file is the y coordinate of particle p
the third line of the file is the z coordinate of particle p
etc.
*******************/
void starting_positions()
{
    int p;
    FILE * startingpositions;//any .txt file
    startingpositions = fopen("startingpositions.txt", "r");
    for(p=0;p<N;p++)
    {
        fscanf(startingpositions,"%lf %lf %lf\n", &particles[p].x[0],&particles[p].x[1],&particles[p].x[2]);
    }
    return;
}

/******************
output_to_file lives up to its name. It ouputs two things:
the coordinates of the particles at every frame into an .xyz file;
the total energy of the system at every frame
*******************/
void output_to_file()
{
    FILE * positions;//this is the .xyz file
    positions = fopen("positions.xyz","a");
    int p,n;
    for(p=0;p<N;p++)
    {
        fprintf(positions,"H %lf %lf %lf\n", particles[p].x[0], particles[p].x[1], particles[p].x[2]);
        //hoping this works lmao
    }
    return;
}

int main()
{
    FILE * positions;//this is the .xyz file
    positions = fopen("positions.xyz","w");
    double cpe, npe;
    FILE * energies; //this is the .txt file
    energies = fopen("energies.dat","w");
    bool guess; 
    int c;//count
    int m;//maximum number of tries
    printf("How many tries do you want to do?\n");//user-directed!
    scanf("%d", &m);
    starting_positions();
    fprintf(positions, "%d \n\n",N);
    cpe = PEfinder();//"current potential energy"
    for(c=0;c<m;c++)
    {
        cpe = cpe;
        rand_p_mover();//monte carlo step one
        npe = PEfinder();//"new potential"
        guess = E_checker(cpe,npe);//sets guess to bool, step two of monte carlo
        if (guess == true)//if new energy is less than the random number
        {
            output_to_file();//does positions
            cpe = npe;
            fprintf(energies,"%d %f\n", c, cpe);
            // if the new energy is accepted, it becomes the 
            // "current" energy for the next loop 
        }
        else //hopefully self-explanatory
        {
            continue;
            cpe = cpe;
        }
        if(c == (.1*m))
        {
            printf("10 percent of the way there!");
        }
        if(c == (.25*m))
        {
            printf("25 percent of the way there!");
        }
        if(c == (.5*m))
        {
            printf("50 percent of the way there!");
        }
        if(c == (.75*m))
        {
            printf("75 percent of the way there!");
        }
    }
    fclose(energies);
    fclose(positions);
    printf("Done! Hope it worked out.\n"); //it never does
    return 0;
}
