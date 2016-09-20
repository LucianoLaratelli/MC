#include <stdio.H>
#include <math.H>
#include <stdlib.H>
#include <stdbool.H>

struct particle
{
    double x[3];
};

double sigma, epsilon, m;

#define N 5 //number of particles
struct particle particles[N];
/******************************************************************************
TO DO:
1. Figure out how starting_positions will work ???????????????? maybe done???????
2. Translate the math for E_checker into code
3. Work out how much the random displacement should be for rand_p_mover
4. Write code for output to file
5. ???
6. profit
******************************************************************************/



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
    epsilon = 1; //this can be changed to experimental values
    sigma = 1; / same as above
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
                for(i=0;i<3;i++)
                {
                    pe += 0.5 * (4 * epsilon * ((sor12)-sor6));
                    // each sum does half the PE per pair
                    // because the way this loop is written 
                    // has each pair included twice
                    // typing a "0.5 *" is easier than thinking about loops
                }
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
    FILE * startingpositions;//any .txt file
    startingpositions = fopen("startingpositions.txt", "r");
    for(p=0;p<N;p++)
    {
        fscanf("%d %d %d", &particles[p].x[0],&particles[p].x[1],&particles[p].x[2])
        //not sure this will work, TEST!!! TEST!!! TEST!!!
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
    return;
}
int main()
{
    FILE * positions;//this is the .xyz file
    FILE * energies; //this is the .txt file
    bool guess; 
    int c;//count
    int m;//maximum number of tries
    starting_positions();
    for(c=0;c<m;c++)
    {
        pe = PEfinder();
        rand_p_mover();
        guess = E_checker();
        if (guess == 1)
        {
            output_to_file;
        }
        else
        {
            continue;
        }
    }
    return 0;
}
