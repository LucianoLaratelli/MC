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




****************************************************************************/
#define N 5 //number of particles
#define T 300//kelvin
#define L 20//length of one side of the cube, 20L = 20 atomic radii I guess

struct particle particles[N];
/******************************************************************************
TO DO:
1. Figure out how starting_positions will work ???????????????? maybe done???????
5. Write code that outputs starting positions into 
   a text file for starting_positions
6. ???
7. profit
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
    //below, the two "random variables" are chosen so they are integers
    //between 0 and the box length. They are then divided to give a 
    //double, which (I think) allows for more states of the system
    for(p=0;p<N;p++)
    {
        for(i=0;i<3;i++)
        {
            int delta,gamma;
            delta = rand() % L;//random variable for a random value
            gamma = rand() % L;//see above
            double dog = delta/gamma;//"delta over gamma"
            particles[p].x[i] += dog;
            while(fabs(particles[p].x[i]) >= L) 
            //we can't leave the box (nor do we want to be at the surface,)
            //so if dog makes p go out of our bounds
            //we re-roll until we get an acceptable value
            {
                delta = rand() % L;
                gamma = rand() % L;
                double dog = delta/gamma;
                particles[p].x[i] += dog;
            }
        }
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
    double prob,beta,k,deltaE;
    deltaE = cpe - npe; //delta between the "current" PE and the "new" PE
    double guess=((double)rand()/(double)RAND_MAX);
    k = 1.3806485279 * pow(10,-23);//joules per kelvin
    beta = 1/(k*T);//temperature as defined above
    prob = exp(-1*beta*deltaE);
    if(prob < guess)
    {
        return True;
    }
    else
    {
        return False;
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
    FILE * startingpositions;//any .txt file
    startingpositions = fopen("startingpositions.txt", "r");
    for(p=0;p<N;p++)
    {
        fscanf("%d %d %d\n", &particles[p].x[0],&particles[p].x[1],&particles[p].x[2]);
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
    int p,n;
    for(p=0;p<N;p++)
    {
        fprintf(positions,"H %d %d %d\n", particles[p].x[0], particles[p].x[1], particles[p].x[2]);
        //hoping this works lmao
    }
    return;
}

int main()
{
    FILE * positions;//this is the .xyz file
    FILE * energies; //this is the .txt file
    positions = fopen("positions.xyz","w");
    energies = fopen("energies.txt","w");
    bool guess; 
    int c;//count
    int m;//maximum number of tries
    printf("How many tries do you want to do?");//user-directed!
    scanf("%d", &m);
    starting_positions();
    cpe = PEfinder();//"current potential energy"
    for(c=0;c<m;c++)
    {
        rand_p_mover();//monte carlo step one
        npe = PEfinder();//"new potential"
        guess = E_checker(cpe,npe);//sets guess to bool, step two of monte carlo
        if (guess == 1)//if new energy is less than the random number
        {
            output_to_file;//does positions and energy for this frame
            cpe = npe;
            fprintf(energies,"%d %d", c, cpe);
            // if the new energy is accepted, it becomes the 
            // "current" energy for the next loop 
        }
        else //hopefully self-explanatory
        {
            continue;
            cpe = cpe;
        }
    }
    return 0;
}
