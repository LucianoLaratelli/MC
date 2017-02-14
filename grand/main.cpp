/*******************************************************************************
This is a monte carlo method simulator of the grand canonical ensemble,
also known as the "mu-V-T" ensemble.  This program takes a system composed of
any number of particles and does one of three moves on it:

1. moves a random particle a random distance
2. adds a particle at a random positions
3. removes a random particle

It outputs: the number of particles per frame, the average number of particles,
and the standard deviation associated with the average
the energies (in the same manner as the number of particles)
the chemical potential (per frame)
a calculated QST per-frame and as an average over all frames
*******************************************************************************/


#include "MonteCarlo.h"


int main(int argc, char *argv[])
{
    //sys holds system variables and anything that needs to be global 
    GCMC_System sys;
    
    int step,
        n;
    MoveType move_type;
    double currentPE,
           newPE,
           sumenergy,
           sumparticles;

    if (argc != 4)
    {
            printf("This program takes three arguments:\n the type of "\
                    "particle, the desired number of iterations,"\
                    " and the desired temperature, in that order.\n");
            exit(EXIT_FAILURE);
    }
    sscanf(argv[1], "%s", sys.particle_type);
    sscanf(argv[2], "%d", &sys.maxStep);//number of iterations
    sscanf(argv[3], "%lf", &sys.system_temp);//kelvin

    input(&sys);

    srandom(42);//time(NULL));
    sys.positions = fopen("positions.xyz", "w");
    sys.energies = fopen("energies.dat", "w");//we'll close this later
    sys.unweightedradial = fopen("unweightedradialdistribution.txt", "w");
    sys.weightedradial = fopen("weightedradialdistribution.txt", "w");
    sys.particlecount = fopen("particlecount.txt", "w");

    currentPE = calculate_PE(&sys);//energy at step 0
    fprintf(sys.energies, "0 %lf\n", currentPE);

    n = sys.particles.size();  // get particle count
    sumenergy = currentPE;//the sum has to include initial starting energy
    sumparticles = n;

    for(step = 1; step<sys.maxStep; step++)
    {
            //first step of the monte carlo algorithm
            //printf("~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
            //printf("BEGINNING MONTE CARLO MOVE\n");
            move_type = make_move(&sys); 
            
            newPE = calculate_PE(&sys);
            //printf("\n\nMAIN NPE = %lf\n\n",newPE);
            if( move_accepted(currentPE, newPE,\
                        move_type, n, &sys,step))
            {
                    n = sys.particles.size();
                    output(&sys,newPE,step,n);
                    currentPE = newPE;//updates the current energy
                    sumenergy += currentPE;//adds to the running total
                    sumparticles += n;
                    if(step > step*0.5)//we let the system equilibrate a bit
                    {
                        radialDistribution(&sys, n,step);
                    }
            }
            else // Move rejected
            {
                    undo_move(&sys, move_type);
                    n = sys.particles.size();
                    output(&sys,currentPE,step,n);
                    sumenergy += currentPE;
                    sumparticles += n;
                    if(step > step*0.5)
                    {
                        radialDistribution(&sys, n,step);
                    }
            }
            //printf("LEAVING MONTE CARLO MOVE\n");
            //printf("~~~~~~~~~~~~~~~~~~~~~~~~\n");
    }
    printf("Run complete! Have a nice day.\n");
    fclose(sys.positions);
    fclose(sys.energies);
    fclose(sys.unweightedradial);
    fclose(sys.weightedradial);
    fclose(sys.particlecount);
    return 0;
}
