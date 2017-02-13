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
    //sys holds variables related to the system:
    //particle number, stored moves, LJ sigma and epsilon
    GCMC_System sys;
    
    //for input
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

    srandom(1);//time(NULL));
    FILE * positions;
    positions = fopen("positions.xyz", "w");
    fclose(positions);
    FILE * energies;
    energies = fopen("energies.dat", "w");//we'll close this later
    FILE * unweightedradial;
    unweightedradial = fopen("unweightedradialdistribution.txt", "w");
    fclose(unweightedradial);
    FILE * weightedradial;
    weightedradial = fopen("weightedradialdistribution.txt", "w");
    fclose(weightedradial);
    FILE * chemicalpotential;
    chemicalpotential = fopen("chemicalpotential.txt","w");
    fclose(chemicalpotential);

    currentPE = calculate_PE(&sys);
    fprintf(energies, "0 %lf\n", currentPE);
    fclose(energies);

    n = sys.particles.size();  // get particle count
    sumenergy = currentPE;//the sum has to include initial starting energy
    sumparticles = n;

    int largest_number_of_particles = n;

    for(step = 1; step<sys.maxStep; step++)
    {
            //first step of the monte carlo algorithm
            move_type = make_move(&sys); 
            
            newPE = calculate_PE(&sys);

            if( move_accepted(currentPE, newPE,\
                        move_type, n, &sys))
            {
                    n = sys.particles.size();
                    output(&sys,newPE,step);
                    currentPE = newPE;//updates the current energy
                    sumenergy += currentPE;//adds to the running total
                    sumparticles += n;
                    //qst_calc(sumparticles, sumenergy, step, sys.system_temp);
                    radialDistribution(&sys, n,step);
                    if(n>largest_number_of_particles)
                    {
                        largest_number_of_particles = n;
                        printf("Highest number of particles = %d\n",\
                                largest_number_of_particles);
                    }
            }
            else // Move rejected
            {
                    n = sys.particles.size();
                    undo_move(&sys, move_type);
                    output(&sys,currentPE,step);
                    sumenergy += currentPE;
                    sumparticles += n;
                    //qst_calc(sumparticles, sumenergy, step, sys.system_temp);
                    radialDistribution(&sys, n,step);
                    if(n>largest_number_of_particles)
                    {
                        largest_number_of_particles = n;
                        printf("Highest number of particles = %d\n",\
                                largest_number_of_particles);
                    }
            }
    }
    printf("Run complete! Have a nice day.\n");
    return 0;
}
