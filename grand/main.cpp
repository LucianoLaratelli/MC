#include "MonteCarlo.h"

int main(int argc, char *argv[])
{
    GCMC_System sys;
    
    int step,
        n;

    MoveType move_type;
    double currentPE,
           newPE;

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

    input(&sys);//set particle type
    srandom(time(NULL));//seed for random is current time

    sys.positions = fopen("positions.xyz", "w");
    sys.energies = fopen("energies.dat", "w");
    sys.unweightedradial = fopen("unweightedradialdistribution.txt", "w");
    sys.weightedradial = fopen("weightedradialdistribution.txt", "w");
    sys.particlecount = fopen("particlecount.dat", "w");
    sys.average_energies = fopen("average_energy.dat","w");

    if(sys.positions == NULL|| \
       sys.energies == NULL|| \
       sys.unweightedradial == NULL|| \
       sys.weightedradial == NULL|| \
       sys.particlecount == NULL|| \
       sys.average_energies == NULL)
    {
        printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        printf("Error! One of the file pointers in main() is broken!\n");
        printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        exit(EXIT_FAILURE);
    }

    sys.start_time = clock();
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("BEGINNING MONTE CARLO RUN\n");
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    currentPE = calculate_PE(&sys);//energy at first step 
    fprintf(sys.energies, "0 %lf\n", currentPE);

    n = sys.particles.size(); //particle count 

    sys.sumenergy = currentPE;
    sys.sumparticles = n;

    for(step = 1; step<sys.maxStep; step++)
    {
            move_type = make_move(&sys); 
            
            newPE = calculate_PE(&sys);
            if( move_accepted(currentPE, newPE,\
                        move_type, &sys))
            {
                    n = sys.particles.size();
                    output(&sys,newPE,step);
                    currentPE = newPE;//updates energy
                    sys.sumenergy += currentPE;
                    sys.sumparticles += n;
                    if(step > step*0.5)//allow system to equilibrate before g(r)
                    {
                        radialDistribution(&sys,step);
                    }
            }
            else // Move rejected
            {
                    undo_move(&sys, move_type);
                    n = sys.particles.size();
                    output(&sys,currentPE,step);
                    sys.sumenergy += currentPE;
                    sys.sumparticles += n;
                    if(step > step*0.5)
                    {
                        radialDistribution(&sys,step);
                    }
            }
    }
    double cycles_till_now = (double)(clock()-sys.start_time),
           time_till_now = cycles_till_now/CLOCKS_PER_SEC;

    printf("This run took %lf seconds.\n Have a nice day!\n"\
            ,time_till_now);//always good to have manners

    fclose(sys.positions);
    fclose(sys.energies);
    fclose(sys.unweightedradial);
    fclose(sys.weightedradial);
    fclose(sys.particlecount);
    fclose(sys.average_energies);

    return 0;
}
