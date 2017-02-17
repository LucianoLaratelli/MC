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
    sys.sigma_squared = sys.sigma*sys.sigma;
    sys.sigma_sixth = sys.sigma_squared * sys.sigma_squared * sys.sigma_squared;
    sys.sigma_twelfth = sys.sigma_sixth * sys.sigma_squared;

    srandom(time(NULL));//seed for random is current time

    sys.volume = box_side_length * box_side_length * box_side_length;

    sys.positions = fopen("positions.xyz", "w");
    sys.energies = fopen("energies.dat", "w");
    sys.unweightedradial = fopen("unweightedradialdistribution.txt", "w");
    sys.weightedradial = fopen("weightedradialdistribution.txt", "w");
    sys.particlecount = fopen("particlecount.dat", "w");

    if(sys.positions == NULL|| \
       sys.energies == NULL|| \
       sys.unweightedradial == NULL|| \
       sys.weightedradial == NULL|| \
       sys.particlecount == NULL) 
    {
        printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        printf("Error! One of the file pointers in main() is broken!\n");
        printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
        exit(EXIT_FAILURE);
    }

    sys.start_time = clock();
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("|                      STARTING  GCMC                      |\n");
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    currentPE = calculate_PE(&sys);//energy at first step 
    fprintf(sys.energies, "0 %lf\n", currentPE);

    n = sys.particles.size(); //particle count 
    sys.sumenergy = currentPE;
    sys.sumparticles = n;
    sys.cutoff = box_side_length * .5;

    for(step = 1; step<sys.maxStep; step++)
    {
            move_type = make_move(&sys); 
            
            newPE = calculate_PE(&sys);
            if(step % (sys.maxStep/10) == 0)
            {
                double cycles_till_now = (double)(clock()-sys.start_time),
                       time_till_now = cycles_till_now/CLOCKS_PER_SEC;
                printf("  %.0f%% of iteration steps done. Time elapsed:"\
                        " %.2lf seconds.\n",\
                        ((double)step/(double)sys.maxStep)*100, time_till_now);
            }
            if(move_accepted(currentPE, newPE,\
                        move_type, &sys))
            {
                    currentPE = newPE;//updates energy
                    n = sys.particles.size();
                    sys.sumenergy += newPE;
                    sys.sumparticles += n;
                    //output(&sys,newPE,step);
                    radialDistribution(&sys,step);
            }
            else // Move rejected
            {
                    undo_move(&sys, move_type);
                    n = sys.particles.size();
                    sys.sumenergy += currentPE;
                    sys.sumparticles += n;
                    //output(&sys,currentPE,step);
                    radialDistribution(&sys,step);
                    
            }
    }
    double cycles_till_now = (double)(clock()-sys.start_time),
           time_till_now = cycles_till_now/CLOCKS_PER_SEC;
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("|                      GCMC  COMPLETE                      |\n");
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("  This run took %f seconds.\nHave a nice day!\n"\
            ,time_till_now);//always good to have manners

    fclose(sys.positions);
    fclose(sys.energies);
    fclose(sys.unweightedradial);
    fclose(sys.weightedradial);
    fclose(sys.particlecount);

    return 0;
}
