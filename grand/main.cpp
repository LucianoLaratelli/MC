#include "MonteCarlo.h"

int main(int argc, char *argv[])
{
    GCMC_System sys;
    
    int step,
        n;
    
    MoveType move_type;

    double currentPE,
           newPE;

    if(argc < 5)
    {
            printf("This program takes four arguments:\n the type of "\
                    "particle, the desired number of iterations,"\
                    "the length of one side of the box,"\
                    " and the desired temperature, in that order.\n");
            printf("It also takes one of two options:\n"\
                   "-ideal: simulates an ideal gas exactly by setting pe = 0\n"\
                   "-energy: outputs energy to file\n");
            exit(EXIT_FAILURE);
    }

    int count =1;

    sys.ideal_flag=false;
    sys.energy_output_flag=false;
    
    //take flags if specified
    //no, you can't use both at the same time
    //if that doesn't make sense, think about what an ideal gas is and what
    //calculate_pe is calculating
    for(int i = 1;i < 2;i++)
    {

        if(strcmp(argv[i], "-ideal")==0)//makes calculate__pe return 0
        {
            sys.ideal_flag = true;
            count++;
            break;
        }
        else if(strcmp(argv[i],"-energy")==0)//output energy
        {
            sys.energy_output_flag = true;
            count++;
            break;
        }
        else
        {
            break;
        }
    }

    sscanf(argv[count], "%s", sys.particle_type);                     count++;
    sscanf(argv[count], "%d", &sys.maxStep);/*number of iterations*/  count++;
    sscanf(argv[count], "%lf", &sys.box_side_length);                 count++;
    sscanf(argv[count], "%lf", &sys.system_temp);//kelvin

    printf("particle type = %s\nmax step = %d\nbox side = %.2lf\ntemp = %.2lf\n",sys.particle_type,sys.maxStep,sys.box_side_length,sys.system_temp);

    sys.nBins = sys.box_side_length/sys.BinSize;
    sys.boxes = (double*)(calloc(sys.nBins,sizeof(double)));

    input(&sys);//set particle type

    sys.sigma_squared = sys.sigma*sys.sigma;
    sys.sigma_sixth = sys.sigma_squared * sys.sigma_squared * sys.sigma_squared;
    sys.sigma_twelfth = sys.sigma_sixth * sys.sigma_squared;

    srandom(time(NULL));//seed for random is current time

    sys.positions = fopen("positions.xyz", "w");
    if(sys.energy_output_flag)
    {
        sys.energies = fopen("energies.dat", "w");
    }
    sys.unweightedradial = fopen("unweightedradialdistribution.txt", "w");
    sys.weightedradial = fopen("weightedradialdistribution.txt", "w");
    sys.particlecount = fopen("particlecount.dat", "w");

    if(sys.positions == NULL|| \
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
    sys.volume = sys.box_side_length * sys.box_side_length * sys.box_side_length;
    sys.cutoff = sys.box_side_length * .5;

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
                    sys.sumenergy += newPE;
                    output(&sys,newPE,step);
                    if(step>=sys.maxStep*.5)
                    {
                        n = sys.particles.size();
                        sys.sumparticles += n;
                        radialDistribution(&sys,step);
                    }
            }
            else // Move rejected
            {
                    undo_move(&sys, move_type);
                    sys.sumenergy += currentPE;
                    output(&sys,currentPE,step);
                    if(step>=sys.maxStep*.5)
                    {
                        n = sys.particles.size();
                        sys.sumparticles += n;
                        radialDistribution(&sys,step);
                    }
            }
    }
    double cycles_till_now = (double)(clock()-sys.start_time),
           time_till_now = cycles_till_now/CLOCKS_PER_SEC;
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("|                      GCMC  COMPLETE                      |\n");
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("  This run took %f seconds.\nHave a nice day!\n"\
            ,time_till_now);//always good to have manners

    free(sys.boxes);
    fclose(sys.positions);
    if(sys.energy_output_flag)
    {
        fclose(sys.energies);
    }
    fclose(sys.unweightedradial);
    fclose(sys.weightedradial);
    fclose(sys.particlecount);

    return 0;
}
