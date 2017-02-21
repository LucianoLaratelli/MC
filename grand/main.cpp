#include "MonteCarlo.h"

int main(int argc, char *argv[])
{
    GCMC_System sys;
    
    int step,
        n;
    
    MoveType move_type;

    double currentPE,
           newPE;

    if(argc < 5 || argc > 9)
    {
            printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
                   "|               INPUT ERROR               |\n"
                   "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
            printf("This program takes at least four arguments,\n"
                    "in the following order:\n"
                    "\tthe type of particle,\n"
                    "\tthe desired number of iterations,\n"\
                    "\tthe length of one side of the box,\n"\
                    "\tand the desired temperature.\n");
            printf("It also takes four (optional) flags:\n"
                   "\t-ideal     : simulates an ideal gas\n"
                   "\t-energy    : outputs energy to a file\n"
                   "\t-pos       : outputs positions to a file\n"
                   "\t-stockmayer: simulates a stockmayer fluid\n");
            exit(EXIT_FAILURE);
    }

    int arg_count = 1;

    sys.ideal_flag=false;
    sys.energy_output_flag=false;
    
    //take flags if specified
    //no, you can't use both at the same time
    //if that doesn't make sense, think about what an ideal gas is and what
    //calculate_pe is calculating
    for(int i = 1;i < argc;i++)
    {

        if(strcmp(argv[i], "-ideal")==0)//makes calculate__pe return 0
        {
            sys.ideal_flag = true;
            arg_count++;
            break;
        }
        else if(strcmp(argv[i],"-energy")==0)//output energy
        {
            sys.energy_output_flag = true;
            arg_count++;
            break;
        }
        else if (strcmp(argv[i],"-stockmayer")==0)
        {
            sys.stockmayer_flag = true;
            arg_count++;
            break;
        }
        else if (strcmp(argv[i],"-pos")==0)
        {
            sys.positions_output_flag = true;
            arg_count++;
            break;
        }
        else
        {
            break;
        }
    }

    sscanf(argv[arg_count], "%s", sys.particle_type);               arg_count++;
    sscanf(argv[arg_count], "%d", &sys.maxStep);/*iteration number*/arg_count++;
    sscanf(argv[arg_count], "%lf", &sys.box_side_length);           arg_count++;
    sscanf(argv[arg_count], "%lf", &sys.system_temp);//kelvin


    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("|                     SYSTEM VARIABLES                     |\n");
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("                   PARTICLE TYPE     = %s                   \n"
           "                   NUM OF ITERATIONS = %d                   \n"
           "                   BOX SIDE LENGTH   = %.0lf                \n"
           "                   TEMPERATURE       = %.0lf                \n",  
           sys.particle_type,sys.maxStep,sys.box_side_length,sys.system_temp);
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    sys.nBins = sys.box_side_length/sys.BinSize;
    sys.boxes = (double*)(calloc(sys.nBins,sizeof(double)));

    input(&sys);//set particle type

    sys.polarizability = 2;
    sys.sigma_squared = sys.sigma*sys.sigma;
    sys.sigma_sixth = sys.sigma_squared * sys.sigma_squared * sys.sigma_squared;
    sys.sigma_twelfth = sys.sigma_sixth * sys.sigma_squared;

    srandom(time(NULL));//seed for random is current time

    if(sys.energy_output_flag)
    {
        sys.energies = fopen("energies.dat", "w");
    }
    if(sys.positions_output_flag)
    {
        sys.positions = fopen("positions.xyz", "w");
    }

    sys.unweightedradial = fopen("unweightedradialdistribution.txt", "w");
    sys.weightedradial = fopen("weightedradialdistribution.txt", "w");


    sys.start_time = clock();
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("|                      STARTING  GCMC                      |\n");
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    currentPE = calculate_PE(&sys);//energy at first step 

    if(sys.energy_output_flag)
    {
        fprintf(sys.energies, "0 %lf\n", currentPE);
    }


    n = sys.particles.size(); //particle count 
    sys.sumenergy = currentPE;
    sys.sumparticles = n;
    sys.volume = sys.box_side_length * sys.box_side_length * sys.box_side_length;
    sys.cutoff = sys.box_side_length * .5;
    particle added;

    added.x[0] = 0;
    added.x[1] = 0;
    added.x[2] = 0;
    sys.particles.push_back(added);

    particle added2;
    added2.x[0] = 4;
    added2.x[1] = 0;
    added2.x[2] = 0;
    sys.particles.push_back(added2);


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
                    output(&sys,newPE);
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
                    output(&sys,currentPE);
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
    printf("This run took %f seconds.\nHave a nice day!\n"\
            ,time_till_now);//always good to have manners

    free(sys.boxes);

    
    if(sys.energy_output_flag)
    {
        fclose(sys.energies);
    }
    if(sys.energy_output_flag)
    {
        fclose(sys.positions);
    }
    fclose(sys.unweightedradial);
    fclose(sys.weightedradial);

    return 0;
}
