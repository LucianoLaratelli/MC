#include "MonteCarlo.h"

int main(int argc, char *argv[])
{
    GCMC_System sys;
    
    int n;
    
    MoveType move_type;

    double currentPE,
           newPE;

    if(argc < 5 || argc > 8)
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
            printf("It also takes three (optional) flags:\n"
                   "\t-ideal     : simulates an ideal gas\n"
                   "\t-energy    : outputs energy to a file\n"
                   "\t-pos       : outputs positions to a file\n");
            exit(EXIT_FAILURE);
    }

    int arg_count = 1;

    sys.ideal_flag=false;
    sys.energy_output_flag = false;
    sys.output_flag = false;
    sys.stockmayer_flag = false;
    
    //take flags if specified
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
        else if (strcmp(argv[i],"-output")==0)
        {
            sys.output_flag = true;
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
    if(sys.output_flag)
    {
        sys.output = fopen("output.txt", "w");
    }

    sys.unweightedradial = fopen("unweightedradialdistribution.txt", "w");
    sys.weightedradial = fopen("weightedradialdistribution.txt", "w");


    sys.start_time = clock();
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("|                      STARTING  GCMC                      |\n");
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    sys.step = 0;
    sys.cutoff = sys.box_side_length * .5;

    particle added;

    added.x[0] = 0;
    added.x[1] = 0;
    added.x[2] = 0;
    added.dipole[0] = 1/85.10597636;
    added.dipole[1] = 0/85.10597636;
    added.dipole[2] = 0/85.10597636;
    sys.particles.push_back(added);

    particle added2;
    added2.x[0] = 4;
    added2.x[1] = 0;
    added2.x[2] = 0;
    added2.dipole[0] = 1/85.10597636;
    added2.dipole[1] = 0/85.10597636;
    added2.dipole[2] = 0/85.10597636;
    sys.particles.push_back(added2);

    currentPE = calculate_PE(&sys);//energy at first step 

    if(sys.energy_output_flag)
    {
        fprintf(sys.energies, "0 %lf\n", currentPE);
    }


    n = sys.particles.size(); //particle count 
    sys.sumenergy = currentPE;
    sys.sumparticles = n;
    sys.volume = sys.box_side_length * sys.box_side_length * sys.box_side_length;

    for(sys.step = 1; sys.step<sys.maxStep; sys.step++)
    {
            move_type = make_move(&sys); 
            
            newPE = calculate_PE(&sys);
            if(sys.step % (sys.maxStep/10) == 0)
            {
                double cycles_till_now = (double)(clock()-sys.start_time),
                       time_till_now = cycles_till_now/CLOCKS_PER_SEC;
                printf("  %.0f%% of iteration steps done. Time elapsed:"\
                        " %.2lf seconds.\n",\
                        ((double)sys.step/(double)sys.maxStep)*100, time_till_now);
            }
            if(move_accepted(currentPE, newPE,\
                        move_type, &sys))
            {
                    currentPE = newPE;//updates energy
                    sys.sumenergy += newPE;
                    output(&sys,newPE);
                    if(sys.step>=sys.maxStep*.5)
                    {
                        n = sys.particles.size();
                        sys.sumparticles += n;
                        radialDistribution(&sys,sys.step);
                    }
            }
            else // Move rejected
            {
                    undo_move(&sys, move_type);
                    sys.sumenergy += currentPE;
                    output(&sys,currentPE);
                    if(sys.step>=sys.maxStep*.5)
                    {
                        n = sys.particles.size();
                        sys.sumparticles += n;
                        radialDistribution(&sys,sys.step);
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
    if(sys.output_flag)
    {
        fclose(sys.output);
    }
    fclose(sys.unweightedradial);
    fclose(sys.weightedradial);

    return 0;
}
