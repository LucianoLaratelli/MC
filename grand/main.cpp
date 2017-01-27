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
	GCMC_System sys;
	
	char particle_type[128] = { 0 };
	int step,
		n,
		maxStep;
	MoveType move_type;
	double particle_mass,
		system_temp,
		currentPE,
		newPE,
		sumenergy,
		sumparticles;


	if (argc != 7)
	{
		printf("This program takes five arguments: the type of\
                        particle, its mass, a LJ sigma, LJ epsilon,the\
                        temperature of the system, and the number of\
                        desired iterations, in that order.");
		exit(EXIT_FAILURE);
	}
	sscanf(argv[1], "%s", particle_type);
	sscanf(argv[2], "%lf", &particle_mass);
	sscanf(argv[3], "%lf", &sys.sigma);
	sscanf(argv[4], "%lf", &sys.epsilon);
	sscanf(argv[5], "%lf", &system_temp);
	sscanf(argv[6], "%d", &maxStep);


	srand(time(NULL));
	FILE * positions;
	FILE * energies;
	FILE * qsts;
	FILE * unweightedradial;
	FILE * weightedradial;
	FILE * stats;
	stats = fopen("stats.txt", "a");
	unweightedradial = fopen("unweightedradialdistribution.txt", "w");
	weightedradial = fopen("weightedradialdistribution.txt", "w");
	positions = fopen("positions.xyz", "w");
	energies = fopen("energies.dat", "w");
	qsts = fopen("qsts.dat", "w");
	fclose(weightedradial);
	fclose(unweightedradial);
	fclose(positions);
	fclose(qsts);





	currentPE = calculate_PE(&sys);
	fprintf(energies, "0 %f\n", currentPE);
	fclose(energies);

	clock_t begin = clock();
	n = sys.particles.size();  // get particle count
	sumenergy = currentPE;//the sum has to include initial starting energy
	sumparticles = n;

	//printf("Starting run now\n");
	for( step = 1; step<maxStep; step++)
	{
		//first step of the monte carlo algorithm
		move_type = make_move(&sys); 
		
		newPE = calculate_PE(&sys);

		if( move_accepted(currentPE, newPE, step,\
                            move_type, n, particle_mass, system_temp))
		{
			n = sys.particles.size();
			output(&sys, particle_type);
			currentPE = newPE;//updates the current energy
			sumenergy += currentPE;//adds to the running total
			sumparticles += n;//see above
			qst_calc(sumparticles, sumenergy, step, system_temp);
		}
		else // Move rejected
		{
			n = sys.particles.size();
			undo_move(&sys, move_type);
			output(&sys, particle_type);
			sumenergy += currentPE;
			sumparticles += n;
			qst_calc(sumparticles, sumenergy, step, system_temp);
		}
	}

	radialDistribution(&sys, n);
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;	
	fclose(stats);
	printf("Done! This run took %f seconds.\
                Have a nice day!\n", time_spent);
	return 0;
}
