#include "MonteCarlo.h"

// if any particle tries to escape the box, it disappears
// and is replaced on the opposite side of the box
// by a more obedient particle
bool positionchecker(GCMC_System *sys, int particleID )
{
	for (int i = 0; i<3; i++)
	{
		if (sys->particles[particleID].x[i] >= box_side_length)
		{
			sys->particles[particleID].x[i] -= box_side_length;
			return false;
		}
		else if (sys->particles[particleID].x[i] < 0)
		{
			sys->particles[particleID].x[i] += box_side_length;
			return false;
		}
	}
	return true;
}


//returns a random float between 0 and 1
double randomish()
{
	double r = rand();
	return r / ((double)RAND_MAX + 1);
}

//undoes a rejected displacement
void unmove_particle(GCMC_System *sys )
{
	int pick = sys->move.pick;
	double phi = sys->move.phi,
		gamma = sys->move.gamma,
		delta = sys->move.delta;
	sys->particles[pick].x[0] -= phi;
	sys->particles[pick].x[1] -= gamma;
	sys->particles[pick].x[2] -= delta;
	return;
}

//moves a random particle a random distance
void move_particle(GCMC_System *sys, int pick)
{
	bool i = false;
	while( !i )
	{
                //phi,gamma, and delta are random floats between 0 and L
		double phi = (randomish())* box_side_length,
			gamma = (randomish())*box_side_length,
			delta = (randomish())*box_side_length;
                // store move in case it needs to be undone
		sys->move.pick = pick;
		sys->move.phi = phi;
		sys->move.gamma = gamma;
		sys->move.delta = delta;
		//make the moves (finally!)
		sys->particles[pick].x[0] += phi;
		sys->particles[pick].x[1] += gamma;
		sys->particles[pick].x[2] += delta;
		//check that the particle hasn't moved out of the box 
		i = positionchecker(sys, pick);
		if(i == false)
		{
			unmove_particle(sys);
		}
	}
	return;
}


//insert particle at random location
void create_particle( GCMC_System *sys)
{
	particle added; //making a struct to push_back to the vector
	added.x[0] = randomish()*box_side_length;
	added.x[1] = randomish()*box_side_length;
	added.x[2] = randomish()*box_side_length;
        //store coordinates of inserted particle in case of rejection
	sys->creator.phi = added.x[0];
	sys->creator.gamma = added.x[1];
	sys->creator.delta = added.x[2];
	sys->particles.push_back(added);
	sys->creator.pick = sys->particles.size() - 1; //index of new particle
	return;
}


//particle deletion
void destroy_particle(GCMC_System *sys, int pick)
{
	sys->destroy.pick = pick;
	sys->destroy.phi   = sys->particles[pick].x[0];
	sys->destroy.gamma = sys->particles[pick].x[1];
	sys->destroy.delta = sys->particles[pick].x[2];
        //"begin" (below) points to the address of the zeroth item 
	// we add "pick" amount of bytes 
	// to get to the address of the "pickth" item 
	sys->particles.erase(sys->particles.begin() + pick);
        return;
}


//make a random move based on a random number. The more random the better!
MoveType make_move(GCMC_System *sys)
{
	MoveType move;
	int pool = sys->particles.size();

	if (pool == 0)//if there are 0 particles we must add some to get started
	{
		create_particle(sys);
		pool = sys->particles.size();
		move = CREATE_PARTICLE;
	}
	else
	{
		double pick = rand() % pool,//picks random particle
		       choice = randomish(); //random float between 0 and 1
		fflush(stdout);
		if (choice<(1.0/3.0))
		{
			move_particle(sys, pick);
			move = TRANSLATE;
		}
		else if (choice >= (2.0/3.0))
		{
			create_particle(sys);
			pool = sys->particles.size();
			move = CREATE_PARTICLE;
		}
		else
		{
			destroy_particle(sys, pick);
			move = DESTROY_PARTICLE;
		}
	}
	return move;
}

//undoes rejected moves based on move types
void undo_move(GCMC_System *sys, MoveType move)
{
	if (move == TRANSLATE)
	{
		unmove_particle(sys);
	}
	else if (move == CREATE_PARTICLE)
	{
		destroy_particle( sys, sys->creator.pick);
	}
	else
	{
		particle added;
		added.x[0] = sys->creator.phi;
		added.x[1] = sys->creator.gamma;
		added.x[2] = sys->creator.delta;
		sys->particles.push_back(added);
	}
	return;
}

//finds distance between every pair of particles for use in 
//potential calculation
double distfinder(GCMC_System *sys, int id_a, int id_b)
{
	double dist,
               delta,
               delta2 = 0.00,
               cutoff = box_side_length * 0.5;
        //distance is done coordinate by coordinate
        //(remember x[0] = x position etc
	for(int i = 0; i<3; i++)
	{
		delta = sys->particles[id_a].x[i] - sys->particles[id_b].x[i];
                //the following conditionals are the minimum image convention
                //that make up part of periodic boundary conditions
		if (delta >= cutoff)
		{
			delta -= box_side_length;
		}
		else if (delta <= -cutoff)
		{
			delta += box_side_length;
		}
                //distance requires sum of squares of distance in each
                //coordinate direction
		delta2 += delta * delta;
	}
	dist = sqrt(delta2);
	return dist;
}

//calculates potential of the system
double calculate_PE(GCMC_System *sys)//returns energy in Kelvin
{
	int b, //counter variables
            c,
            pool = sys->particles.size();
	double r,
	       pe = 0.00,
               cutoff = box_side_length * 0.5;
        //powers of sigma
	double s2,
               s6,
               s12;
        //inverse powers of distance between two particles
	double rinv,
               rinv2,
               rinv6,
               rinv12;
	double sor6,//powers of sigma over r for the potential equation
               sor12;
	s2 = sys->sigma * sys->sigma;
	s6 = s2 * s2 * s2;
	s12 = s6 * s6;
	for (b = 0; b < pool - 1; b++)
	{
            for (c = b + 1; c < pool; c++)
            {
                r = distfinder(sys, b, c);
                if (r<cutoff)
                {
                        continue;
                }
                rinv = 1.0 / r;
                rinv2 = rinv * rinv;
                rinv6 = rinv2 * rinv2 * rinv2;
                rinv12 = rinv6 * rinv6;
                sor6 = s6 * rinv6;
                sor12 = s12 * rinv12;
                pe += 4.0 * sys->epsilon * (sor12 - sor6);
            }
	}
	return pe;
}

/*******************************************************************************
* based on page 130 of UMS(2.ed.) by Frenkel and Smit.                 
*******************************************************************************/
bool move_accepted(double cpe, double npe, int c, MoveType move_type,\
                   int n, GCMC_System *sys)
{
	FILE * energies;
        FILE * chemicalpotential;
	energies = fopen("energies.dat", "a");
        chemicalpotential = fopen("chemicalpotential.txt", "a");
        double delta = npe - cpe,
               random = randomish(),
               pi = M_PI,
               beta = 1 / (k * sys->system_temp),//thermodynamic beta
               // lambda is the de broglie thermal wavelength
               lambda = (h) / (sqrt(2 * pi*sys->particle_mass*k*\
                           sys->system_temp)),
               lambdacubed = lambda * lambda *  lambda,
               volume = box_side_length * box_side_length *\
                        box_side_length,
               particle_density = n / box_side_length, 
               //mu is the chemical potential, NEEDS FIX
               mu = k * sys->system_temp * log(lambdacubed) * particle_density;
	int pool = n, //size of the system BEFORE the move we are making
            poolplus = pool + 1;//size of the system AFTER particle insertion 
        fprintf(chemicalpotential, "%d %f \n",c, mu);
        fclose(chemicalpotential);
        //always accept a move that lowers the energy:
	if (delta < 0)
	{
		fprintf(energies, "%d %f \n", c, npe);
		fclose(energies);
		return true;
	}
	else if (move_type == TRANSLATE )//if we MOVED a particle 
	{
		double acceptance = exp(-1 * beta * delta);
		if (acceptance > random)
		{
			fprintf(energies, "%d %f \n", c, npe);
			fclose(energies);
			return true;
		}
		else
		{
			fprintf(energies, "%d %f \n", c, npe);
			fclose(energies);
			return false;
		}
	}
	else if (move_type == CREATE_PARTICLE)//if we CREATED a particle
	{
		double volume_term = volume / (lambdacubed*poolplus),
                       goes_in_exponential = beta* (mu - npe + cpe),
                       acceptance = volume_term*exp(goes_in_exponential);
		if (acceptance > random)
		{
			fprintf(energies, "%d %f \n", c, npe);
			fclose(energies);
			return true;
		}
		else
		{
			fprintf(energies, "%d %f \n", c, npe);
			fclose(energies);
			return false;
		}
	}
	else if (move_type == DESTROY_PARTICLE )//if we DESTROYED a particle
	{
		double volume_term = (lambdacubed*pool) / volume,
			goes_in_exponential = -beta * (mu + delta),
			acceptance = volume_term * exp(goes_in_exponential);
		if (acceptance > random)
		{
			fprintf(energies, "%d %f \n", c, npe);
			fclose(energies);
			return true;
		}
		else
		{
			fprintf(energies, "%d %f \n", c, npe);
			fclose(energies);
			return false;
		}
	}
	return true;
}

double qst_calc(int N, double energy, int c, double system_temp)
{
	FILE * qsts;
	qsts = fopen("qsts.dat", "a");
	double T = system_temp,
               average_N = N / c,//<N>
               average_N_all_squared = average_N * average_N,//<N>^2
               average_energy = energy / c,//<U>
               N_squared = average_N * average_N,
               average_of_N_squared = N_squared * N_squared,//<N^2>
               particles_by_energy = average_N * average_energy,
               average_of_particles_by_energy = particles_by_energy / c, //<NU>
               numerator = (average_of_particles_by_energy)-\
                           (average_N * average_energy),
               denominator = average_of_N_squared - average_N_all_squared,
               QST = k * T - (numerator / denominator);
	fprintf(qsts, "%d %lf\n", c, QST);
	fclose(qsts);
	return QST;
}

double sphere_volume(double diameter)
{
	double radius = diameter / 2.0;//this is important
	return (4.0 / 3.0) * M_PI * radius * radius * radius;
}

void radialDistribution(GCMC_System *sys, int n)
{
	FILE * weightedradial;
	FILE * unweightedradial;
	unweightedradial = fopen("unweightedradialdistribution.txt", "a");
	weightedradial = fopen("weightedradialdistribution.txt", "a");
	const int nBins = (box_side_length / 2 + 1) * 10; //total number of bins
	double BinSize = box_side_length / (2.0 * (double)nBins),
		num_density = (double)n / (double)(box_side_length* \
                               box_side_length*box_side_length),
		expected_number_of_particles,
		boxes[nBins] = {0},
		current_shell,
		previous_shell,
		shell_volume_delta,
		dist,
                cutoff = box_side_length * 0.5;
	int IK;
	for (int I = 0; I<n - 1; I++)
	{
		for (int K = I + 1; K<n; K++)
		{
			dist = distfinder(sys, I, K);
			if (dist < cutoff)
			{
				IK = int(dist / BinSize);
				boxes[IK] += 2;
			}
		}
	}
	printf("number of particles: %d\n", n);
	previous_shell = 0;
	for (int I = 1; I <= nBins; I++)
	{
		fprintf(unweightedradial, "%lf\n", boxes[I - 1]);
		current_shell = I;
		shell_volume_delta = (sphere_volume(current_shell) -\
                                      sphere_volume(previous_shell));
		expected_number_of_particles = shell_volume_delta * num_density;
		boxes[I - 1] /= expected_number_of_particles;
		fprintf(weightedradial, "%lf\n", boxes[I - 1]);
		previous_shell = current_shell;
	}
	fclose(unweightedradial);
	fclose(weightedradial);
	return;
}

/*******************************************************************************
 * input reads in the particle type and stores its corresponding mass (AMU),
 * LJ epsilon (K), and LJ sigma(angstroms) in the system struct
 * ****************************************************************************/
void input(GCMC_System *sys, char *particle_type)
{
   char argon[] = "Ar",
        helium[] = "He",
        neon[] = "Ne",
        krypton[] = "Kr",
        xenon[] = "Xe",
        water[] = "Water";
   if(strcmp(particle_type,argon) == 0)
   {
       sys->sigma = 3.371914;
       sys->epsilon = 128.326802;
       sys->particle_mass = 39.948;
   }
   else if(strcmp(particle_type,helium) == 0)
   {
       sys->sigma = 2.653089;
       sys->epsilon = 9.071224;
       sys->particle_mass = 4.0026;
   }
   else if(strcmp(particle_type,neon) == 0)
   {
       sys->sigma = 2.785823;
       sys->epsilon = 36.824138;
       sys->particle_mass = 20.1797;
   }
   else if(strcmp(particle_type,krypton) == 0)
   {
       sys->sigma = 3.601271;
       sys->epsilon = 183.795833;
       sys->particle_mass = 83.798;
   }
   else if(strcmp(particle_type,xenon) == 0)
   {
       sys->sigma = 3.956802;
       sys->epsilon = 237.985247;
       sys->particle_mass = 131.293;
   }
   //water is for Stockmeyer fluids only
   else if(strcmp(particle_type,water) == 0)
   {
       sys->sigma = 3.15100;
       sys->epsilon = 76.42000;
       sys->particle_mass = 18.016;
   }
   else
   {
       printf("Not a supported chemical species!\nAllowed values for Lennard-"\
               "Jones are:\nAr\nNe\nHe\nKr\nXe\nAllowed values for "\
               "Stockmeyer are:\nWater\nPlease try again!\n");
       exit(EXIT_FAILURE);
   }
   return;
}     


void output(GCMC_System *sys, char *particle_type)
{
        char water[] = "Water",
             oxygen[] = "O";
        if(strcmp(particle_type,water)==0)
        {
            strcpy(particle_type,oxygen);
        }
	FILE * positions;
	positions = fopen("positions.xyz", "a");
	int p,
            pool;
	pool = sys->particles.size();
	fprintf(positions, "%d \n\n", pool);
	for (p = 0; p<pool; p++)
	{
		fprintf(positions, "%s %lf %lf %lf\n",\
                        particle_type, sys->particles[p].x[0],\
                        sys->particles[p].x[1], sys->particles[p].x[2]);
	}
	fclose(positions);
	return;
}
