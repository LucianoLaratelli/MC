#include "MonteCarlo.h"


// if any particle tries to escape the box, it disappears
// and is replaced on the opposite side of the box by a more obedient particle
bool positionchecker(GCMC_System *sys, int particleID )
{
	for (int i = 0; i<3; i++)
	{
		if (sys->particles[particleID].x[i] >= L)
		{
			sys->particles[particleID].x[i] -= L;
			return false;
		}
		else if (sys->particles[particleID].x[i] < 0)
		{
			sys->particles[particleID].x[i] += L;
			return false;
		}
	}
	return true;
}


//randomish returns a random float between 0 and 1
double randomish()
{
	double r = rand();
	return r / ((double)RAND_MAX + 1);
}

void unmove_particle( GCMC_System *sys )
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

//particlemover moves a random particle a random distance
void move_particle(GCMC_System *sys, int pick)
{
	//phi,gamma, and delta are random floats between 0 and L
	bool i = false;
	while( !i )
	{
		double phi = (randomish())*L,
			gamma = (randomish())*L,
			delta = (randomish())*L;
		//these next few lines store our moves in a struct case
                //we need to undo the move
		sys->move.pick = pick;
		sys->move.phi = phi;
		sys->move.gamma = gamma;
		sys->move.delta = delta;
		//we make the moves (finally!)
		sys->particles[pick].x[0] += phi;
		sys->particles[pick].x[1] += gamma;
		sys->particles[pick].x[2] += delta;
		//lastly, we check that the particle hasn't moved out of the box 
		i = positionchecker(sys, pick);
		if(i == false)
		{
			unmove_particle(sys);
		}
	}
	return;
}


void create_particle( GCMC_System *sys)
{
	particle added; //making a struct to push_back to the vector
	added.x[0] = randomish()*L;
	added.x[1] = randomish()*L;
	added.x[2] = randomish()*L;
	sys->creator.phi = added.x[0];
	sys->creator.gamma = added.x[1];
	sys->creator.delta = added.x[2];
	sys->particles.push_back(added);
        //the index number of the new particle:
	sys->creator.pick = sys->particles.size() - 1; 
	return;
}


void destroy_particle(GCMC_System *sys, int pick)
{
	sys->destroy.pick = pick;
	sys->destroy.phi   = sys->particles[pick].x[0];
	sys->destroy.gamma = sys->particles[pick].x[1];
	sys->destroy.delta = sys->particles[pick].x[2];
        //"begin" (below) points to the address of the zeroth item 
	// we add "pick" amount of bytes 
	// to get to the next pointer
	sys->particles.erase(sys->particles.begin() + pick);
        return;
}


MoveType make_move( GCMC_System *sys)
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
			choice = randomish() * 3; //random float between 0 and 3
		fflush(stdout);
		if (choice<1.0)
		{
			move_particle(sys, pick);
			move = TRANSLATE;
		}
		else if (choice >= 2.0)
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

double distfinder(GCMC_System *sys, int id_a, int id_b)
{
	double dist,
		delta,
		delta2 = 0.00;
	for (int i = 0; i<3; i++)
	{
		delta = sys->particles[id_a].x[i] - sys->particles[id_b].x[i];
		if (delta >= cutoff)
		{
			delta -= L;
		}
		if (delta <= -cutoff)
		{
			delta += L;
		}
		delta2 += delta * delta;
	}
	dist = sqrt(delta2);
	//printf("%lf \n", dist);
	return dist;
}

double calculate_PE(GCMC_System *sys)//returns energy in Kelvin
{
	int b, //counter variables
		c,
		pool = sys->particles.size();
	double r,
		pe = 0.00;
	double s2,//powers of sigma
		s6,
		s12;
	double rinv,//inverse powers of distance between two particles
		r2,
		r6,
		r12;
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
			rinv = 1 / r;
			r2 = rinv * rinv;
			r6 = r2 * r2 * r2;
			r12 = r6 * r6;
			sor6 = s6 * r6;
			sor12 = s12 * r12;
			pe += 4 * sys->epsilon * (sor12 - sor6);
		}
	}
	return pe;
}

/*******************************************************************************
* based on page 130 of UMS(2.ed.) by Frenkel and Smit.                 
*******************************************************************************/
bool move_accepted(double cpe, double npe, int c, MoveType move_type,\
                   int n, double particle_mass, double system_temp)
{
	FILE * energies;
	energies = fopen("energies.dat", "a");
        double delta = npe - cpe,
               random = randomish(),
               pi = M_PI,
               beta = 1 / (k * system_temp),//thermodynamic beta
               // lambda is the de broglie thermal wavelength
               lambda = (h) / (sqrt(2 * pi*particle_mass*k*system_temp)),
               lambdacubed = lambda * lambda * lambda,
               volume = L * L * L,
               particle_density = n / L, //  (n * 1.88) / L,
               //mu is the chemical potential, NEEDS FIX
               mu = k * system_temp * log(lambdacubed) * particle_density;
	int pool = n, //size of the system BEFORE the move we are making
            poolplus = pool + 1;//size of the system AFTER particle insertion 
	if (delta < 0)//always accept a move that lowers the energy
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
		if (acceptance> random)
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
	double radius = diameter / 2;//this is important
	return 4.0 / 3.0 * M_PI * radius * radius * radius;
}

void radialDistribution(GCMC_System *sys, int n)
{
	FILE * weightedradial;
	FILE * unweightedradial;
	unweightedradial = fopen("unweightedradialdistribution.txt", "a");
	weightedradial = fopen("weightedradialdistribution.txt", "a");
	const int nBins = (L / 2 + 1) * 10; //total number of bins
	double BinSize = L / (2.0 * (double)nBins),
		num_density = (double)n / (double)(L*L*L),
		expected_number_of_particles,
		boxes[nBins] = { 0 },
		current_shell,
		previous_shell,
		shell_volume_delta,
		dist;
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

void output(GCMC_System *sys, char *particle_type)
{
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
