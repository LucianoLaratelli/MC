#include "MonteCarlo.h"

/*******************************************************************************
 * input reads in the particle type and stores its corresponding mass (AMU),
 * LJ epsilon (K), and LJ sigma(angstroms) in the system struct
 * ****************************************************************************/
void input(GCMC_System *sys)
{
   char argon[] = "Ar",
        helium[] = "He",
        neon[] = "Ne",
        krypton[] = "Kr",
        xenon[] = "Xe",
        water[] = "Water";
   char* particle_type= sys->particle_type;
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
       char oxygen[] = "O";
       strcpy(sys->particle_type,oxygen);
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

double calculate_PE(GCMC_System *sys)
{
        int pool = sys->particles.size();
	double dist,
	       pe = 0.00,
               inverse_distance,//inverse distance
	       s_over_d_6,//powers of sigma over dist
               s_over_d_12;
	for (int a = 0; a < pool - 1; a++)
	{
            for (int b = a + 1; b < pool; b++)
            {
                dist = distfinder(sys, a, b);//in one dimension
                inverse_distance = 1.0 / dist;
                double dist_squared = inverse_distance*inverse_distance,
                       dist_sixth = dist_squared*dist_squared*dist_squared,
                       dist_twelfth = dist_sixth*dist_squared;
                s_over_d_12 = sys->sigma_twelfth * dist_twelfth; 
                s_over_d_6 = sys->sigma_sixth * dist_sixth; 
                pe += 4.0 * sys->epsilon * (s_over_d_12 - s_over_d_6);
            }
	}
	return pe;//in KELVIN
}

double distfinder(GCMC_System *sys, int id_a, int id_b)
{
	double dist,
               delta,
               delta2 = 0.00;

	for(int i = 0; i<3; i++)
	{
		delta = sys->particles[id_a].x[i] - sys->particles[id_b].x[i];
                //the following conditionals are the minimum image convention
		if (delta >= sys->cutoff)
		{
			delta -= sys->box_side_length;
		}
		else if (delta <= -sys->cutoff)
		{
			delta += sys->box_side_length;
		}
                //distance requires sum of squares of distance in each
                //coordinate direction
		delta2 += delta * delta;
	}
	dist = sqrt(delta2);
	return dist;
}

//returns a random float between 0 and 1
double randomish()
{
	double r = random();
	return r / ((double)RAND_MAX + 1);
}

MoveType make_move(GCMC_System *sys)
{
        //MoveType is an enum in MonteCarlo.h 
	MoveType move;
	int pool = sys->particles.size();

	if (pool == 0)
	{
		create_particle(sys);
		pool = sys->particles.size();
		move = CREATE_PARTICLE;
	}
	else
	{
		double pick = random() % pool,//picks random particle
		       choice = randomish(); //random float between 0 and 1
		fflush(stdout);
		if (choice<0.33333)
		{
			create_particle(sys);
			move = CREATE_PARTICLE;
		}
		else if (choice >= (.666667))
		{
			move_particle(sys, pick);
			move = TRANSLATE;
		}
		else
		{
			destroy_particle(sys, pick);
			move = DESTROY_PARTICLE;
		}
	}
	return move;
}


//insert particle at random location
void create_particle( GCMC_System *sys)
{
        //we make a struct of type "particle"
	particle to_be_inserted; 
        //we create random coordinates for the particle 
	to_be_inserted.x[0] = randomish()*sys->box_side_length;
	to_be_inserted.x[1] = randomish()*sys->box_side_length;
	to_be_inserted.x[2] = randomish()*sys->box_side_length;
        //we add the particle to the vector that holds all our particles
 	sys->particles.push_back(to_be_inserted);
	return;
}


//Displace a random particle a random distance
void move_particle(GCMC_System *sys, int pick)
{
        double phi = (randomish()-0.5)* sys->box_side_length,
               gamma = (randomish()-0.5)*sys->box_side_length,
               delta = (randomish()-0.5)*sys->box_side_length;
        //store displacement in case it needs to be undone
        sys->move.pick = pick;
        sys->move.phi = phi;
        sys->move.gamma = gamma;
        sys->move.delta = delta;
        //make the moves 
        sys->particles[pick].x[0] += phi;
        sys->particles[pick].x[1] += gamma;
        sys->particles[pick].x[2] += delta;
        //check that the particle hasn't moved out of the box 
        for(int I=0;I<3;I++)
        {
            if(sys->particles[pick].x[I] > sys->box_side_length)
            {
                sys->particles[pick].x[I] -= sys->box_side_length;
            }
            else if(sys->particles[pick].x[I] < 0)
            {
                sys->particles[pick].x[I] += sys->box_side_length;
            }
        }
	return;
}

//particle deletion
void destroy_particle(GCMC_System *sys, int pick)
{
	sys->destroy.phi = sys->particles[pick].x[0];
	sys->destroy.gamma = sys->particles[pick].x[1];
	sys->destroy.delta = sys->particles[pick].x[2];
        //"begin" (below) points to the address of the zeroth item 
	// we move forward "pick" addresses 
	// to get to the address of the item we want
	sys->particles.erase(sys->particles.begin()+pick);
        return;
}


bool move_accepted(double cpe, double npe, MoveType move_type,\
                   GCMC_System *sys)
{
        double delta = npe - cpe,
               e = M_E,
               beta = 1.0 / (k * sys->system_temp),//thermodynamic beta
               volume = sys->volume,
               boltzmann_factor = pow(e,(-beta*delta)),
               acceptance;
	int pool = sys->particles.size(),
            poolplus = pool + 1;
	if (move_type == TRANSLATE )
	{
                acceptance = boltzmann_factor;
		if (acceptance > randomish())
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else if (move_type == CREATE_PARTICLE)
	{
                acceptance = boltzmann_factor* volume * conv_factor / \
                             (sys->system_temp * (double)pool);
		if (acceptance > randomish())
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else if (move_type == DESTROY_PARTICLE )//if we DESTROYED a particle
	{
                acceptance = boltzmann_factor * sys->system_temp * \
                             (double)poolplus / (volume * conv_factor);
		if (acceptance > randomish())
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	return true;
}

//undoes rejected moves based on move types
void undo_move(GCMC_System *sys, MoveType move)
{
	if (move == CREATE_PARTICLE)
	{
	        undo_insertion(sys);	
	}
	else if (move == TRANSLATE)
	{
		unmove_particle(sys);
	}
	else
	{
		particle added;
		added.x[0] = sys->destroy.phi;
		added.x[1] = sys->destroy.gamma;
		added.x[2] = sys->destroy.delta;
		sys->particles.push_back(added);
	}
	return;
}

void undo_insertion(GCMC_System *sys)
{
    sys->particles.pop_back();
}


//undoes a rejected displacement
void unmove_particle(GCMC_System *sys)
{
	int pick = sys->move.pick;
	double phi = sys->move.phi,
		gamma = sys->move.gamma,
		delta = sys->move.delta;
	sys->particles[pick].x[0] -= phi;
	sys->particles[pick].x[1] -= gamma;
	sys->particles[pick].x[2] -= delta;
        //check the particle is in the box
        for(int I=0;I<3;I++)
        {
            if(sys->particles[pick].x[I] > sys->box_side_length)
            {
                sys->particles[pick].x[I] -= sys->box_side_length;
            }
            else if(sys->particles[pick].x[I] < 0)
            {
                sys->particles[pick].x[I] += sys->box_side_length;
            }
        }
	return;
}

double sphere_volume(GCMC_System *sys,double diameter)
{
        diameter = diameter*sys->BinSize;//weighting
	double radius = diameter;//this is important!
        return (4.0 / 3.0) * M_PI * radius * radius * radius;
}


void radialDistribution(GCMC_System *sys,int step)
{
	const int nBins = sys->nBins; //total number of bins
	int IK,
            n = sys->particles.size(),
            num_pairs = 0;
	double  BinSize = sys->BinSize,
		expected_number_of_particles,
		diameter_of_current_sphere,
	        diameter_of_previous_sphere,
		shell_volume,
		dist;
	for (int I = 0; I<n - 1; I++)
	{
		for (int K = I + 1; K<n; K++)
		{
                    dist = distfinder(sys, I, K);
                    if(dist>sys->cutoff)
                    {
                        continue;
                    }
                    IK = int(dist / BinSize);
                    num_pairs+=1;
                    sys->boxes[IK] += 2;
		}
	}
	diameter_of_previous_sphere = 0;
        
        if(step==sys->maxStep-1)
        {
            double real_density_hours = (sys->sumparticles/(sys->maxStep*.5))/sys->volume;
            for (int I = 1; I <= nBins; I++)
            {
                    sys->boxes[I-1] /=  (sys->maxStep*.5);
                    sys->boxes[I-1] /= sys->sumparticles/(sys->maxStep*.5);
                    fprintf(sys->unweightedradial, "%lf\t%lf\n",BinSize*(I-1),\
                            sys->boxes[I - 1]);
                    diameter_of_current_sphere = I;
                    shell_volume = sphere_volume(sys,diameter_of_current_sphere)
                               - sphere_volume(sys,diameter_of_previous_sphere);
                    expected_number_of_particles = shell_volume * real_density_hours;
                    sys->boxes[I - 1] /= (expected_number_of_particles);
                    fprintf(sys->weightedradial, "%lf\t%lf\n",BinSize*(I-1),\
                            sys->boxes[I - 1]);
                    diameter_of_previous_sphere++;//update sphere diameter
            }
        }
	return;
}

void output(GCMC_System *sys, double accepted_energy, int step)
{
	int pool = sys->particles.size();
        if(pool==0)
        {
            return;
        }

	fprintf(sys->positions, "%d \n\n", pool);
	for (int p = 0; p<pool; p++)
	{
		fprintf(sys->positions, "%s %lf %lf %lf\n",\
                        sys->particle_type, sys->particles[p].x[0],\
                        sys->particles[p].x[1], sys->particles[p].x[2]);
	}

        double average_num_particles = sys->sumparticles/(step);
        fprintf(sys->particlecount, "%d %lf\n",step,average_num_particles);
        
        double average_energy = sys->sumenergy/step;
        fprintf(sys->energies, "%d %lf\
                %lf\n",step,accepted_energy,average_energy);
	return;
}
