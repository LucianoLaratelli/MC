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
        o2[] = "O2",
        water[] = "Water";
   char * particle_type= sys->particle_type;
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
   else if(strcmp(particle_type,o2) == 0)//thanks Doug!
    {
        sys->sigma = 0.8;//angstroms
        sys->epsilon = 4.0;//kcal / mol, CONVERT TO PROPER UNITS
        sys->particle_mass = 15.9994;
    }
   //water is for Stockmeyer fluids only
   else if(strcmp(particle_type,water) == 0)
   {
       char oxygen[] = "O";//water is modeled by a single-point oxygen
       strcpy(sys->particle_type,oxygen);
       sys->sigma = 3.15100;
       sys->epsilon = 76.42000;
       sys->particle_mass = 18.016;
       sys->dipole_magnitude = 1.86*85.10597636;
       sys->polarizability = 0;
       sys->stockmayer_flag = true;//set the flag so we can never forget it
   }
   else
   {
       printf("Not a supported chemical species!\nAllowed values for Lennard-"\
               "Jones are:\nAr\nNe\nHe\nKr\nXe\nO2\nAllowed values for "\
               "Stockmeyer are:\nWater\nPlease try again!\n");
       exit(EXIT_FAILURE);
   }
   return;
}     

double calculate_PE(GCMC_System *sys)
{
        if(sys->ideal_flag)
        {
            return 0;
        }
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
        if(sys->stockmayer_flag)
        {
            double ** matrix = matrix_madness(sys);
            double correction = 0;
            for(int i = 0;i<pool;i++)
            {
                for(int j=i+1;j<pool;j++)
                {
                    if(i==j)
                    {
                        continue;
                    }
                    for(int p = 0;p<3;p++)
                    {
                        for(int q = 0;q<3;q++)
                        {
                            correction += sys->particles[i].dipole[p]*
                                          matrix[(3*i)+p][(3*j)+q]*
                                          sys->particles[j].dipole[q];
                        }
                    }
                }
            }
            printf("dipole dipole = %lf\njust lj = %lf\n",correction,pe);
            pe += correction; 
            free(matrix);
        }
	return pe;//in KELVIN
}

double distfinder(GCMC_System *sys, int id_a, int id_b)
{
	double dist = 0.0,
               delta = 0.0,
               delta2 = 0.00;

	for(int i = 0; i<3; i++)
	{
		delta = sys->particles[id_a].x[i] - sys->particles[id_b].x[i];
                //the following conditionals are the minimum image convention
		if (delta >= sys->cutoff)
		{
			delta -= sys->box_side_length;
		}
		else if (delta <= (-1*sys->cutoff))
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

//return a random double between min and max (thanks StackOverflow!)
double random_range(double min, double max)
{
	return min + (random() / ((double)RAND_MAX) * (max - min));
}

MoveType make_move(GCMC_System *sys)
{
        //MoveType is an enum in MonteCarlo.h 
	MoveType move;
	int pool = sys->particles.size();
        if(sys->NVT_flag)
        {
            double pick = random() % pool;
            move_particle(sys,pick);
            move = TRANSLATE;
        }
        else
        {
            if (pool == 0)
            {
                    create_particle(sys);
                    pool = sys->particles.size();
                    move = CREATE_PARTICLE;
            }
            else
            {
                double pick = random() % pool,//picks random particle
                       choice = random_range(0,1);//random float between 0 and 1
                fflush(stdout);
                if (choice<0.33333)
                {
                        create_particle(sys);
                        move = CREATE_PARTICLE;
                }
                else if (choice >= (.666667))
                {
                        move_particle(sys,pick);
                        move = TRANSLATE;
                }
                else
                {
                        destroy_particle(sys, pick);
                        move = DESTROY_PARTICLE;
                }
            }
        }
	return move;
}


double ** matrix_madness(GCMC_System *sys)
{
    int n = sys->particles.size();
    double ** matrix =(double**) malloc(sizeof(double*)*3*n);
    for(int i = 0;i<3*n;i++)
    {
        matrix[i] = (double*)calloc(3*n, sizeof(double));
    }
    for(int a = 0;a < 3 * n;a++)
    {
        if(sys->polarizability == 0)
        {
            matrix[a][a]= 1000000000000000000000000000000.f;
        }
        else
        {
            matrix[a][a] = 1/sys->polarizability;
        }
    }
    for(int i = 0;i<n-1;i++)
    {
        for(int j = i+1;j<n;j++)
        {
            double r = distfinder(sys,i,j),
                   rinv = 1/r,
                   rinv2 = rinv * rinv,
                   rinv3 = rinv2 * rinv,
                   rinv5 = rinv2 * rinv3,
                   del_x = sys->particles[i].x[0] - sys->particles[j].x[0],
                   del_y = sys->particles[i].x[1] - sys->particles[j].x[1],
                   del_z = sys->particles[i].x[2] - sys->particles[j].x[2];
                   if(del_x > sys->cutoff)
                   {
                       del_x -= sys->cutoff;
                   }
                   if(del_y > sys->cutoff)
                   {
                       del_y -= sys->cutoff;
                   }
                   if(del_z > sys->cutoff)
                   {
                       del_z -= sys->cutoff;
                   }
            double deltas[3] = {del_x, del_y, del_z};
            for(int p = 0; p<3;p++)
            {
                for(int q = 0;q<3;q++)
                {
                    matrix[(3*i)+p][(3*j)+q] = -3*rinv5*deltas[p]*deltas[q];
                    if(p==q)
                    {
                        matrix[(3*i)+p][(3*j)+q] += rinv3;
                    }
                    matrix[(3*j)+p][(3*i)+q] = matrix[(3*i)+p][(3*j)+q];
                }
            }
        }
    }
    return matrix;
    free(matrix);
}



//this function picks a random point on the surface of a unit sphere
double * pick_dipole_direction(GCMC_System *sys)
{
    double theta = random_range(0,2*M_PI),
           z = random_range(-1,1),
           x = sqrt(1-(z*z)) * cos(theta),
           y = sqrt(1-(z*z)) * cos(theta);
    double * dipole;
    dipole = (double*)malloc(sizeof(double) * 3);
    dipole[0] = x * sys->dipole_magnitude;
    dipole[1] = y * sys->dipole_magnitude;
    dipole[2] = z * sys->dipole_magnitude;
    return dipole;
}


//insert particle at random location
void create_particle( GCMC_System *sys)
{
    //we make a struct of type "particle"
    particle to_be_inserted; 
    //we create random coordinates for the particle 
    to_be_inserted.x[0] = random_range(0,sys->box_side_length);
    to_be_inserted.x[1] = random_range(0,sys->box_side_length);
    to_be_inserted.x[2] = random_range(0,sys->box_side_length);

    if(sys->stockmayer_flag)
    {
        double * dipole = pick_dipole_direction(sys);
        to_be_inserted.dipole[0] = dipole[0];
        to_be_inserted.dipole[1] = dipole[1];
        to_be_inserted.dipole[2] = dipole[2];
        free(dipole);
    }
    //we add the particle to the vector that holds all our particles
    sys->particles.push_back(to_be_inserted);
    return;
}


//Displace a random particle a random distance
void move_particle(GCMC_System *sys, int pick)
{
        double negative_half_box = -0.5 * sys->box_side_length;
        double phi = random_range(negative_half_box,sys->cutoff),
               gamma = random_range(negative_half_box,sys->cutoff), 
               delta = random_range(negative_half_box,sys->cutoff);
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
        if(sys->stockmayer_flag)
        {
            sys->move.dipole[0] = sys->particles[pick].dipole[0];
            sys->move.dipole[1] = sys->particles[pick].dipole[1];
            sys->move.dipole[2] = sys->particles[pick].dipole[2];

            double * dipole = pick_dipole_direction(sys);
            sys->particles[pick].dipole[0] = dipole[0];
            sys->particles[pick].dipole[1] = dipole[1];
            sys->particles[pick].dipole[2] = dipole[2];
            free(dipole);
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
               acceptance,
               random = random_range(0,1);
	int pool = sys->particles.size(),
            poolplus = pool + 1;
	if (move_type == TRANSLATE )
	{
                acceptance = boltzmann_factor;
		if (acceptance > random)
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
                acceptance =  boltzmann_factor* volume * conv_factor / \
                             (sys->system_temp * (double)pool);
		if (acceptance > random)
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
                acceptance =  boltzmann_factor * sys->system_temp * \
                             (double)poolplus / (volume * conv_factor);
		if (acceptance > random)
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
        if(sys->stockmayer_flag)
        {
            sys->particles[pick].dipole[0] = sys->move.dipole[0];
            sys->particles[pick].dipole[1] = sys->move.dipole[1];
            sys->particles[pick].dipole[2] = sys->move.dipole[2];
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

void output(GCMC_System *sys, double accepted_energy)
{
        if(sys->energy_output_flag)
        {
            fprintf(sys->energies, "%lf \n",accepted_energy);
        }
	int pool = sys->particles.size();
        if(pool==0)
        {
            return;
        }

        if(sys->output_flag)
        {
            fprintf(sys->output,"ITEM: TIMESTEP\n"
                                "%d\n",sys->step);
            fprintf(sys->output,"ITEM: NUMBER OF ATOMS\n"
                                "%d\n",(int)sys->particles.size());
            fprintf(sys->output,"ITEM: BOX BOUNDS pp pp pp\n"
                               "0 %lf\n0 %lf\n0 %lf\n",\
                               sys->box_side_length,sys->box_side_length,\
                               sys->box_side_length);
            fprintf(sys->output,"ITEM: ATOMS id type x y z mux muy muz\n");
            for (int p = 0; p<pool; p++)
            {
                    fprintf(sys->output, "%d 6 %lf %lf %lf %lf %lf %lf\n",\
                            p,sys->particles[p].x[0], sys->particles[p].x[1],\
                            sys->particles[p].x[2],\
                            sys->particles[p].dipole[0]/85.10597636,\
                            sys->particles[p].dipole[1]/85.10597636,\
                            sys->particles[p].dipole[2]/85.10597636);
            }
        }
	return;
}
