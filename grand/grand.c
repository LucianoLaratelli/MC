/*************************************************************************************************************
This is a monte carlo method simulator of the grand canonical ensemble, also known as the "mu-V-T" ensemble.
This program takes a system composed of any number of particles and does one of three moves on it:

1. moves a random particle a random distance
2. adds a particle at a random positions
3. removes a random particle

It outputs: the number of particles per frame, the average number of particles,
and the standard deviation associated with the average
the energies (in the same manner as the number of particles)
the chemical potential (per frame)
a calculated QST per-frame and as an average over all frames
**************************************************************************************************************/

#include <vector>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <random>
#include <stdlib.h>

typedef struct particle_particle
{
    double x[3];
} particle;

struct move_values
{
    int pick;
    double phi,gamma,delta;
};

struct create_values
{
    int pick;
    double phi, gamma, delta;
};

struct remove_values
{
    int pick;
    double phi, gamma, delta;
};
const int L = 22; //length of one side of the box
const double k = 1.0; //boltzmann factor
const double h = 6.626 * pow(10,-34);//planck constant

std::vector <particle> particles;
struct move_values move;
struct create_values creator;
struct remove_values destroy;

// positionchecker imposes periodic boundary conditions on a particle N by its id.  
// if any particle tries to escape the box, it disappears and is replaced on the 
// opposite side of the box by a more obedient particle
bool positionchecker(int particleID) 
{ 
    for(int i=0;i<3;i++) 
    { 
        if(particles[particleID].x[i] >= L) 
        { 
            particles[particleID].x[i] -= L;
            return false;
        }
        else if(particles[particleID].x[i] < 0)
        {
            particles[particleID].x[i] += L;
            return false;
        }
    }
    return true;
}

//randomish returns a random float between 0 and 1
double randomish()
{
    static int first_time_through = 1;
    if(first_time_through)
    {
        srandom(time(NULL));
        first_time_through = 0;
    }
    double r = random();
    return r/((double)RAND_MAX + 1);
}


void particleunmover()
{
    int pick = move.pick;
    double phi = move.phi, 
           gamma = move.gamma,
           delta = move.delta;
    particles[pick].x[0] -=phi;
    particles[pick].x[1] -= gamma;
    particles[pick].x[2] -= delta;
    return;
}


//particlemover moves a random particle a random distance
void particlemover(int pick)
{
    //phi,gamma, and delta are random floats between -L/2 and L/2
    //now we use a "bool counter" to make sure we don't move a particle
    //out of the box. Not sure how well this will work but oh well
    int i = false;
    while(i == false)
    {
        double phi   = (randomish()-0.5)*L,
               gamma = (randomish()-0.5)*L,
               delta = (randomish()-0.5)*L;
        //these next few lines store our moves in a struct in case they suck and we
        //need to undo the move
        move.pick = pick;
        move.phi = phi;
        move.gamma = gamma;
        move.delta = delta;
        //we make the moves (finally!)
        particles[pick].x[0] += phi;
        particles[pick].x[1] += gamma;
        particles[pick].x[2] += delta;
        //lastly, we check that the particle hasn't moved out of the box with 
        //positionchecker
        i = positionchecker(pick);
        if(i==false)
        {
            particleunmover();
        }
    }
    return;
}


//the creator adds particles to the system if it is chosen to do so by move_chooser
//below. The coordinates of the new particle are random based on the length 
//of the box
void thecreator()
{
    particle added; //making a struct to push_back to the vector
    added.x[0] = randomish()*L;
    added.x[1] = randomish()*L;
    added.x[2] = randomish()*L;
    creator.phi = added.x[0];
    creator.gamma = added.x[1];
    creator.delta = added.x[2];
    particles.push_back(added);
    int pool = particles.size()-1; //the index number of the new particle
    creator.pick = pool;
    return;
}

void thedestroyer(int pick)
{
    destroy.pick = pick;
    destroy.phi = particles[pick].x[0];
    destroy.gamma = particles[pick].x[1];
    destroy.delta = particles[pick].x[2];
    particles.erase(particles.begin() + pick);//begin refers to the pointer of the
                                              //zeroth item; 
                                              // we add "pick" amount of bytes 
                                              // to get to the next pointer
    return;
}


//move chooser picks a random float between 0 and three and uses it to call one of
//three other functions, which impose the monte carlo method onto the system.
int move_chooser()
{
    int flag,
        pool = particles.size();
    if(pool == 0)//if there are 0 particles we must add some to get started
    {
        thecreator();
        pool = particles.size();
        flag = 0;
    }
    else
    {
        double pick = random()%pool,//picks random particle from those available
               choice = randomish()* 3; //choice is a random float between 0 and 3
        fflush(stdout);
        if(choice<1.0)
        {
            particlemover(pick);
            flag = 0;
        }
        else if(choice>2.0)
        {
            thecreator();
            pool = particles.size();
            flag = 1;
        }
        else
        {
            thedestroyer(pick);
            flag = 2;
        }
    }
    return flag;
}

void move_undoer(int flag)
{
    if(flag == 0)
    {
        particleunmover();
    }
    else if(flag == 1)
    {
        thedestroyer(creator.pick);
    }
    else
    {
        particle added;
        added.x[0] = creator.phi;
        added.x[1] = creator.gamma;
        added.x[2] = creator.delta;
        particles.push_back(added);
    }
    return;
}

double distfinder(int id_a, int id_b)
{
    double dist,
           delta,
           delta2 = 0.00,
           cutoff = 0.5 * L;
    for(int i=0;i<3;i++)
    {
        delta = particles[id_a].x[i] - particles[id_b].x[i];
        if(fabs(delta) >= cutoff)
        {
            if(delta>0)
            {
                delta -= L;
            }
            if(delta<0)
            {
                delta += L;
            }
        }
        delta2 += delta * delta;
    }
    dist = sqrt(delta2);
    printf("%lf \n", dist);
    return dist;
}

double pecalc(double sigma, double epsilon)//returns energy in Kelvin
{
    int b, //counter variables
        c,
        pool = particles.size();
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
    double cutoff = L * 0.5;
    s2 = sigma * sigma;
    s6 = s2 * s2 * s2;
    s12 = s6 * s6;
    for(b = 0; b < pool-1; b++)
    {
        for(c = b + 1; c < pool; c++)
        {
            r = distfinder(b,c);
            if(r<cutoff)
            {
                continue;
            }
            rinv = 1 / r;
            r2 = rinv * rinv;
            r6 = r2 * r2 * r2;
            r12 = r6 * r6;
            sor6 = s6 * r6;
            sor12 = s12 * r12; 
            pe += (4 * epsilon * (sor12 - sor6));
        }
    }
    return pe;
}

bool move_acceptor(double cpe, double npe, int c, int flag, int n, double particle_mass,double system_temp)
{
   double  T = system_temp,
           delta,
           guess,
           probability,
           mu,
           beta = 1 / (k * T);
    FILE * energies;
    energies = fopen("energies.dat","a");
    delta = npe - cpe;
    double relativemu,
          m,
          lambda,
          lambdacubed;
    // NEED TO CALCULATE MU PER MOVE (IT SHOULD BE MINIMIZED AFTER N MOVES)
    int pool = n;
    int poolplus = pool + 1;
    double volume = L * L * L;
    double density = pool / volume;
    double pi = M_PI;
    m = particle_mass;//AMU
    lambda = (h)/(sqrt(2*pi*m*k*T)) ;
    lambdacubed = lambda * lambda * lambda;
    mu = k * T * log(lambdacubed) * density;
    relativemu = mu -( k * T * log(lambdacubed));
    if(delta < 0)
    {
        fprintf(energies,"%d %f \n", c,npe);
        fclose(energies);
        return true;
    }
    else if(flag == 0)//if we MOVED a particle 
    {
        guess = ((double)random()/(double)RAND_MAX);
        probability  = exp(-1 * beta * delta);
        if(probability > guess)
        {
            fprintf(energies,"%d %f \n", c,npe);
            fclose(energies);
            return true;
        }
        else
        {
            fprintf(energies,"%d %f \n", c,npe);
            fclose(energies);
            return false;
        }
    }
    else if(flag == 1)//if we CREATED a particle
    {
        double acceptance = -(beta*delta) + (beta*relativemu) +(log(volume/poolplus));
        double condition = exp(acceptance);
        if(condition < 1)
        {
            fprintf(energies,"%d %f \n", c,npe);
            fclose(energies);
            return true;
        }
        else
        {
            fprintf(energies,"%d %f \n", c,npe);
            fclose(energies);
            return false;
        }
    }
    else if(flag == 2)//if we DESTROYED a particle
    {
        double oneterm = -(beta*delta)-(beta*relativemu),
               acceptance = exp(oneterm),
               condition = pool/ (volume * acceptance);
        if(condition < 1)
        {
            fprintf(energies,"%d %f \n", c,npe);
            fclose(energies);
            return true;
        }
        else
        {
            fprintf(energies,"%d %f \n", c,npe);
            fclose(energies);
            return false;
        }
    }
    return true;
}

double qst_calc(int N, double energy, int c,double system_temp)//sorry
{
    FILE * qsts;
    qsts = fopen("qsts.dat","a");
    double T = system_temp,
           average_N = N/c,//<N>
           average_N_all_squared = average_N * average_N,//<N>^2
           average_energy = energy / c,//<U>
           N_squared = average_N * average_N,
           average_of_N_squared = N_squared * N_squared,//<N^2>
           particles_by_energy = average_N * average_energy,
           average_of_particles_by_energy = particles_by_energy / c, //<NU>
           numerator = (average_of_particles_by_energy) - (average_N * average_energy),
           denominator = average_of_N_squared - average_N_all_squared,
           QST = k * T - ( numerator / denominator);
    fprintf(qsts, "%d %lf\n", c, QST);
    fclose(qsts);
    return QST;
}

double sphere_volume(int radius)
{
    return 1.3 * M_PI * radius * radius * radius;
}

void histogram(int n)
{
    double num_density = (double)n / (double)(L*L*L);
    double expected_number_of_particles;
    double boxes[L] = {0};//initialize your fucking variables 
    for(int I = 0;I<n-1;I++)
    {
        for(int K = I + 1;K<n;K++)
        {
            double dist = distfinder(I,K);
            for(int M = 1;M<L;M++)
            {
                if(dist < (L*.5) && dist < M && dist>=(M-1))
                {
                    boxes[M-1]++;
                }
            }
        }
    }
    printf("number of particles: %d\n", n);
    double previous_sphere_volume = 0;
    for(int I = 1;I<=L;I++)
    {
        printf("%lf ",boxes[I-1]);
        double current_sphere_volume = sphere_volume(I);
        double shell_size = current_sphere_volume - previous_sphere_volume;
        previous_sphere_volume = current_sphere_volume;
        expected_number_of_particles = shell_size * num_density;
        boxes[I-1] /= expected_number_of_particles;
        printf("%f ", boxes[I-1]);
    }
    return;
}

void output(char *particle_type)
{
    FILE * positions;
    positions = fopen("positions.xyz","a");
    int p,
        pool;
    pool = particles.size();
    fprintf(positions, "%d \n\n",pool);
    for(p=0;p<pool;p++)
    {
        fprintf(positions, "%s %lf %lf %lf\n",particle_type, particles[p].x[0],particles[p].x[1], particles[p].x[2]);
    }
    fclose(positions); 
    return;
}

int main(int argc, char *argv[])
{
    if(argc != 6)
    {
        printf("This program takes five arguments: the type of particle, its mass, a LJ sigma, LJ epsilon, and the temperature of the system, in that order.");
        exit(EXIT_FAILURE);
    }
    bool guess;
    char particle_type[128] = {0};
    int c,
        n,
        max,
        flag;
    double particle_mass,
           sigma,
           epsilon,
           system_temp,
           cpe,
           npe,
           sumenergy,
           sumparticles;
    sscanf(argv[1],"%s",particle_type);
    sscanf(argv[2],"%lf",&particle_mass);
    sscanf(argv[3],"%lf",&sigma);
    sscanf(argv[4],"%lf",&epsilon);
    sscanf(argv[5],"%lf",&system_temp);
    FILE * positions;
    FILE * energies;
    FILE * qsts;
    positions = fopen("positions.xyz","w");
    energies = fopen("energies.dat","w");
    qsts = fopen("qsts.dat","w");
    fclose(positions);
    fclose(qsts);
    cpe = pecalc(sigma, epsilon);//we find the starting potential and call it our "current" one
    fprintf(energies,"0 %f\n", cpe);
    fclose(energies);
    printf("How many tries do you want?\n");
    clock_t begin = clock();
    fflush(stdout);//forces printf out of buffer
    scanf("%i",&max);
    n = particles.size();
    sumenergy = cpe;//the sum has to include the initial starting energy
    sumparticles = n;
    for(c=1;c<max;c++)
    {
        while(n>=0)//we can't have negative particles
        {
            flag = move_chooser();//first step of the monte carlo, also stores which move we did in case we need to undo it
            npe = pecalc(sigma,epsilon);//"new potential energy" in K
            guess = move_acceptor(cpe,npe,c,flag,n,particle_mass,system_temp);//move acceptor rejects or accepts the move
            if(guess == true)
            {
                n = particles.size();
                output(particle_type);
                cpe = npe;//updates the current energy
                sumenergy += cpe;//adds to the running total
                sumparticles += n;//see above
                qst_calc(sumparticles,sumenergy,c,system_temp);
                break;//breaks out of the while loop, increments c
            }
            else
            {
                n = particles.size();
                move_undoer(flag);//uses the flag to undo a move if necessary
                output(particle_type);
                sumenergy += cpe;//adds the old energy to the average
                sumparticles += n;
                qst_calc(sumparticles,sumenergy,c,system_temp);
                break;//breaks out of the while loop, increments c
            }
        }
    }
    histogram(n);
    clock_t end = clock();
    double time_spent = (double)(end-begin) / CLOCKS_PER_SEC;
    printf("Done! This run took %f seconds. Have a nice day!\n", time_spent);
    return 0;
}
