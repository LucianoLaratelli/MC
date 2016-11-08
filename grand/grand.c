/*****
This is a monte carlo method simulator of the grand canonical ensemble, also known as the "mu-V-T" ensemble. This program takes a system composed of any number of particles and does one of three moves on it:

1. moves a random particle a random distance
2. adds a particle at a random positions
3. removes a random particle

It outputs : the number of particles per frame, the average number of particles,
and the standard deviation associated with the average
the energies (in the same manner as the number of particles)
the chemical potential (per frame)
a calculated QST per-frame and as an average over all frames
*****/

#include <vector>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <random>

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

//we start with a box with sides of Length L
const int L = 2;
const int T = 1;
std::vector <particle> particles;
struct move_values move;
struct create_values creator;
struct remove_values destroy;

//positionchecker imposes periodic boundary conditions on a particle N by its id. 
//if any particle tries to escape the box, it disappears and is replaced on the
//opposite side of the box by a more obedient particle
bool positionchecker(int particleID)
{
    bool a;
    int i;
    for(i=0;i<3;i++)
    {
        if(particles[particleID].x[i] >= L)
        {
            particles[particleID].x[i] -= L;
            a = false;
        }
        if(particles[particleID].x[i] < 0)
        {
            particles[particleID].x[i] += L;
            a = false;
        }
        else
        {
            a = true;
        }
    }
    return a;
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
    int pick;
    double phi, gamma, delta;
    pick = move.pick;
    phi = move.phi;
    gamma = move.gamma;
    delta = move.gamma;
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
        double phi   = (randomish()-0.5)*L;
        double gamma = (randomish()-0.5)*L;
        double delta = (randomish()-0.5)*L;
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
        pool;
    double choice,
           pick;
    pool = particles.size();//amount of particles currently available
    if(pool == 0)
    {
        thecreator();
        pool = particles.size();
        flag = 1;
    }
    else
    {
        pick = rand()%pool;//picks random particle from those available
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
           pool = particles.size();
    int i;
    for(i=0;i<3;i++)
    {
        delta = particles[id_a].x[i] - particles[id_b].x[i];
        delta2 += delta * delta;
    }
    dist = sqrt(delta2);
    if(dist>0.5*L)
    {
        dist -= L;
    }
    return dist;
}

double pecalc()
{
    int b, //counter variables
        c;
    double pool = particles.size();
    double sigma = 1.00,//variables for the lennard-jones potential
           epsilon = 1.00,
           r,
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
    s2 = sigma * sigma;
    s6 = s2 * s2 * s2;
    s12 = s6 * s6;
    for(b = 0; b < pool-1; b++)
    {
        for(c = b + 1; c < pool; c++)
        {
            r = distfinder(b,c);
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

bool move_acceptor(double cpe, double npe, int c)
{
    double delta,
           guess,
           probability,
           beta,
           k;
    FILE * energies;
    energies = fopen("energies.dat","a");
    delta = npe - cpe;
    if(delta < 0)
    {
        printf("npe ifforce=%f\n",npe);
        fprintf(energies,"%d %f \n", c,npe);
        fclose(energies);
        return true;
    }
    guess = ((double)rand()/(double)RAND_MAX);
    k = 1;
    beta = 1 /(k * T);
    probability  = exp(-1 * beta * delta);
    printf("npe=%f\n",npe);
    if(probability > guess)
    {
        printf("npe if=%f\n",npe);
        fprintf(energies,"%d %f \n", c,npe);
        fclose(energies);
        return true;
    }
    else
    {
    printf("cpe else=%f\n",cpe);
        fprintf(energies,"%d %f \n", c,npe);
        fclose(energies);
        return false;
    }
}
/*
double average(double sum, int frames)
{
    return average;
}
*/
void output()
{
    FILE * positions;
    positions = fopen("positions.xyz","a");
        
    return;
}

int main()
{
    bool guess;
    double cpe,
           npe,
           sum,
           average;
    int c,
        m,
        max,
        flag;
    FILE * positions;
    FILE * energies;
    FILE * qsts;
    positions = fopen("positions.xyz","w");
    energies = fopen("energies.dat","w");
    qsts = fopen("qsts.dat","w");
    fclose(positions);
    fclose(qsts);
    cpe = pecalc();//we find the starting potential and call it our "current" one
    fprintf(energies,"0 %f\n", cpe);
    fclose(energies);
    printf("How many tries do you want?\n");
    clock_t begin = clock();
    fflush(stdout);//forces printf out of buffer
    scanf("%i",&max);
    m = particles.size();
    sum = cpe;//the sum has to include the initial starting energy
    for(c=1;c<max;c++)
    {
        while(m>=0)//we can't have negative particles
        {
            flag = move_chooser();//first step of the monte carlo, also stores which move we did in case we need to undo it
            npe = pecalc();//"new potential energy"
            guess = move_acceptor(cpe,npe,c);//move acceptor rejects or accepts the move
            if(guess == true)
            {
                output();
                cpe = npe;//updates the current energy
                sum += cpe;//adds to the average
                break;//breaks out of the while loop, increments c
            }
            else
            {
                move_undoer(flag);//uses the flag to undo a move if necessary
                output();
                sum += cpe;//adds the old energy to the average
                break;//breaks out of the while loop, increments c
            }
        }
    }
    clock_t end = clock();
    double time_spent = (double)(end-begin) / CLOCKS_PER_SEC;
    printf("Done! This run took %f seconds. Have a nice day!\n", time_spent);
    return 0;
}
