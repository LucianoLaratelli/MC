/*****
This is a monte carlo method simulator of the grand canonical ensemble, also known as the "mu-V-T" ensemble. This program takes a system composed of any number of particles and does one of three moves on it:

1. moves a random particle a random distance
2. adds a particle at a random positions
3. removes a random particle

It outputs :
the number of particles (both per-frame and a moving average)
the energies (both per-frame and as a moving average)
*****/

#include <vector>

std::vector<double> particle;



void random_p_mover()
{
    return;
}

void particle_maker()
{
    return;
}

void particle_remover()
{
    return;
}

//move chooser picks a random float between 0 and three and uses it to call one of three other functions: random_p_mover, particle_maker, or particle_remover. move chooser is called by main (tentative)
void move_chooser()
{
    double choice;
    //choice is a random float between 0 and 3
    if(choice<1)
    {
        //one of them
    }
    if(choice>2)
    {
        //another one
    }
    else
    {
        //the last one
    }
    return;
}
int main()
{
    return 0;
}
