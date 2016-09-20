#include <stdio.H>
#include <math.H>
#include <stdlib.H>

/*******************************************************************************
This is a simple program that places p atoms in a cube of volume m^3 with spacing
s. It writes these values to a text file startingpositions.txt that is read
by a monte carlo or other molecular simulation program in the same folder

values of N should be numbers whose cube root is an integer, sorry about that.
********************************************************************************/


int main()
{
    FILE * startingpositions;
    startingpositions = fopen("startingpositions.txt", "w");
    int N,s;//particles, spacing, max (length), volume
    double V,l;
    int i,j,k;//loop variables
    printf("How many particles do you want?\n");
    scanf("%d", &N);
    V = pow(N,.3333333333333);
    printf("How far should each particle be away from its neighbors?\n");
    scanf("%d", &s);
    for(i=0;i<V;i++)
    {
        for(j=0;j<V;j++)
        {
            for(k=0;k<V;k++)
            {
                int is = i*s;
                int js = j*s;
                int ks = k*s;
                fprintf(startingpositions,"%d %d %d\n", is, js, ks);
            }
        }
    }
    fclose(startingpositions);
    return 0;
}

