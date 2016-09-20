#include <stdio.H>
#include <math.H>
#include <stdlib.H>

/*******************************************************************************
This is a simple program that places p atoms in a cube of volume m^3 with spacing
s. It writes these values to a text file startingpositions.txt that is read
by a monte carlo or other molecular simulation program in the same folder
********************************************************************************/


int main()
{
    FILE * startingpositions;
    startingpositions = fopen("startingpositions.txt", "w");
    int N,s;//particles, spacing, max (length), volume
    double V,l;
    int i,p,j,k;//loop variables
    printf("How many particles do you want?\n");
    scanf("%d", &N);
    printf("What is the volume of the box?\n");
    scanf("%lf", &V);
    l = pow(V,1/3);
    printf("l = %lf",l);
    printf("How far should each particle be away from its neighbors?\n");
    scanf("%d", &s);
    for(p=0;p<N;)    
    {
        for(i=0;i<l;i++)
        {
            for(j=0;j<l;j++)
            {
                for(k=0;k<l;k++)
                {
                    int is = i*s;
                    int js = j*s;
                    int ks = k*s;
                    fprintf(startingpositions,"%d %d %d\n", is, js, ks);
                    p++;
                }
            }
        }
    }
    return 0;
}

