#include <vector>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

std::vector<struct> particles;

struct particle
{
    double x[3]
};

int main()
{
    int i;
    for(i=0;i<10;i++)
    {
        double q,r;
        q = rand();
        r = q*i;
        v.push_back(r);
        printf("This is the %d element of v: %lf\n",i,v.at(i));
    }
    return 0;
}

