#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*******************************************************************************
This is a simple program that places p atoms in a cube of volume m^3 with
spacing s. It writes these values to a text file startingpositions.txt that is
read by a monte carlo or other molecular simulation program in the same folder

values of N should be numbers whose cube root is an integer, sorry about that.
********************************************************************************/

int main() {
  FILE *startingpositions;
  startingpositions = fopen("startingpositions.txt", "w");
  int N; // particles, spacing, max (length), volume
  double V, l, s;
  int i, j, k; // loop variables
  printf("How many particles do you want?\n");
  scanf("%d", &N);
  V = pow(N, .3333333333333);
  printf("How far should each particle be away from its neighbors?\n");
  scanf("%lf", &s);
  if (N == 2) {
    for (i = 0; i < N; i++) {
      int is = i * s;
      fprintf(startingpositions, "%d %d %d\n", is, is, is);
    }
  } else {
    for (i = 0; i < V; i++) {
      for (j = 0; j < V; j++) {
        for (k = 0; k < V; k++) {
          double is = i * s;
          double js = j * s;
          double ks = k * s;
          fprintf(startingpositions, "%f %f %f\n", is, js, ks);
        }
      }
    }
  }
  fclose(startingpositions);
  return 0;
}
