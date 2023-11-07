
#include "omp.h"
#include "stdio.h"

void printMP()
{
    omp_set_num_threads(4);
#pragma omp parallel
    printf("Hello from openMP\n");
}