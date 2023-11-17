#include "stdlib.h"
#include "stdio.h"
#include "cilk/cilk.h"

void printCilk()
{
    cilk_for(int i = 0; i < 10; i++)
    {
        printf("%d\n", i);
    }
}