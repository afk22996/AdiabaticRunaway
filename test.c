#include <stdio.h>
#include <stdlib.h>

int main(void){
	double *x = (double *)malloc(2*sizeof(double));

	printf("%.15f, %.15f\n", x[0], x[1]);
	return 0;
}