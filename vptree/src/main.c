#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "../inc/vptree.h"



int main()
{
double  time;
struct timeval startwtime, endwtime;


int n,d;
printf("Press the number of points\n");
scanf("%d",&n);
printf("Press the number of dimensions\n");
scanf("%d",&d);




double *array=(double *)malloc(n*d*sizeof(double));
for(int i=0;i<n;i++){
    for(int j=0;j<d;j++){
        array[i*d+j]=(double)rand()/(double)RAND_MAX;
   //   printf("%lf ",array[i*d+j]);
    }
    printf("\n");
}






 gettimeofday (&startwtime, NULL);
 vptree *tree ;



tree=buildvp(array,n,d);



 gettimeofday (&endwtime, NULL);
    time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
  		      + endwtime.tv_sec - startwtime.tv_sec);

    printf("\t\tExecution Time: %f sec\n",time);






//printf("TELEIOSES!\n");
deleteTree(tree);
free(array);
    return 0;
}


