

#include "../inc/vptree.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

#define THRESHOLD 100000
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

/*
 is a subtree. Features of this are coordinates, median distance, index of dataset, the left subtree-child and the right subtree
child 
*/




/*
 This funcion swaps rows of an 2D integer matrix.
*/



void swap_rows(int *index, int x, int y){
    int temp;

    temp=index[x];
    index[x]=index[y];
    index[y]=temp;
   
}


/*
This function is a component of quick select algorithm as shown in https://en.wikipedia.org/wiki/Quickselect .
It places the numbers which are smaller than the number on pivotIndex in the left-hand positions. The largest numbers are on the right. 
The numbers are not sorted.
*/

int partition_array(double *list, int left, int right, int pivotIndex,int *index,int d){
    double pivotValue=list[pivotIndex];
    double temp=list[pivotIndex];
    list[pivotIndex]=list[right];
    list[right]=temp;

    int itemp=index[pivotIndex];
    index[pivotIndex]=index[right];
    index[right]=itemp;


    int storeIndex=left;
    for(int i=left;i<right;i++){
        if(list[i]<pivotValue){
             temp=list[storeIndex];
            list[storeIndex]=list[i];
            list[i]=temp;

            itemp=index[storeIndex];
            index[storeIndex]=index[i];
            index[i]=itemp;


            storeIndex++;
        }
    }
     temp=list[right];
    list[right]=list[storeIndex];
    list[storeIndex]=temp;

    itemp=index[right];
    index[right]=index[storeIndex];
    index[storeIndex]=itemp;


    return storeIndex;
}

/*
This function is a component of quick select algorithm as shown in https://en.wikipedia.org/wiki/Quickselect.
Returns the k-th smallest element of list within left..right inclusive.
*/
double selectFunction(double *list,int left, int right,int k,int *index,int d){
    if(left==right){
        return list[left];
    }
    int pivotIndex=left + floor(rand() % (right - left + 1));
    pivotIndex=partition_array(list,left,right,pivotIndex,index,d);
    if(k==pivotIndex){
        return list[k];
    }else if(k<pivotIndex){
        return selectFunction(list,left,pivotIndex-1,k,index,d);
    }else{
        return selectFunction(list,pivotIndex+1,right,k,index,d);
    }
}


/*
It is a retrospective function used to create the tree from top to leaf.
 From an array of size n, select the last element of the tree as root. 
 The median distance is equal to the median of the distance vector if the number of the table is odd. 
 Otherwise the  n / 2 -1  element is selected.
 This is happening because the internal points remains the same while the computational cost is reduced compared to the median.
*/
 vptree *recursevp(double *X, int *index ,int n, int d){


   if(n==0){
         vptree *tree = malloc(sizeof *tree);
        tree=NULL;
        return tree;
    }
    if(n==1){
        vptree *tree = malloc(sizeof *tree);
        tree->vp=(double *)malloc(d*sizeof(double));
 
        for(int i=0;i<d;i++){
            tree->vp[i]=X[index[n-1]*d+i];
  
        }
        tree->md=0;

        tree->inner=NULL;
        tree->outer=NULL;
        tree->idx=index[0];
 
        return tree;

    }


       vptree *tree = malloc(sizeof *tree);
 
   tree->vp=(double *)malloc(d*sizeof(double));

   for(int i=0;i<d;i++){
    tree->vp[i]=X[index[n-1]*d+i];
    tree->idx=index[n-1];


   }


   double *distances=(double *)malloc((n-1)*sizeof(double));
	if(n*d>THRESHOLD){
   #pragma omp shared(distances,X,index) private(i,sum)
	{
		#pragma omp parallel for schedule(static)
	
   for (int i=0;i<n-1;i++){
        double sum=0;
        for(int j=0;j<d;j++){
            sum+=(X[index[i]*d+j]-X[index[n-1]*d+j])*(X[index[i]*d+j]-X[index[n-1]*d+j]);
        }

        distances[i]=sqrt(sum);

   }
   
	}
}else{
	for (int i=0;i<n-1;i++){
        double sum=0;
        for(int j=0;j<d;j++){
            sum+=(X[index[i]*d+j]-X[index[n-1]*d+j])*(X[index[i]*d+j]-X[index[n-1]*d+j]);
        }

        distances[i]=sqrt(sum);

   }
	
}
   
	
   double median=0;
  median=selectFunction(distances,0,n-2, floor((n-2)/2),index,d);
 /*
    for(int i=0;i<n-1;i++){
        printf("Distances[i]= %lf \n",distances[i]);
    }
    */
    tree->md=median;
	free(distances);
    
   // printf("%lf |\n",median);



if(n*d>THRESHOLD){
	
   #pragma omp parallel sections
	{
		#pragma omp section
		tree->inner=recursevp(X,index,(int)floor(n/2),d);
		#pragma omp section
   tree->outer=recursevp(X,&index[(int)floor(n/2)],(int)floor((n-1)/2),d);
	
	}
	
	}else{
	tree->inner=recursevp(X,index,(int)floor(n/2),d);
   	tree->outer=recursevp(X,&index[(int)floor(n/2)],(int)floor((n-1)/2),d);
	
}

    return tree;
    free(tree->vp);
    free(tree);

}

/*
This function generates a vantage point tree of a dataset with n points and d dimensions.
*/
 vptree *buildvp(double *X,int n, int d){
    int *index=(int *)malloc(n*sizeof(int));
    for(int i=0;i<n;i++){
        index[i]=i;
    }

   return recursevp(X,index,n,d);
    free(index);
};

//! Return vantage-point subtree with points inside radius
/*!

  \param node   A vantage-point tree
  \return The vantage-point subtree
*/
 vptree *getInner(vptree *tree){
    return tree->inner;
};

//! Return vantage-point subtree with points outside radius
/*!

  \param tree   A vantage-point tree
  \return The vantage-point subtree
*/
 vptree *getOuter(vptree *tree){
    return tree->outer;
};

//! Return median of distances to vantage point 
/*!

  \param tree   A vantage-point tree
  \return The median distance
*/
double getMD(vptree *tree){
    return tree->md;
}

//! Return the coordinates of the vantage point
/*!

  \param tree   A vantage-point tree
  \return The coordinates [d-dimensional vector]
*/
double *getVP( vptree *tree){
    return tree->vp;
}

//! Return the index of the vantage point
/*!

  \param tree   A vantage-point tree
  \return The index to the input vector of data points
*/
int getIDX( vptree *tree){
    return tree->idx;
}

//! Delete the vantage point tree
/*!

  \param tree   A vantage-point tree
  
*/
void deleteTree(vptree *tree)
{
    if (tree == NULL) return;


    deleteTree(tree->inner);
    deleteTree(tree->outer);


    free(tree);
}



