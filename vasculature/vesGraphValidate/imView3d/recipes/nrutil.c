/*************************************************************************
NAME: 		nrutil.c

DATE WRITTEN:	October 11, 1989

WRITTEN BY: 	Dimitrios M. Emiris (from NUMERICAL RECIPES IN C).

PURPOSE:	Here, a number of subroutines are included in order to customize
		certain trivial operations. nrerror(), provides the error 
		messages, ivector(),vector(), and dvector(), allocate memory 
		space for integer, float and double vectors, respectively, as 
		imatrix(), matrix(), and dmatrix(), do for matrices. submatrix()
		creates a submatrix of an original matrix, while the free_()
		subroutines free the allocated memory space..

ARGUMENTS:	indices for the matrices and vectors; vary depending upon the
		specific subroutine.

FUNCTIONS CALLED:  none

REMARKS: 	none 

***************************************************************************/

#include <math.h>
#include <malloc.h>
#include <stdio.h>
#include <math.h>

void nrerror(error_text)
char error_text[];
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}



float *vector(nl,nh)
int nl,nh;
{
	float *v;

	v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}

int *ivector(nl,nh)
int nl,nh;
{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}

double *dvector(nl,nh)
int nl,nh;
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}



float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	int i;
	float **m;

	m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	int i,**m;

	m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}



float **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
float **a;
int oldrl,oldrh,oldcl,oldch,newrl,newcl;
{
	int i,j;
	float **m;

	m=(float **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(float*));
	if (!m) nrerror("allocation failure in submatrix()");
	m -= newrl;

	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

	return m;
}



void free_vector(v,nl,nh)
float *v;
int nl,nh;
{
	free((char*) (v+nl));
}

void free_ivector(v,nl,nh)
int *v,nl,nh;
{
	free((char*) (v+nl));
}

void free_dvector(v,nl,nh)
double *v;
int nl,nh;
{
	free((char*) (v+nl));
}



void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
int nrl,nrh,ncl,nch;
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl,nrh,ncl,nch;
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}



void free_submatrix(b,nrl,nrh,ncl,nch)
float **b;
int nrl,nrh,ncl,nch;
{
	free((char*) (b+nrl));
}



float **convert_matrix(a,nrl,nrh,ncl,nch)
float *a;
int nrl,nrh,ncl,nch;
{
	int i,j,nrow,ncol;
	float **m;

	nrow=nrh-nrl+1;
	ncol=nch-ncl+1;
	m = (float **) malloc((unsigned) (nrow)*sizeof(float*));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m -= nrl;
	for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
	return m;
}



void free_convert_matrix(b,nrl,nrh,ncl,nch)
float **b;
int nrl,nrh,ncl,nch;
{
	free((char*) (b+nrl));
}

double ***dstruct(n1l,n1h,n2l,n2h,n3l,n3h)
int n1l,n1h,n2l,n2h,n3l,n3h;
{
	int i,j;
	double ***m;

	m=(double ***) malloc((unsigned) (n1h-n1l+1)*sizeof(double**));
	if (!m) nrerror("allocation failure 1 in dstruct()");
	m -= n1l;

	for (i=n1l;i<=n1h;i++){
		m[i]=(double **) malloc((unsigned) (n2h-n2l+1)*sizeof(double*));
		if (!m[i]) nrerror("allocation failure 2 in dstruct()");
		m[i] -= n2l;
		}

	for (i=n1l;i<=n1h;i++){
	   for (j=n2l;j<=n2h;j++){
	   	   m[i][j]=(double *) malloc((unsigned) (n3h-n3l+1)*sizeof(double));
		   if (!m[i][j]) nrerror("allocation failure 3 in dstruct()");
		   m[i][j] -= n3l;
		   }
	   }

	return m;
}

void free_dstruct(m,n1l,n1h,n2l,n2h,n3l,n3h)
double ***m;
int n1l,n1h,n2l,n2h,n3l,n3h;
{
	int i,j;

	for (i=n1h;i>=n1l;i--)
		for (j=n2h;j>=n2l;j--)
			free ((char*) (m[i][j]+n3l));

	for (i=n1h;i>=n1l;i--)
		free ((char*) (m[i]+n2l));

	free ((char*) (m+n1l));
}


