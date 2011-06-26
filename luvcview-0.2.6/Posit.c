/*****************************************************************************/

/* 

Code for computing object pose using POSIT 

Inputs: 1) a file containing number of points and coordinates of object points;
        2) a file containing integer coordinates of image points;
Important: Image point at position i in image file should correspond to object point at 
position i in object file

Outputs: 1) Rotation of scene with respect to camera
         2) Translation vector from projection center of camera
	    to FIRST point of scene (object) file.
            
Reference: D. DeMenthon and L.S. Davis, "Model-Based Object Pose in 25 Lines of Code", 
International Journal of Computer Vision, 15, pp. 123-141, June 1995. 

*/

/*****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "muhaha.h"
/*****************************************************************************/

const int maxCount = 30;/* exit iteration with a message after this many loops */
const int nbObjectCoords = 3; /* x, y, z */
const int nbImageCoords = 2; /* x, y */

/*****************************************************************************/

/* Structures */

/*****************************************************************************/

/* Prototypes: */

double **InitDoubleArray(int nbRows, int nbCols);
int **InitIntArray(int nbRows, int nbCols);
void FreeDoubleArray(double **anArray, int nbRows);
void FreeIntArray(int **anArray, int nbRows);
void PrintMatrix(double **mat, int nbRows, int nbCols);
void PrintVector(double vect[], int nbElems);
double **ReadObject(char *name, int *nbPts);
int **ReadImage(char *name, int imageCenter[2], int nbPts);
void GetObjectAndImage(char *objectName, char *imageName, TObject *object, TImage *image);
void POS(TObject object, TImage image, TCamera *camera);
void POSIT(TObject object, TImage image, TCamera *camera);
long GetImageDifference(TImage image);
int  PseudoInverse(double **A,int N, double **B);
void PrintRotation(double rotation[][3]);
void nrerror(char error_text[]);
int main();

/*****************************************************************************/


int debug1 = 1;/* 0 for no printout */

/*****************************************************************************/

double **ReadObject(char *name, int *nbPts)
/* 
Read object feature points from a file 
*/
{
FILE *fp;
int i, j;
float coord;
double **objectPts;
		
	fp = fopen(name, "r");
	if (fp != NULL){
                if(debug1) printf("Object data read from file:\n");
		int fu = fscanf(fp, "%d", nbPts);
		if(debug1) printf("%6d\n", *nbPts);
		
		/* Initialize objectPts array */
		objectPts = InitDoubleArray(*nbPts, nbObjectCoords);
		
		for (i=0; i< *nbPts;i++){
 			for (j=0; j<nbObjectCoords; j++){
				int fu = fscanf(fp, "%f", &coord);
				objectPts[i][j] = coord;
				if(debug1) printf("%.2f ",coord);
			}
			if(debug1) printf("\n");
		}
		if(debug1) printf("\n");
	}
	else nrerror("could not open file of object points");
	return objectPts;
}

/*****************************************************************************/

int **ReadImage(char *name, int imageCenter[2], int nbPts)
/*
read coords of image points, which are INTEGER
there should be as many image points as there were object points
and the order should be the same order as that of the corresponding 
object points of the file read by ReadObject.
*/
{
	FILE *fp;
	int i, j;
	int coord;
	int **imagePts;
	
	fp = fopen(name, "r");
	if (fp != NULL){
		imagePts = InitIntArray(nbPts, nbImageCoords);
                if(debug1) printf("Image data read from file:\n");
		for (i=0; i<nbPts;i++){
			for (j=0; j<2; j++){
				int fu = fscanf(fp, "%d", &coord);
				imagePts[i][j] = coord - imageCenter[j];
				if(debug1) printf("%d ", coord);
			}
			if(debug1) printf("\n");
		}
		if(debug1) printf("\n");
	}
	else printf("could not open file ");
	return imagePts;
}

/*****************************************************************************/

void GetObjectAndImage(char *objectName, char *imageName, TObject *object, TImage *image)
/* 
Initialize data structures for object and image data and read data from files 
*/
{
int i, j, nbPoints;
/* Object data */

	object->objectPts = ReadObject(objectName, &nbPoints);
	object->nbPts = nbPoints;

	object->objectVects = InitDoubleArray(object->nbPts, nbObjectCoords);
	object->objectCopy = InitDoubleArray(object->nbPts, nbObjectCoords);
		
	for (i=0;i<object->nbPts;i++){
		for(j=0;j<nbObjectCoords;j++){
			object->objectVects[i][j] = object->objectCopy[i][j] =
					object->objectPts[i][j] - object->objectPts[0][j];
		}
	}

	object->objectMatrix = InitDoubleArray(nbObjectCoords, object->nbPts);

/* Image data */
	
	image->nbPts = object->nbPts;
	image->imagePts = ReadImage(imageName, image->imageCenter, image->nbPts);

	image->imageVects = InitDoubleArray(image->nbPts, nbImageCoords);
	image->oldImageVects = InitDoubleArray(image->nbPts, nbImageCoords);
	image->epsilon = (double *)malloc(image->nbPts * sizeof(double));
}

	
/*****************************************************************************/

void POS(TObject object, TImage image, TCamera *camera)
/*
POS function 
(Pose from Orthography and Scaling, a scaled orthographic proj. approximation).
Returns one translation and one rotation.
*/
{
double    I0[3], J0[3], row1[3], row2[3], row3[3];
double    I0I0, J0J0;
int  i, j;
double scale, scale1, scale2;


	/*Computing I0 and J0, the vectors I and J in TRs listed above */
	for (i=0;i<nbObjectCoords;i++){
    	I0[i]=0;
    	J0[i]=0;
    	for (j=0;j<object.nbPts;j++){
			I0[i]+=object.objectMatrix[i][j]*image.imageVects[j][0];
			J0[i]+=object.objectMatrix[i][j]*image.imageVects[j][1];
    	}
	}

	I0I0=I0[0]*I0[0] + I0[1]*I0[1] + I0[2]*I0[2];
	J0J0=J0[0]*J0[0] + J0[1]*J0[1] + J0[2]*J0[2];

	scale1 = sqrt(I0I0);
	scale2 = sqrt(J0J0);
	scale = (scale1 + scale2) / 2.0;
	
	/*Computing TRANSLATION */
	camera->translation[0] = image.imagePts[0][0]/scale;
	camera->translation[1] = image.imagePts[0][1]/scale;
	camera->translation[2] = camera->focalLength/scale;

	/* Computing ROTATION */
	for (i=0;i<3;i++){
 	   row1[i]=I0[i]/scale1;
 	   row2[i]=J0[i]/scale2;
	}
	row3[0]=row1[1]*row2[2]-row1[2]*row2[1];/* Cross-product to obtain third row */
	row3[1]=row1[2]*row2[0]-row1[0]*row2[2];
	row3[2]=row1[0]*row2[1]-row1[1]*row2[0];
	for (i=0;i<3;i++){ 
 	   camera->rotation[0][i]=row1[i];
 	   camera->rotation[1][i]=row2[i];
 	   camera->rotation[2][i]=row3[i];
	}
}

/*****************************************************************************/

long GetImageDifference(TImage image)
/*
Get sum of differences between coordinates in lists
of old image points and new image points
*/
{
int i, j;
long sumOfDiffs = 0;
	
	for (i=0;i<image.nbPts;i++){
		for (j=0;j<2;j++){
			sumOfDiffs += abs(floor(0.5+image.imageVects[i][j])-floor(0.5+image.oldImageVects[i][j]));
		}
	}
	return sumOfDiffs;
}

/*****************************************************************************/

void POSIT(TObject object, TImage image, TCamera *camera)
/* 
Iterate over results obtained by the POS function;
see paper "Model-Based Object Pose in 25 Lines of Code", IJCV 15, pp. 123-141, 1995.
*/
{
int i, j, iCount;
int converged;
long imageDiff;

/* Starting point for iteration loop */
	for(iCount=0;iCount<maxCount;iCount++){
		if(iCount==0){
			for (i=0;i<image.nbPts;i++){
				for(j=0;j<nbImageCoords;j++){
   				image.imageVects[i][j]=(double)(image.imagePts[i][j]-image.imagePts[0][j]);
				}
   	 	}
   	 }
   	 else{/* iCount>0 */
			/* Compute new imageVects */
			for (i=0;i<image.nbPts;i++){
				image.epsilon[i] = 0.0;
				for (j=0;j<3;j++){
 	  				image.epsilon[i] += object.objectVects[i][j] * camera->rotation[2][j]; /*dot product M0Mi.k*/
 	  			}
 	  			image.epsilon[i] /= camera->translation[2]; /* divide by Z0 */
 	  		}
 			/* Corrected image vectors */	
			for (i=0;i<image.nbPts;i++){
				for(j=0;j<nbImageCoords;j++){
	    			image.imageVects[i][j]=(double)image.imagePts[i][j]*(1+image.epsilon[i])-image.imagePts[0][j];
				}
 			}
			imageDiff = GetImageDifference(image);/*using pts gives same result */
			printf("imageDiff %ld\n",imageDiff);
		}
			
		/* Remember old imageVects */
		for (i=0;i<image.nbPts;i++){
    		image.oldImageVects[i][0]=image.imageVects[i][0];
   	 	image.oldImageVects[i][1]=image.imageVects[i][1];
		}
		
		/* Compute rotation and translation */
		POS(object,image, camera);
	
 		converged = (iCount>0 && imageDiff==0);
 		if(converged) break;

		if(iCount==maxCount && !converged){
			printf("POSIT did not converge \n");
		}
	}/*end for*/;
}

/*****************************************************************************/

void PrintRotation(double rotation[][3])
{
int i, j;

	for (i=0; i<3;i++){
		for (j=0; j<3; j++){
			printf("%.3f ", rotation[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

/*****************************************************************************/



#if 0
int main()
/*
Read file of image points, file of object points,
compute object matrix by pseudoinverse operation
then get the pose
*/
{
int fLength = 760;/* Focal Length in PIXELS */
int imageCenterX = 0;/* In this test, origin of image coordinates is at image center */
int imageCenterY = 0;/* More often, imageCenterX=256, imageCenterY=256 */
int isMaxRank;
TObject testObject;
TImage testImage;
TCamera testCamera;

	testImage.imageCenter[0] = imageCenterX;
	testImage.imageCenter[1] = imageCenterY;

	/* Init structures and read object and image data:*/	
	GetObjectAndImage("TestCube", "CubeImage", &testObject, &testImage);
		
	testCamera.focalLength = fLength;

	/* Find object matrix as pseudoInverse of matrix of object vector coordinates */	
	isMaxRank = PseudoInverse(testObject.objectCopy, testObject.nbPts, testObject.objectMatrix);
	/* Check if matrix rank for object matrix is less than 3: */
	if(!isMaxRank){
		nrerror("object is too flat; another method is required ");/* exit */
	}
	else{
		/* Find object pose by POSIT, Pose from Orthography, Scaling, and ITerations */
		POSIT(testObject, testImage, &testCamera);

		/* Print results */
		printf("\n");
		printf("Rotation matrix of scene in camera reference frame\n");
		PrintRotation(testCamera.rotation);
		printf("\n");
		printf("Translation vector of scene in camera reference frame\n" 
                        "from camera projecction center to FIRST point of scene (object) file\n");
		PrintVector(testCamera.translation, nbObjectCoords);
		printf("\n");
	}
	/* Do a bit of memory cleanup: */
	FreeDoubleArray(testObject.objectPts, testObject.nbPts);
	FreeDoubleArray(testObject.objectVects, testObject.nbPts);
	FreeDoubleArray(testObject.objectCopy, testObject.nbPts);
	FreeDoubleArray(testObject.objectMatrix, nbObjectCoords);
	FreeIntArray(testImage.imagePts, testImage.nbPts);
	FreeDoubleArray(testImage.imageVects, testImage.nbPts);
	FreeDoubleArray(testImage.oldImageVects, testImage.nbPts);
	
        return 0;	
}
#endif

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/****************************************************************************/

/* 
Refer to "Numerical Recipes"
W.H. Press, B.P. Flannery, S.A. Teukolsky, et W.T. Vetterling, 
Cambridge University Press
1993
*/

/****************************************************************************/

static double sqrarg;

void nrerror(char error_text[]);
double *vector(int nl, int nh);
void free_vector(double *v, int nl);
double **InitDoubleArray(int nbRows, int nbCols);
int **InitIntArray(int nbRows, int nbCols);
void FreeDoubleArray(double **anArray, int nbRows);
void FreeIntArray(int **anArray, int nbRows);
double pythag(double a, double b);
void PrintMatrix(double **mat, int nbRows, int nbCols);
void PrintVector(double vect[], int nbElems);
void svdcmp(double **a, int m, int n, double *w, double **v);
int  PseudoInverse(double **A,int N, double **B);

/****************************************************************************/

static double at,bt,ct;
static double maxarg1,maxarg2;

const int nbCoords = 3;

int debug2 = 1; /* 0 to avoid printing */

#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg * sqrarg)

#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
	(maxarg1) : (maxarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/****************************************************************************/

void nrerror(char error_text[]) 
{

printf("%s\n",error_text);
printf("...now exiting to system...\n");
exit(1);
}

/****************************************************************************/

double *vector(int nl, int nh) /*allocation of memory*/
{
double *v;
v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
if (!v) nrerror("allocation failure in vector()");
return v-nl;
}

/****************************************************************************/

void free_vector(double *v, int nl) /*deallocation*/
{
free((char*) (v+nl));
}

/****************************************************************************/

double pythag(double a, double b)
{
double absa, absb;

	absa = fabs(a);
	absb = fabs(b);
	if (absa > absb) return absa * sqrt(1.0 + SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa/absb)));
}

/****************************************************************************/

double **InitDoubleArray(int nbRows, int nbCols)
{
/* 
Allocate memory for a 2D array of double 
*/
	int i;
	double **anArray;

	anArray = (double **)malloc(nbRows * sizeof(double *));
	if (!anArray) nrerror("allocation failure in InitDoubleArray");
	for(i=0;i<nbRows;i++){
		anArray[i]=(double *)malloc(nbCols * sizeof(double));
		if (!anArray[i]) nrerror("allocation failure in InitDoubleArray");
	}
	return anArray;
}

/****************************************************************************/

int **InitIntArray(int nbRows, int nbCols)
{
/* 
Allocate memory for a 2D array of int 
*/
int i;
int **anArray;

	anArray = (int **)malloc(nbRows * sizeof(int *));
	if (!anArray) nrerror("allocation failure in InitIntArray");
	for(i=0;i<nbRows;i++){
		anArray[i]=(int *)malloc(nbCols * sizeof(int));
		if (!anArray[i]) nrerror("allocation failure in InitIntArray");
	}
	return anArray;
}

/****************************************************************************/

void FreeDoubleArray(double **anArray, int nbRows)
{
int i;

	for(i=0;i<nbRows;i++){
		free(anArray[i]);
	}
}

/****************************************************************************/

void FreeIntArray(int **anArray, int nbRows)
{
int i;

	for(i=0;i<nbRows;i++){
		free(anArray[i]);
	}
}

/****************************************************************************/

void PrintMatrix(double **mat, int nbRows, int nbCols)
{
int i, j;
	for (i=0; i<nbRows;i++){
		for (j=0; j<nbCols; j++){
			printf("%.5f ", mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

/****************************************************************************/

void PrintVector(double vect[], int nbElems)
{
int i;
	for (i=0; i<nbElems;i++){
		printf("%.2f ", vect[i]);
	}
	printf("\n");
}

	
/****************************************************************************/


void svdcmp(double **a, int m, int n, double w[], double **v) 
/*
Given a matrix a with dimension m x n, get singular value decomposition, 
a=U.W.(V)t. Matrix U replaces a in output. Diagonal matrix W is provided as vector w.
m must be larger than n or equal to it. If not, user must pad with zeroes.
*/

{
int flag,i,its,j,jj,k,l,nm;
double c,f,h,s,x,y,z;
double anorm=0.0,g=0.0,scale=0.0;
double *rv1,*vector();
double myNb = 1;

  if (m < n) nrerror("SVDCMP: You must augment A with extra zero rows");
  rv1=vector((int)(1),n);
  for (i=1;i<=n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m) {
      for (k=i;k<=m;k++) scale += fabs(a[k-1][i-1]);
      if (scale) {
	for (k=i;k<=m;k++) {
	  a[k-1][i-1] /= scale;
	  s += a[k-1][i-1]*a[k-1][i-1];
	}
	f=a[i-1][i-1];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i-1][i-1]=f-g;
	if (i != n) {
	  for (j=l;j<=n;j++) {
	    for (s=0.0,k=i;k<=m;k++) s += a[k-1][i-1]*a[k-1][j-1];
	    f=s/h;
	    for (k=i;k<=m;k++) a[k-1][j-1] += f*a[k-1][i-1];
	  }
	}
	for (k=i;k<=m;k++) a[k-1][i-1] *= scale;
      }
    }
    w[i-1]=scale*g;
    g=s=scale=0.0;
    if (i <= m && i != n) {
      for (k=l;k<=n;k++) scale += fabs(a[i-1][k-1]);
      if (scale) {
	for (k=l;k<=n;k++) {
	  a[i-1][k-1] /= scale;
	  s += a[i-1][k-1]*a[i-1][k-1];
	}
	f=a[i-1][l-1];
	g = -SIGN(sqrt(s),f);
	h=f*g-s;
	a[i-1][l-1]=f-g;
	for (k=l;k<=n;k++) rv1[k]=a[i-1][k-1]/h;
/*	if (i != m) { */
	  for (j=l;j<=m;j++) {
	    for (s=0.0,k=l;k<=n;k++) s += a[j-1][k-1]*a[i-1][k-1];
	    for (k=l;k<=n;k++) a[j-1][k-1] += s*rv1[k];
	  }
/*	}*/
	for (k=l;k<=n;k++) a[i-1][k-1] *= scale;
      }
    }
    anorm=MAX(anorm,(fabs(w[i-1])+fabs(rv1[i])));
  }
  for (i=n;i>=1;i--) {
    if (i < n) {
      if (g) {
	for (j=l;j<=n;j++)
	  v[j-1][i-1]=(a[i-1][j-1]/a[i-1][l-1])/g;
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=l;k<=n;k++) s += a[i-1][k-1]*v[k-1][j-1];
	  for (k=l;k<=n;k++) v[k-1][j-1] += s*v[k-1][i-1];
	}
      }
      for (j=l;j<=n;j++) v[i-1][j-1]=v[j-1][i-1]=0.0;
    }
    v[i-1][i-1]=1.0;
    g=rv1[i];
    l=i;
  }
  for (i=n;i>=1;i--) {
    l=i+1;
    g=w[i-1];
    if (i < n)
      for (j=l;j<=n;j++) a[i-1][j-1]=0.0;
    if (g) {
      g=1.0/g;
      if (i != n) {
	for (j=l;j<=n;j++) {
	  for (s=0.0,k=l;k<=m;k++) s += a[k-1][i-1]*a[k-1][j-1];
	  f=(s/a[i-1][i-1])*g;
	  for (k=i;k<=m;k++) a[k-1][j-1] += f*a[k-1][i-1];
	}
      }
      for (j=i;j<=m;j++) a[j-1][i-1] *= g;
    }
    else {
    for (j=i;j<=m;j++) a[j-1][i-1]=0.0;
    }
    ++a[i-1][i-1];
  }
  for (k=n;k>=1;k--) {
    for (its=1;its<=30;its++) {
      flag=1;
      for (l=k;l>=1;l--) {
	nm=l-1;
	if (fabs(rv1[l])+anorm == anorm) {
	  flag=0;
	  break;
	}
	if (fabs(w[nm])+anorm == anorm) break;
      }
      if (flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f=s*rv1[i];
	  if (fabs(f)+anorm != anorm) {
	    g=w[i-1];
	    h=pythag(f,g);
	    w[i-1]=h;
	    h=1.0/h;
	    c=g*h;
	    s=(-f*h);
	    for (j=1;j<=m;j++) {
	      y=a[j-1][nm-1];
	      z=a[j-1][i-1];
	      a[j-1][nm-1]=y*c+z*s;
	      a[j-1][i-1]=z*c-y*s;
	    }
	  }
	}
      }
      z=w[k-1];
      if (l == k) {
	if (z < 0.0) {
	  w[k-1] = -z;
	  for (j=1;j<=n;j++) v[j-1][k-1]=(-v[j-1][k-1]);
	}
	break;
      }
      if (its == 50) nrerror("No convergence in 50 SVDCMP iterations");
      x=w[l-1];
      nm=k-1;
      y=w[nm-1];
      g=rv1[nm];
      h=rv1[k];
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g=pythag(f,1.0);
      f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
      c=s=1.0;
      for (j=l;j<=nm;j++) {
	i=j+1;
	g=rv1[i];
	y=w[i-1];
	h=s*g;
	g=c*g;
	z=pythag(f,h);
	rv1[j]=z;
	c=f/z;
	s=h/z;
	f=x*c+g*s;
	g=g*c-x*s;
	h=y*s;
	y=y*c;
	for (jj=1;jj<=n;jj++) {
	  x=v[jj-1][j-1];
	  z=v[jj-1][i-1];
	  v[jj-1][j-1]=x*c+z*s;
	  v[jj-1][i-1]=z*c-x*s;
	}
	z=pythag(f,h);
	w[j-1]=z;
	if (z) {
	  z=1.0/z;
	  c=f*z;
	  s=h*z;
	}
	f=(c*g)+(s*y);
	x=(c*y)-(s*g);
	for (jj=1;jj<=m;jj++) {
	  y=a[jj-1][j-1];
	  z=a[jj-1][i-1];
	  a[jj-1][j-1]=y*c+z*s;
	  a[jj-1][i-1]=z*c-y*s;
	}
      }
      rv1[l]=0.0;
      rv1[k]=f;
      w[k-1]=x;
    }
  }
  free_vector(rv1,(int)(1));
}

/****************************************************************************/

int  PseudoInverse(double **A, int N, double **B) 
/*
B returns the pseudoinverse of A (A with dimension N rows x 3 columns).
It is the matrix v.[diag(1/wi)].(u)t (cf. svdcmp())
Function returns True if B has maximum rank, False otherwise
*/
{
  void svdcmp();
  double **V,temp[3][3];
  double W[3];
  double WMAX;
  double TOL = 0.01;
  int i,j,k;
  int isMaxRank = 1;/* stays true if no singular value under tolerance level

  /*allocations*/
  V = InitDoubleArray(nbCoords, nbCoords);
  
  /*Singular value decomposition*/
  PrintMatrix(A, N, nbCoords);
  svdcmp(A,N,nbCoords,W,V);

  if(debug2){
  	PrintMatrix(A, N, nbCoords);
  	PrintMatrix(V, nbCoords, nbCoords);
  	PrintVector(W, nbCoords);
  }
  /*Getting largest singular value*/
  WMAX=0.0;
  for (i=0;i<nbCoords;i++){
      if (W[i]>WMAX) WMAX=W[i];
  }

  /*Checking for signular values smaller than TOL times the largest one*/
  for (i=0;i<nbCoords;i++){
      if (W[i]<TOL*WMAX){
      	W[i]=0;
      	isMaxRank = 0;
      	return isMaxRank;
      }
  }
  
  if(isMaxRank){
  	/*Computing B*/     
  	for (i=0;i<3;i++){
      	for (j=0;j<3;j++){
	  		temp[i][j]=V[i][j]/W[j];
		}
    }
  	for (i=0;i<3;i++){
		for (j=0;j<N;j++){
			B[i][j]=0.0;
	  		for (k=0;k<3;k++){
	  			B[i][j]+=temp[i][k]*A[j][k];
			}
    	}
	}
	if(debug2) PrintMatrix(B, nbCoords, N);

	/*desallocations*/
	for (i=0;i<nbCoords;i++) free(V[i]);

	return isMaxRank;
  }
}

/****************************************************************************/

