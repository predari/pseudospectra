/*******************************************************************************
*  Copyright (C) 2009-2012 Intel Corporation. All Rights Reserved.
*  The information and material ("Material") provided below is owned by Intel
*  Corporation or its suppliers or licensors, and title to such Material remains
*  with Intel Corporation or its suppliers or licensors. The Material contains
*  proprietary information of Intel or its suppliers and licensors. The Material
*  is protected by worldwide copyright laws and treaty provisions. No part of
*  the Material may be copied, reproduced, published, uploaded, posted,
*  transmitted, or distributed in any way without Intel's prior express written
*  permission. No license under any patent, copyright or other intellectual
*  property rights in the Material is granted to or conferred upon you, either
*  expressly, by implication, inducement, estoppel or otherwise. Any license
*  under such intellectual property rights must be express and approved by Intel
*  in writing.
*
* 
* Same as my_conec_exam.c but passing in conrec instead of double *xaxis and double *
* yaxis double _Complex z.
* 
* 
*  A static example with an grcar array A of dimensions  (M,N) = (8,8).
*  The main data structure in the implementation is a single pointer.
*  Compile with: gcc -o con my_conrec_exam_2.c paulslib.c -lm bitmaplib.c -llapacke -llapack -lrefblas -L/usr/lib/gcc/x86_64-linux-gnu/4.4 -lgfortran
*  Should be compared with the matlab function @tesla ~/predari/Documents/thesis/pseudospectra/grcar_example.m
* 
********************************************************************************
*/
/*
   LAPACKE_zgesvd Example.
   =======================

   Program computes the singular value decomposition of a general
   rectangular complex matrix A:

   (  5.91, -5.69) (  7.09,  2.72) (  7.78, -4.06) ( -0.79, -7.21)
   ( -3.15, -4.08) ( -1.89,  3.27) (  4.57, -2.07) ( -3.88, -3.30)
   ( -4.89,  4.20) (  4.10, -6.70) (  3.28, -3.84) (  3.84,  1.19)

   Description.
   ============

   The routine computes the singular value decomposition (SVD) of a complex
   m-by-n matrix A, optionally computing the left and/or right singular
   vectors. The SVD is written as

   A = U*SIGMA*VH

   where SIGMA is an m-by-n matrix which is zero except for its min(m,n)
   diagonal elements, U is an m-by-m unitary matrix and VH (V conjugate
   transposed) is an n-by-n unitary matrix. The diagonal elements of SIGMA
   are the singular values of A; they are real and non-negative, and are
   returned in descending order. The first min(m, n) columns of U and V are
   the left and right singular vectors of A.

   Note that the routine returns VH, not V.

   Example Program Results.
   ========================

 LAPACKE_zgesvd (row-major, high-level) Example Program Results

 Singular values
  17.63  11.61   6.78

 Left singular vectors (stored columnwise)
 ( -0.86,  0.00) (  0.40,  0.00) (  0.32,  0.00)
 ( -0.35,  0.13) ( -0.24, -0.21) ( -0.63,  0.60)
 (  0.15,  0.32) (  0.61,  0.61) ( -0.36,  0.10)

 Right singular vectors (stored rowwise)
 ( -0.22,  0.51) ( -0.37, -0.32) ( -0.53,  0.11) (  0.15,  0.38)
 (  0.31,  0.31) (  0.09, -0.57) (  0.18, -0.39) (  0.38, -0.39)
 (  0.53,  0.24) (  0.49,  0.28) ( -0.47, -0.25) ( -0.15,  0.19)
*/

#include "math.h"
#include "paulslib.h"
#include "bitmaplib.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <lapacke.h>


#define min(a,b) ((a)>(b)?(b):(a))

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda );
extern void print_rmatrix( char* desc, lapack_int m, lapack_int n, double* a, lapack_int lda );
void prtdat(int nx, int ny, double * u1, char *fnam);


/****
 * 
 *  THESE ARE EXTRA LINES FOR DRAWING THE PSEUDOSPECTRA 
 * 
 * ****/
/* Scaling the image for better resolution */
//#define SCALE 5
/* Two dimensional array of data for conrec*/
double **data;

/* Arrays for the x axis and y axis coordinates */
/****
 * XAXIS AND YAXIS CAN BE BROUGHT FROM THE REAL PART AND IMAGINARY 
 * PART OF z ARRAY CORRESPONDINGLY THAT STORES COMPLEX NUMBERS. 
 * **/
double _Complex *z;

/* Array for the contour levels */
/***
 * FOR THE TIME BEING ONLY ONE CORRESPONDING TO e=0.1
 * **/
#define NCONTOUR 1
double contours[NCONTOUR];

/* Image on which the contours will be drawn, see bitmaplib.c */
BITMAP4 *image;

/* Prototype for CONREC and the line drawing function */
void CONREC(double **,int,int,int,int,double _Complex *,int,double *,
       void (*drawline)(double,double,double,double,double));
void drawline(double,double,double,double,double);

/* Debugging - count the number of line segments drawn */
int vectorsdrawn = 0;
 FILE *fp;

/* Parameters */
#define M 50
#define N 50
#define NGRID 50
#define LDA N
#define LDU M
#define LDVT N

/***
 * I WOULD LIKE XMAX TO BE DEFINED AS DOUBLE
 ***/
#define XMAX      4                /* maximum boundary of x-axis of the domain */
#define XMIN     -2                 /*minimum*/
#define YMAX      4                 /* maximum boundary of y-axis of the domain */
#define YMIN      -4                    /*minimum*/

/* Main program */
int main(){
        /* not Locals */
        lapack_complex_double *a, *temp, * u, *vt;
        lapack_int m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT, info;
		
        /* Local arrays */
		//void prtdat();
		double *s;
        double *superb;
	int svd_count=0;
		int i, j ,ix ,iy, index, ii, jj ;    
		double  x_min,
				x_max,
				y_min,
				y_max,
				stepx,						/* step size for finding gridpoints coordinates in x and y dimension.*/
				stepy;
		double e=0.1;        
		/* Array used for the ploting of
		* grid, as an input to the 
		* draw_pseudospectra function. */
		double *plot;
		//double plot[n][n]; 
		COLOUR colour;
		BITMAP4 col,grey = {128,128,128,0};
		
		       
        /* Memory alocations*/
		temp = malloc((lda*m)*sizeof(lapack_complex_double));
		a = malloc((lda*m)*sizeof(lapack_complex_double));
		u = malloc((ldu*m)*sizeof(lapack_complex_double));
		vt = malloc((ldvt*n)*sizeof(lapack_complex_double));
		s = malloc(m*sizeof(double));
		superb = malloc(min(m,n)*sizeof(double));
		plot = malloc((NGRID*NGRID)*sizeof(double));
		z = malloc((NGRID*NGRID)*sizeof(double _Complex));
		
	
		//allocating the 2D array data.
	  if ((data = malloc(NGRID*sizeof(double *))) == NULL) {
      fprintf(stderr,"Failed to malloc space for the data\n");
      exit(-1);
   }
   for (i=0;i<NGRID;i++) {
      if ((data[i] = malloc(NGRID*sizeof(double))) == NULL) {
         fprintf(stderr,"Failed to malloc space for the data\n");
         exit(-1);
      }
   }
   for (i=0;i<NGRID;i++){
      for (j=0;j<NGRID;j++){
         data[i][j] = 0;
	 //   printf("%f\t",data[i][j]);
	 }
		  }
		
/*
		printf("-------------------------------------------------\n");
		printf("        ---------------------------------           \n");
		printf ("Starting Computing Pseudopsecta of grcar Matrix\n");
		printf("Give the doundaries of the 2-dimenional domain\n");
		printf("Insert the minimum value of x-axis\n");
		clearerr(stdin);
		scanf("%lf",&x_min);
		//getchar();
		printf("Insert the maximum value of x-axis\n");
		scanf("%lf",&x_max);
		printf("Insert the minimum value of y-axis\n");
		scanf("%lf",&y_min);
		printf("Insert the maximun value of y-axis\n");
		scanf("%lf",&y_max);
		//printf("Give the grid size you want:\n");
		//scanf("%d",&n);
*/ 	  
		/*if (x_min==0.0)*/  x_min=XMIN;
		/*if (x_max==0.0)*/  x_max=XMAX;
		/*if (y_min==0.0)*/  y_min=YMIN;
		/*if (y_max==0.0)*/  y_max=YMAX;
	
		/* Initialize grid */
		printf("The size of the domain is: X=[%f-%f]  Y=[%f-%f] \n",x_min,x_max,y_min,y_max);
	  
		stepx=(double)abs(x_max-x_min)/(NGRID-1);
		stepy=(double)abs(y_max-y_min)/(NGRID-1);
		printf("To stepx einai %f\n",stepx);
		printf("To stepy einai %f\n",stepy);	
	   

	
	   for (i =0; i <NGRID*NGRID; i++){
			z[i]=x_min+(i/n * stepx)+(y_min + (i%n * stepy))*I;
		 // z[i]=lapack_make_complex_double( i/n,i%n); just for testing
		//**	printf( " (%6.2f,%6.2f)", lapack_complex_double_real(z[i]), lapack_complex_double_imag(z[i]) );
		}

	   memset(temp,0,(lda*m)*sizeof(*temp));
	   memset(a,0,(lda*m)*sizeof(*a));
	   memset(u,0,(ldu*m)*sizeof(*u));
	   memset(vt,0,(ldvt*m)*sizeof(*vt));
		
		
	   j=0;
	   for (i = 0; i < lda*m ; i=i+n ){	
			if(i==0){
				a[i]=lapack_make_complex_double( 1,0);
				a[i+1]=lapack_make_complex_double( 1,0);
				a[i+2]=lapack_make_complex_double( 1,0);
				a[i+3]=lapack_make_complex_double( 1,0);
			}
			else if(i == (n-3)*n ){
				a[i+j]=lapack_make_complex_double( -1,0);
				a[i+(j+1)]=lapack_make_complex_double( 1,0);
				a[i+(j+2)]=lapack_make_complex_double( 1,0);
				a[i+(j+3)]=lapack_make_complex_double( 1,0);
				j++;
			}
			else if(i == (n-2)*n ){
				a[i+j]=lapack_make_complex_double( -1,0);
				a[i+(j+1)]=lapack_make_complex_double( 1,0);
				a[i+(j+2)]=lapack_make_complex_double( 1,0);
				j++;
			}
			else if(i == (n-1)*n ){
				a[i+j]=lapack_make_complex_double( -1,0);
				a[i+(j+1)]=lapack_make_complex_double( 1,0);
				j++;
			}
			else{
				a[i+j]=lapack_make_complex_double( -1,0);
				a[i+(j+1)]=lapack_make_complex_double( 1,0);
				a[i+(j+2)]=lapack_make_complex_double( 1,0);
				a[i+(j+3)]=lapack_make_complex_double( 1,0);
				a[i+(j+4)]=lapack_make_complex_double( 1,0);
				j++;
			}
		} 

		//print_matrix("Entry Matrix A", m, n, a, lda );
		for (iy = 0; iy < NGRID*NGRID; iy++){   
			 //printf("temp size %d, a size %d",(lda*m)*sizeof(*temp),(lda*m)*sizeof(*a));
			memcpy(temp, a ,(lda*m)*sizeof(*temp));
			 //~ print_matrix( "Entry Matrix Temp just after memcopy", m, n, temp, lda );
			 //~ print_matrix( "Entry Matrix A just after memcopy", m, n, a, lda );
			// printf( "To  z[%d](%6.4f,%6.4f)\n",iy,lapack_complex_double_real(z[iy]),lapack_complex_double_imag(z[iy]) );
			for (i = 0; i < lda*m ; i=i+(n+1)){	
				//~ printf("%d",i);
				//~ printf( "To  a[%d](%6.2f,%6.2f)\t",i, lapack_complex_double_real(a[i]), lapack_complex_double_imag(a[i]) );
				//~ printf( "To  z[%d](%6.2f,%6.2f)\n",iy,lapack_complex_double_real(z[iy]),lapack_complex_double_imag(z[iy]) );
				
				temp[i]=a[i]-z[iy];
				//~ temp[index] = lapack_make_complex_double(lapack_complex_double_real(a[index])-lapack_complex_double_real(z[iy]),  lapack_complex_double_imag(a[index])-lapack_complex_double_imag(z[iy])    );
				//~ printf( " temp[%d](%6.2f,%6.2f)", i,lapack_complex_double_real(temp[i]), lapack_complex_double_imag(temp[i]) );
				//~ printf( "\n");
			}
			//printf("GRCAR MATRIX AFTER SUBSTRACTION (%d,%d)\n",iy/n,iy%n);
			//~ print_matrix( "Entry Matrix Temp just before", m, n, temp, lda );
	
			/* Executable statements */
			//~ print_matrix( "AT THE BEGINING OF THE FOR LOOP", m, n, a, lda );
			printf( "LAPACKE_zgesvd (row-major, high-level) Example Program Results(%d,%d)\n",iy/NGRID,iy%NGRID);
			/* Compute SVD */
			info = LAPACKE_zgesvd( LAPACK_ROW_MAJOR, 'N', 'N', m, n, temp, lda, s, NULL, ldu, NULL, ldvt, superb );
			svd_count++;
			//~ 
			//~ print_matrix( "IN THE MIDDLE OF THE FOR LOOP", m, n, a, lda );
			//~ print_matrix( "IN THE MIDDLE OF THE FOR LOOP-TEMP", m, n, temp, lda );
			/* Check for convergence */
			if( info > 0 ) {
				printf( "The algorithm computing SVD failed to converge.\n" );
				exit( 1 );
			}
			/* Print singular values */
			if( info == 0){
//				printf("Solution\n");	
				for ( i= 0; i< m; i++ ) {
//					printf(" s[ %d ] = %f\n", i, s[ i ] );
				}
			}
			
			if(s[m-1] <= e){
				printf("THIS ELEMENT BELONGS TO PSEUDOSPECTRA (%d,%d):%6.10f\n",(iy/NGRID+1),(iy%NGRID+1),s[m-1]);
				/*to index tis parapanw ektupwshs anaferetai sto index tou antistoixou mhtrwou apo thn synarthsh ths matlab grcar_example.m*/
				//~ plot[iy/n][iy%n]=s[m-1];
				plot[iy]=s[m-1];
			 }
			 //~ else   plot[iy/n][iy%n]=0;
				else plot[iy]=0;
		
	
		
	//~ print_rmatrix( "Singular values", 1, m, s, 1 );
	
	/* Print left singular vectors */
	// print_matrix( "Left singular vectors (stored columnwise)", m, m, u, ldu );
	/* Print right singular vectors */
	// print_matrix( "Right singular vectors (stored rowwise)", m, n, vt, ldvt );
		}
		
		
				prtdat(NGRID, NGRID, plot, "svd_with_disks.data");
		printf("Total number of svd evaluations in the %d,%d grid is:\t %d\n",NGRID,NGRID,svd_count);
		
		
		//giving values to data from plot
		for (i = 0; i<NGRID*NGRID; i++) data[i/NGRID][i%NGRID] = plot[i];
	   /////////////////
    BITMAP4 black = {0,0,0,0};
    Draw_Line(image,NGRID,NGRID,x_min,y_min,x_max,y_min,black);
   //////////////////	
		//~ contours[0] = 0.1;
		//~ contours[1] = 0.01;
		//~ contours[2] = 0.001;
		//~ contours[3] = 0.0001;
		//~ contours[4] = 0.00001;
		if ((image = Create_Bitmap(NGRID,NGRID)) == NULL) {
      fprintf(stderr,"Malloc of bitmap failed\n");
      exit(-1);
   }
   Erase_Bitmap(image,NGRID,NGRID,grey); /* Not strictly necessary */
   for (j=0;j<NGRID;j++) {
      for (i=0;i<NGRID;i++) {
         colour = GetColour(data[i][j],0,0.1,1);                   /////////////////////////////////////////////////////
         col.r = colour.r * 255;
         col.g = colour.g * 255;
         col.b = colour.b * 255;
         Draw_Pixel(image,NGRID,NGRID,(double)i,(double)j,col);
      }
   }


   /* Finally do the contouring */
   CONREC(data,0,NGRID-1,0,NGRID-1,
      z,NCONTOUR,contours,drawline);
   fprintf(stderr,"Drew %d vectors\n",vectorsdrawn);

   /* 
      Write the image as a TGA file 
      See bitmaplib.c for more details, or write "image"
      in your own prefered format.
   */
   if ((fp = fopen("image.tga","w")) == NULL) {
      fprintf(stderr,"Failed to open output image\n");
      exit(-1);
   }
   Write_Bitmap(fp,image,NGRID,NGRID,12);
   fclose(fp);

		
		
		
		exit(0);
} /* End of LAPACKE_zgesvd Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda ) {
        lapack_int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ )  
                        printf( " (%6.2f,%6.2f)", lapack_complex_double_real(a[i*lda+j]), lapack_complex_double_imag(a[i*lda+j]) );
                printf( "\n" );
        }
}

/* Auxiliary routine: printing a vector of integers */
void print_int_vector( char* desc, lapack_int n, lapack_int* a ) {
        lapack_int j;
        printf( "\n %s\n", desc );
        for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
        printf( "\n" );
}


 void prtdat(int nx, int ny, double *u1, char *fnam) {
 int ix, iy;
// FILE *fp;
// 
 fp = fopen(fnam, "w");
 //fprintf(fp, "%d %d\n", nx, ny);
   for (ix = 0; ix < nx*ny; ix++) {
     fprintf(fp, "%6.10f", u1[ix]);
     if (ix % nx == nx-1) 
       fprintf(fp, "\n");
     else
       fprintf(fp, " ");
     }
 fclose(fp);
 }


void drawline(double x1,double y1,double x2,double y2,double z)
{
   BITMAP4 black = {0,0,0,0};

   //~ if (x1 < 0 || x1 >= M || x2 < 0 || x2 > N)
      //~ fprintf(stderr,"Shouldn't get here, x out of bounds: %g %g\n",x1,x2);
   //~ if (y1 < 0 || y1 >= N || y2 < 0 || y2 > N)
      //~ fprintf(stderr,"Shouldn't get here, y out of bounds: %g %g\n",y1,y2);
  
  
   Draw_Line(image,NGRID,NGRID,(int)x1,(int)y1,(int)x2,(int)y2,black);

   
   vectorsdrawn++;
}

/*
   Derivation from CONREC
   d               ! matrix of data to contour
   ilb,iub,jlb,jub ! index bounds of data matrix
   x               ! data matrix column coordinates
   y               ! data matrix row coordinates
   nc              ! number of contour levels
   z               ! contour levels in increasing order
*/
void CONREC(double **d,int ilb,int iub,int jlb,int jub,
   double _Complex *zcomplex,int nc,double *z,
   void (*ConrecLine)(double,double,double,double,double))
{
#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])

   int m1,m2,m3,case_value;
   double dmin,dmax,x1,x2,y1,y2;
   int i,j,k,m;
   double h[5];
   int sh[5];
   double xh[5],yh[5];
   int im[4] = {0,1,1,0},jm[4]={0,0,1,1};
   int castab[3][3][3] = {
     { {0,0,8},{0,2,5},{7,6,9} },
     { {0,3,4},{1,3,1},{4,3,0} },
     { {9,6,7},{5,2,0},{8,0,0} }
   };
   double temp1,temp2;

   for (j=(jub-1);j>=jlb;j--) {
      for (i=ilb;i<=iub-1;i++) {
         temp1 = MIN(d[i][j],d[i][j+1]);
         temp2 = MIN(d[i+1][j],d[i+1][j+1]);
         dmin  = MIN(temp1,temp2);
         temp1 = MAX(d[i][j],d[i][j+1]);
         temp2 = MAX(d[i+1][j],d[i+1][j+1]);
         dmax  = MAX(temp1,temp2);
         if (dmax < z[0] || dmin > z[nc-1])
            continue;
         for (k=0;k<nc;k++) {
            if (z[k] < dmin || z[k] > dmax)
               continue;
            for (m=4;m>=0;m--) {
               if (m > 0) {
                  h[m]  = d[i+im[m-1]][j+jm[m-1]]-z[k];
                  xh[m] =/*x*/ creal(zcomplex[i+im[m-1]]);
                  yh[m] = cimag(zcomplex[j+jm[m-1]]);
               } else {
                  h[0]  = 0.25 * (h[1]+h[2]+h[3]+h[4]);
                  xh[0] = 0.50 * (creal(zcomplex[i]+zcomplex[i+1]));/*x*/
                  yh[0] = 0.50 * (cimag(zcomplex[j]+zcomplex[j+1]));
               }
               if (h[m] > 0.0)
                  sh[m] = 1;
               else if (h[m] < 0.0)
                  sh[m] = -1;
               else
                  sh[m] = 0;
            }

            /*
               Note: at this stage the relative heights of the corners and the
               centre are in the h array, and the corresponding coordinates are
               in the xh and yh arrays. The centre of the box is indexed by 0
               and the 4 corners by 1 to 4 as shown below.
               Each triangle is then indexed by the parameter m, and the 3
               vertices of each triangle are indexed by parameters m1,m2,and m3.
               It is assumed that the centre of the box is always vertex 2
               though this isimportant only when all 3 vertices lie exactly on
               the same contour level, in which case only the side of the box
               is drawn.
                  vertex 4 +-------------------+ vertex 3
                           | \               / |
                           |   \    m-3    /   |
                           |     \       /     |
                           |       \   /       |
                           |  m=2    X   m=2   |       the centre is vertex 0
                           |       /   \       |
                           |     /       \     |
                           |   /    m=1    \   |
                           | /               \ |
                  vertex 1 +-------------------+ vertex 2
            */
            /* Scan each triangle in the box */
            for (m=1;m<=4;m++) {
               m1 = m;
               m2 = 0;
               if (m != 4)
                  m3 = m + 1;
               else
                  m3 = 1;
               if ((case_value = castab[sh[m1]+1][sh[m2]+1][sh[m3]+1]) == 0)
                  continue;
               switch (case_value) {
               case 1: /* Line between vertices 1 and 2 */
                   x1 = xh[m1];
                   y1 = yh[m1];
                   x2 = xh[m2];
                   y2 = yh[m2];
                   break;
               case 2: /* Line between vertices 2 and 3 */
                   x1 = xh[m2];
                   y1 = yh[m2];
                   x2 = xh[m3];
                   y2 = yh[m3];
                   break;
               case 3: /* Line between vertices 3 and 1 */
                   x1 = xh[m3];
                   y1 = yh[m3];
                   x2 = xh[m1];
                   y2 = yh[m1];
                   break;
               case 4: /* Line between vertex 1 and side 2-3 */
                   x1 = xh[m1];
                   y1 = yh[m1];
                   x2 = xsect(m2,m3);
                   y2 = ysect(m2,m3);
                   break;
               case 5: /* Line between vertex 2 and side 3-1 */
                   x1 = xh[m2];
                   y1 = yh[m2];
                   x2 = xsect(m3,m1);
                   y2 = ysect(m3,m1);
                   break;
               case 6: /* Line between vertex 3 and side 1-2 */
                   x1 = xh[m1];
                   y1 = yh[m1];
                   x2 = xsect(m1,m2);
                   y2 = ysect(m1,m2);
                   break;
               case 7: /* Line between sides 1-2 and 2-3 */
                   x1 = xsect(m1,m2);
                   y1 = ysect(m1,m2);
                   x2 = xsect(m2,m3);
                   y2 = ysect(m2,m3);
                   break;
               case 8: /* Line between sides 2-3 and 3-1 */
                   x1 = xsect(m2,m3);
                   y1 = ysect(m2,m3);
                   x2 = xsect(m3,m1);
                   y2 = ysect(m3,m1);
                   break;
               case 9: /* Line between sides 3-1 and 1-2 */
                   x1 = xsect(m3,m1);
                   y1 = ysect(m3,m1);
                   x2 = xsect(m1,m2);
                   y2 = ysect(m1,m2);
                   break;
               default:
                   break;
               }

               /* Finally draw the line */
               ConrecLine(x1,y1,x2,y2,z[k]);
            } /* m */
         } /* k - contour */
      } /* i */
   } /* j */
}
