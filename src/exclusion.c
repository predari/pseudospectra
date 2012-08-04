/**
 * 27/07/2012
 * Trying to imlement exclusion disks.
 * 
 * Same as my_conrec_exam.c but passing in conrec instead of double *xaxis and double *yaxis 
 * double _Complex z.
 * 
 * A static example with an grcar array A of dimensions  (M,N) = (8,8).
 * The main data structure in the implementation is a single pointer.
 * Should be compared with the matlab function @tesla ~/predari/Documents/thesis/pseudospectra/grcar_example.m
 * 
 * Compile with: 
 * 	gcc -o ex my_exclusion2.c paulslib.c -lm bitmaplib.c -llapacke -llapack -lrefblas -L/usr/lib/gcc/x86_64-linux-gnu/4.4 -lgfortran
 */

#include "../lib/paulslib.h"
#include "../lib/bitmaplib.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <lapacke.h>

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
#define SCALE 15
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
/*A (m by n)*/
#define M 32
#define N 32
#define NGRID 50
#define LDA N
#define LDU M
#define LDVT N

/***
 * I WOULD LIKE XMAX TO BE DEFINED AS DOUBLE
 ***/
#define XMAX      4 /*4                /* maximum boundary of x-axis of the domain */
#define XMIN     -1 /*-2                 /*minimum*/
#define YMAX      4 /*4                 /* maximum boundary of y-axis of the domain */
#define YMIN      -4 /*-4                    /*minimum*/

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
		int i, j ,ix ,iy, index, ii, jj, k ;    
		double  x_min,
				x_max,
				y_min,
				y_max,
				stepx,						/* step size for finding gridpoints coordinates in x and y dimension.*/
				stepy;
		double e=0.1;  
		double r;   /*r is the radious of exclusion disks */
		int i_max, j_max, start_point; 
		double gamma;     
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
		superb = malloc(MIN(m,n)*sizeof(double));
		plot = malloc((NGRID*NGRID)*sizeof(double));
		z = malloc((NGRID*NGRID)*sizeof(double _Complex));
		
		
		//allocating the 2D array data.
	  if ((data = malloc(SCALE*NGRID*sizeof(double *))) == NULL) {
      fprintf(stderr,"Failed to malloc space for the data\n");
      exit(-1);
   }
   for (i=0;i<SCALE*NGRID;i++) {
      if ((data[i] = malloc(SCALE*NGRID*sizeof(double))) == NULL) {
         fprintf(stderr,"Failed to malloc space for the data\n");
         exit(-1);
      }
   }
   for (i=0;i<SCALE*NGRID;i++){
      for (j=0;j<SCALE*NGRID;j++){
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
		    printf("%f+%fi\n",creal(z[i]),cimag(z[i]));
		}
       
       //memset(plot,-1,(NGRID*NGRID)*sizeof(double));
	    for (i =0; i <NGRID*NGRID; i++){
		    plot[i]=-1;
		//	printf("%f\t",plot[i]);
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

		print_matrix("Entry Matrix A", m, n, a, lda );
		for (iy = 0; iy < NGRID*NGRID; iy++){   
			

			
			if(plot[iy]==0) continue;
			
			
			 //printf("temp size %d, a size %d",(lda*m)*sizeof(*temp),(lda*m)*sizeof(*a));
			memcpy(temp, a ,(lda*m)*sizeof(*temp));
			
			for (i = 0; i < lda*m ; i=i+(n+1)){	

				
				temp[i]=a[i]-z[iy];

			}
			printf( "LAPACKE_zgesvd (row-major, high-level) Example Program Results(%d,%d)\n",iy/NGRID,iy%NGRID);
			/* Compute SVD */
			info = LAPACKE_zgesvd( LAPACK_ROW_MAJOR, 'N', 'N', m, n, temp, lda, s, NULL, ldu, NULL, ldvt, superb );
			svd_count++;

			/* Check for convergence */
			if( info > 0 ) {
				printf( "The algorithm computing SVD failed to converge.\n" );
				exit( 1 );
			}
			/* Print singular values */
			if( info == 0){
	
				for ( i= 0; i< m; i++ ) {

				}
			}
			
			if(s[m-1] <= e){
				printf("THIS ELEMENT BELONGS TO PSEUDOSPECTRA (%d,%d):%6.10f\n",(iy/NGRID+1),(iy%NGRID+1),s[m-1]);

				plot[iy]=s[m-1];
			 }
				else {

					r = s[m-1]-e;
					
					
					i_max = r/stepx;
					j_max = r/stepy;
					printf("(x=%d,y=%d)\n",i_max,j_max);
					//plot[iy]=0;
					
					/* for x direction */
				//	start_point = iy-j_max*NGRID-i_max;
					//~ 
					//~ for(i=start_point; i<start_point + j_max*2*NGRID + i_max*2 + 1; i++){
						//~ if(i < 0 || (i%NGRID > (iy + i_max )%NGRID) || (i%NGRID < (iy - i_max )%NGRID)) continue;
						//~ else plot[i] = 0;
					//~ }
					
					//~ for(i=start_point; i<start_point + j_max*2*NGRID+1; i=i+NGRID){
							//~ printf("%d\n",i);
						//~ for(j=0; j < 2*i_max+1; j++){
					//~ //		printf("%d\n",i+j);
							//~ if( (i+j) < 0 || (i+j)%NGRID>i_max ) continue;
							//~ 
							//~ else {printf("%d\n",i+j); plot[i+j]=0;}
						//~ }
					//~ }
					int upper_middle = iy-j_max*NGRID;
					int lower_middle = iy + j_max*NGRID;
					double gamma;
					
					for(i = upper_middle; i < lower_middle+1; i=i+NGRID){
					 int indi=j_max;
					// printf("To i einai %d\n",i);
					 for(j=0; j < i_max+1; j++){
						if( (i+j) < 0 || (i+j)/NGRID != i/NGRID ) continue;
						else { /*for now we check every gridpoint of the square if it belongs to the disk*/
							 // printf("%d\n",i+j);
							 gamma = sqrt(pow(abs(indi)*stepy,2) + pow(j*stepx,2));
			   //              printf("The radius in this exclusion disk is %f\n",r);
				//			 printf("to gamma tou %d einai %f\n",i+j,gamma);
							 if(gamma < r) plot[i+j]=0; 
							 }
					    if( (i-j) < 0 || (i-j)/NGRID != i/NGRID ) continue;
					    else { 
					//		printf("The radius in this exclusion disk is %f\n",r);
							gamma = sqrt(pow(abs(indi)*stepy,2) + pow(j*stepx,2));
						//	 printf("to gamma tou %d einai %f\n",i-j,gamma);
							 if(gamma < r) plot[i-j]=0;
							}
			       }
			       indi--;
				}

		//	printf("Iteration:%d\n",iy);


		//	for (i =0; i <NGRID*NGRID; i++){
		//	if(i%NGRID==0) printf("\n");
		//	printf("%f\t",plot[i]);
		//    }
		    
		}
		}
		
		prtdat(NGRID, NGRID, plot, "svd_with_disks.data");
		printf("Total number of svd evaluations in the %d,%d grid is:\t %d\n",NGRID,NGRID,svd_count);
		
		//giving values to data from plot
		for (i = 0; i<NGRID*NGRID; i++)  data[SCALE*(i/NGRID)][SCALE*(i%NGRID)] = plot[i];
	   /////////////////
    BITMAP4 black = {0,0,0,0};
    Draw_Line(image,NGRID,NGRID,x_min,y_min,x_max,y_min,black);
   //////////////////	
		//~ contours[0] = 0.1;
		//~ contours[1] = 0.01;
		//~ contours[2] = 0.001;
		//~ contours[3] = 0.0001;
		//~ contours[4] = 0.00001;
		if ((image = Create_Bitmap(SCALE*NGRID,SCALE*NGRID)) == NULL) {
      fprintf(stderr,"Malloc of bitmap failed\n");
      exit(-1);
   }
 Erase_Bitmap(image,SCALE*NGRID,SCALE*NGRID,grey); /* Not strictly necessary */
   for (j=0;j<SCALE*NGRID;j++) {
      for (i=0;i<SCALE*NGRID;i++) {
         colour = GetColour(data[i][j],0,0.1,1);      /////////////////////////////////////////////
         col.r = colour.r * 255;
        // col.b = colour.b * 255;
       //  Draw_Pixel(image,SCALE*NGRID,SCALE*NGRID,(double)i,(double)j,col);
        //          colour = GetColour(data[i][j],0,0.0001,1);      /////////////////////////////////////////////
       //  col.g = colour.g * 255;
         Draw_Pixel(image,SCALE*NGRID,SCALE*NGRID,(double)i,(double)j,col);
      }
   }

   /* Finally do the contouring */
   CONREC(data,0,SCALE*NGRID-1,0,SCALE*NGRID-1,
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
   Write_Bitmap(fp,image,SCALE*NGRID,SCALE*NGRID,12);
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
  
  
   Draw_Line(image,SCALE*NGRID,SCALE*NGRID,(int)x1,(int)y1,(int)x2,(int)y2,black);

   
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
                  xh[m] =/*x*/ SCALE*creal(zcomplex[i+im[m-1]]);
                  yh[m] = SCALE*cimag(zcomplex[j+jm[m-1]]);
               } else {
                  h[0]  = 0.25 * (h[1]+h[2]+h[3]+h[4]);
                  xh[0] = 0.50 * SCALE*(creal(zcomplex[i]+zcomplex[i+1]));/*x*/
                  yh[0] = 0.50 * SCALE*(cimag(zcomplex[j]+zcomplex[j+1]));
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

