#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <lapacke.h>

#include "../lib/paulslib.h"
#include "../lib/bitmaplib.h"

/* Parameters */
/*  A (m by n) */
#define M		50
#define N		50
#define NGRID	50
#define LDA		N
#define LDU		M
#define LDVT		N
#define XMAX		4	/* maximum boundary of x-axis of the domain */ 
#define XMIN	   -2	/* minimum */
#define YMAX		4	/* maximum boundary of y-axis of the domain */
#define YMIN	   -4	/* minimum */
#define SCALE		1	/* Scale ratio of the data */

/* Constants */
const BITMAP4 
BLACK_COLOR = {0x00, 0x00, 0x00, 0},
WHITE_COLOR = {0xFF, 0xFF, 0xFF, 0},
RED_COLOR   = {0xFF, 0x00, 0x00, 0},
GREEN_COLOR = {0x00, 0xFF, 0x00, 0},
BLUE_COLOR  = {0x00, 0x00, 0xFF, 0},
YELLOW_COLOR= {0xFF, 0xFF, 0x00, 0};

/* Prototypes */
void prtdat(int nx, int ny, double *u1, char *fnam);
void drawline (BITMAP4 *, double, double, double, double, double);
void CONREC(BITMAP4 *image, double **d,int ilb,int iub,int jlb,int jub, double *x,double *y,int nc,double *z,
   void (*ConrecLine)(BITMAP4 *, double,double,double,double,double));

/* DEBUG */
int vectorsdrawn = 0;

int main ( int argc, char **argv ) 
{
  /* Local Variables */
  int i, j, ix, iy, index, ii, jj;

  double x_min, x_max, y_min, y_max, *s, *superb, *plot, **data;
  double stepx, stepy; // Step size for finding gridpoints coordinates in x and y dimension.
  double e=0.1;
  
  double _Complex *z;

  lapack_complex_double *a, *temp, *u, *vt;

  lapack_int 
	info,
	m 	= M, 
	n 	= N,
	lda = LDA,
	ldu = LDU,
	ldvt= LDVT;
  
  /* Memory allocations */
  superb = malloc ( MIN ( m, n ) * sizeof ( double ) );
  temp   = malloc ( ( lda * m ) * sizeof ( lapack_complex_double ) );
  plot   = malloc ( ( NGRID * NGRID ) * sizeof ( double ) );
  vt     = malloc ( ( ldvt * n ) * sizeof ( lapack_complex_double ) );
  a      = malloc ( ( lda * m ) * sizeof ( lapack_complex_double ) );
  u      = malloc ( ( ldu * m ) * sizeof ( lapack_complex_double ) );
  s      = malloc ( m * sizeof ( double ) );
  z      = malloc ( ( NGRID * NGRID ) * sizeof ( double _Complex ) );

  /* Allocating the 2D array data. */
  if ( ( data = malloc ( SCALE * NGRID * sizeof ( double * ) ) ) == NULL ) {
      fprintf ( stderr, "Failed to malloc space for the data\n" );
      exit ( -1 );
  }

  for ( i = 0; i < SCALE * NGRID; i++ ) {
      if ( ( data[i] = malloc ( SCALE * NGRID * sizeof ( double ) ) ) == NULL ) {
	  fprintf ( stderr, "Failed to malloc space for the data\n" );
	  exit ( -1 );
      }
  }

  for ( i = 0; i < SCALE * NGRID; i++ ) {
      for ( j = 0; j < SCALE * NGRID; j++ ) {
	  data[i][j] = 0;
	  //   printf("%f\t",data[i][j]);
      }
  }

  x_min = XMIN;
  x_max = XMAX;
  y_min = YMIN;
  y_max = YMAX;

  /* Initialize grid */
  printf ( "The size of the domain is: X=[%f-%f]  Y=[%f-%f] \n", x_min, x_max, y_min, y_max );

  stepx = ( double ) abs ( x_max - x_min ) / ( NGRID - 1 );
  stepy = ( double ) abs ( y_max - y_min ) / ( NGRID - 1 );

  printf ( "To stepx einai %f\n", stepx );
  printf ( "To stepy einai %f\n", stepy );

  for ( i = 0; i < NGRID * NGRID; i++ ) {
      z[i] = x_min + ( i / n * stepx ) + ( y_min + ( i % n * stepy ) ) * I;
  }

  memset ( temp, 0, ( lda * m ) * sizeof ( *temp ) );
  memset ( a, 0, ( lda * m ) * sizeof ( *a ) );
  memset ( u, 0, ( ldu * m ) * sizeof ( *u ) );
  memset ( vt, 0, ( ldvt * m ) * sizeof ( *vt ) );

  j = 0;
  for ( i = 0; i < lda * m; i = i + n ) {
	if ( i == 0 ) {
	  a[i] = lapack_make_complex_double ( 1, 0 );
	  a[i + 1] = lapack_make_complex_double ( 1, 0 );
	  a[i + 2] = lapack_make_complex_double ( 1, 0 );
	  a[i + 3] = lapack_make_complex_double ( 1, 0 );
	} else if ( i == ( n - 3 ) * n ) {
	  a[i + j] = lapack_make_complex_double ( -1, 0 );
	  a[i + ( j + 1 )] = lapack_make_complex_double ( 1, 0 );
	  a[i + ( j + 2 )] = lapack_make_complex_double ( 1, 0 );
	  a[i + ( j + 3 )] = lapack_make_complex_double ( 1, 0 );
	  j++;
	} else if ( i == ( n - 2 ) * n ) {
	  a[i + j] = lapack_make_complex_double ( -1, 0 );
	  a[i + ( j + 1 )] = lapack_make_complex_double ( 1, 0 );
	  a[i + ( j + 2 )] = lapack_make_complex_double ( 1, 0 );
	  j++;
	} else if ( i == ( n - 1 ) * n ) {
	  a[i + j] = lapack_make_complex_double ( -1, 0 );
	  a[i + ( j + 1 )] = lapack_make_complex_double ( 1, 0 );
	  j++;
	} else {
	  a[i + j] = lapack_make_complex_double ( -1, 0 );
	  a[i + ( j + 1 )] = lapack_make_complex_double ( 1, 0 );
	  a[i + ( j + 2 )] = lapack_make_complex_double ( 1, 0 );
	  a[i + ( j + 3 )] = lapack_make_complex_double ( 1, 0 );
	  a[i + ( j + 4 )] = lapack_make_complex_double ( 1, 0 );
	  j++;
	}
  }

  //print_matrix("Entry Matrix A", m, n, a, lda );
  for ( iy = 0; iy < NGRID * NGRID; iy++ )
  {
	memcpy ( temp, a, ( lda * m ) * sizeof ( *temp ) );

	for ( i = 0; i < lda * m; i = i + ( n + 1 ) ){
	  temp[i] = a[i] - z[iy];
	}

	printf("LAPACKE_zgesvd (row-major, high-level) Example Program Results(%d,%d)\n", iy / NGRID, iy % NGRID );

	/* Compute SVD */
	info = LAPACKE_zgesvd ( LAPACK_ROW_MAJOR, 'N', 'N', m, n, temp, lda, s, NULL, ldu, NULL, ldvt, superb );

	/* Check for convergence */
	if ( info > 0 ) {
	  printf ( "The algorithm computing SVD failed to converge.\n" );
	  exit (1);
	}

	if ( s[m-1] <= e ) {
	  printf( "THIS ELEMENT BELONGS TO PSEUDOSPECTRA (%d,%d):%6.10f\n", ( iy / NGRID + 1 ), ( iy % NGRID + 1 ), s[m - 1] );
	  plot[iy] = s[m - 1];
	}
	else
	  plot[iy] = 0;

  }
  
  prtdat ( NGRID, NGRID, plot, "svd_no_disks.data" );

  /* GENERATE CONTOUR */
  #define NCONTOUR 	3	/* Number of contour levels */
  double contours[NCONTOUR]; 	/* Contains the desiRED_COLOR contour levels */
  double *xaxis,*yaxis;		/* Arrays for the x axis and y axis coordinates */
  BITMAP4 *image; 				/* Image on which the contours will be drawn, see bitmaplib.c */
  COLOUR colour;
  BITMAP4 col, bgcolor = BLACK_COLOR;

  /* Giving values to data from plot */
  for ( i = 0; i < NGRID * NGRID; i++ )
      data[SCALE * ( i / NGRID )][SCALE * ( i % NGRID )] = plot[i];

  /*  Set up the axis coordinates */
   if ((xaxis = malloc(SCALE*NGRID*sizeof(double))) == NULL) {
      fprintf(stderr,"Failed to malloc space for the xaxis data\n");
      exit(-1);
   }
   for (i=0;i<SCALE*NGRID;i++) 
      xaxis[i] = i;
   if ((yaxis = malloc(SCALE*NGRID*sizeof(double))) == NULL) {
      fprintf(stderr,"Failed to malloc space for the yaxis data\n");
      exit(-1);
   }
   for (i=0;i<SCALE*NGRID;i++) 
      yaxis[i] = i;

  /* Set up the contour levels */
  contours[0] = 0.1 / 50000.0;
  //contours[1] = 0.1 / 3000.0;
  contours[1] = 0.1 / 80.0;
  //contours[3] = 0.1 / 5.0;
  contours[2] = 0.1;

  if ( (image = Create_Bitmap ( SCALE * NGRID, SCALE * NGRID ) ) == NULL ) {
      fprintf ( stderr, "Malloc of bitmap failed\n" );
      exit ( -1 );
  }

  Erase_Bitmap ( image, SCALE * NGRID, SCALE * NGRID, bgcolor );	/* Not strictly necessary */

  for ( j = 0; j < SCALE * NGRID; j++ ) {
	for ( i = 0; i < SCALE * NGRID; i++ ) {
	  if ( data[i][j] > 0 ) {
		if ( data[i][j] < contours[1] )
		  col = YELLOW_COLOR;
		else if ( data[i][j] < contours[2] )
		  col = GREEN_COLOR;
		else if ( data[i][j] < contours[3] )
		  col = WHITE_COLOR;
		else
		  col = RED_COLOR;
	  } else 
	  {
		col = bgcolor;
	  }
	  Draw_Pixel ( image, SCALE * NGRID, SCALE * NGRID, ( double ) i, ( double ) j, col );
	}
  }

  /* Finally do the contouring */
  //CONREC(image, data, 0, SCALE*NGRID-1, 0, SCALE*NGRID-1, z, NCONTOUR, contours, drawline);
  
  CONREC(image, data,0,SCALE*NGRID-1,0,SCALE*NGRID-1,
		 xaxis,yaxis,NCONTOUR,contours,drawline);
  
  fprintf(stderr,"Drew %d vectors\n",vectorsdrawn);
  
  FILE *fp = fopen ( "image_NO_SCALE.tga", "w" );
  Write_Bitmap ( fp, image, SCALE*NGRID, SCALE*NGRID, 12 );

  /* Create a bigger image */
  if (0) {
	const int sc = 10;
	BITMAP4 *bigImage = bigImage = Create_Bitmap ( sc*SCALE*NGRID, sc*SCALE*NGRID );

	fp = fopen ( "image_BICUBIC_SCALE.tga", "w" );
	BiCubicScale ( image, SCALE*NGRID, SCALE*NGRID, bigImage, sc*SCALE*NGRID, sc*SCALE*NGRID );
	Write_Bitmap ( fp, bigImage, sc*SCALE*NGRID, sc*SCALE*NGRID, 12 );

	fp = fopen ( "image_GAUSSIAN_SCALE.tga", "w" );
	GaussianScale ( image, SCALE*NGRID, SCALE*NGRID, bigImage, sc*SCALE*NGRID, sc*SCALE*NGRID, 0.5 );
	Write_Bitmap ( fp, bigImage, sc*SCALE*NGRID, sc*SCALE*NGRID, 12 );
  }
  
  fclose ( fp );
  exit ( 0 );

} /* End of LAPACKE_zgesvd Example */

/* Functions */
void print_matrix (char *desc, lapack_int m, lapack_int n, lapack_complex_double * a, lapack_int lda)
{
    lapack_int i, j;
    printf("\n %s\n", desc);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
            printf(" (%6.2f,%6.2f)",
                   lapack_complex_double_real(a[i * lda + j]),
                   lapack_complex_double_imag(a[i * lda + j]));
        printf("\n");
    }
}

void print_int_vector (char *desc, lapack_int n, lapack_int * a)
{
    lapack_int j;
    printf("\n %s\n", desc);
    for (j = 0; j < n; j++)
        printf(" %6i", a[j]);
    printf("\n");
}

void prtdat (int nx, int ny, double *u1, char *fnam)
{
    int ix;
    FILE *fp = fopen(fnam, "w");

    for (ix = 0; ix < nx * ny; ix++) {
        fprintf(fp, "%f", u1[ix]);
        if (ix % nx == nx - 1)
            fprintf(fp, "\n");
        else
            fprintf(fp, " ");
    }

    fclose(fp);
}

void drawline (BITMAP4 *image, double x1, double y1, double x2, double y2, double z)
{
    Draw_Line(image, SCALE * NGRID, SCALE * NGRID, (int)x1, (int)y1, (int)x2, (int)y2, BLUE_COLOR);
    vectorsdrawn++;
}

/**
 *  Derivation from CONREC
 *  d               ! matrix of data to contour
 *  ilb,iub,jlb,jub ! index bounds of data matrix
 *  x               ! data matrix column coordinates
 *  y               ! data matrix row coordinates
 *  nc              ! number of contour levels
 *  z               ! contour levels in increasing order
 */
void CONREC(BITMAP4 *image, double **d,int ilb,int iub,int jlb,int jub,
   double *x,double *y,int nc,double *z,
   void (*ConrecLine)(BITMAP4 *, double,double,double,double,double))
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
                  xh[m] = x[i+im[m-1]];
                  yh[m] = y[j+jm[m-1]];
               } else {
                  h[0]  = 0.25 * (h[1]+h[2]+h[3]+h[4]);
                  xh[0] = 0.50 * (x[i]+x[i+1]);
                  yh[0] = 0.50 * (y[j]+y[j+1]);
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
               ConrecLine(image, x1,y1,x2,y2,z[k]);
            } /* m */
         } /* k - contour */
      } /* i */
   } /* j */
}

