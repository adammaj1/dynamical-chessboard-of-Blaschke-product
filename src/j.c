/*
https://arxiv.org/pdf/1510.00043v2.pdf

Parabolic and near-parabolic renormalizations for local degree three by Fei Yang

B3 is a Blaschke product whose Julia set is the unit circle. The point z = 1
is a 1-parabolic fixed point with two attracting petals. In particular, the unit disk D
and C rb D are two immediate basins of 1
the first column of figure 14 indicate the dynamical chessboards of f and B3.



B3
(%i9) display2d:false;

(%o9) false
(%i10) a;

(%o10) ((z+1/2)/(z/2+1))^3
(%i11) ratsimp(a);

(%o11) (8*z^3+12*z^2+6*z+1)/(z^3+6*z^2+12*z+8)
(%i12) ratexpand(a);

(%o12) (8*z^3)/(z^3+6*z^2+12*z+8)+(12*z^2)/(z^3+6*z^2+12*z+8)
                                 +(6*z)/(z^3+6*z^2+12*z+8)
                                 +1/(z^3+6*z^2+12*z+8)





t B3 is a Blaschke product whose Julia set is the unit circle. The point z = 1
is a 1-parabolic fixed point with two attracting petals. 
In particular, the unit disk D and C rb D are two immediate basins of 1.

the dynamical chessboards of f and B3.

  here are:
  * 1 critical point  z=0.0
  * 1 a nondegenerated 1-parabolic fixed point
  

 

  Adam Majewski
  adammaj1 aaattt o2 dot pl  // o like oxygen not 0 like zero 
  
  
  
  Structure of a program or how to analyze the program 
  
  
  ============== Image X ========================
  
  DrawImageOf -> DrawPointOf -> ComputeColorOf ( FunctionTypeT FunctionType , complex double z) -> ComputeColor
  
  
  check only last function  which computes color of one pixel for given Function Type
  
  

   
  ==========================================

  
  ---------------------------------
  indent d.c 
  default is gnu style 
  -------------------



  c console progam 
  
  export  OMP_DISPLAY_ENV="TRUE"	
  gcc d.c -lm -Wall -march=native -fopenmp
  time ./a.out > b.txt


  gcc d.c -lm -Wall -march=native -fopenmp


  time ./a.out

  time ./a.out >i.txt
  time ./a.out >e.txt
  
  
  
  
  
  
  convert -limit memory 1000mb -limit disk 1gb dd30010000_20_3_0.90.pgm -resize 2000x2000 10.png

  
  
  
*/

#include <stdio.h>
#include <stdlib.h>		// malloc
#include <string.h>		// strcat
#include <math.h>		// M_PI; needs -lm also
#include <complex.h>
#include <omp.h>		// OpenMP
#include <limits.h>		// Maximum value for an unsigned long long int



// https://sourceforge.net/p/predef/wiki/Standards/

#if defined(__STDC__)
#define PREDEF_STANDARD_C_1989
#if defined(__STDC_VERSION__)
#if (__STDC_VERSION__ >= 199409L)
#define PREDEF_STANDARD_C_1994
#endif
#if (__STDC_VERSION__ >= 199901L)
#define PREDEF_STANDARD_C_1999
#endif
#endif
#endif




/* --------------------------------- global variables and consts ------------------------------------------------------------ */

// 


/* 
FunctionType = representing functions
BD = Binary decomposition
MBD  = Modified BD is better, so BD is not used


*/


typedef enum  {FatouBasins = 0, FatouComponents = 2,  LSM = 3, LS2M = 4, Unknown = 5 , BD = 6, MBD = 7 , SAC = 8, DLD = 9, ND = 10 , NP= 11, POT = 12 , Blend = 13, DEM = 14, IBD = 15, ParabolicCheckerboard = 16, ParabolicCheckerboard2 = 17
		
} FunctionTypeT; 
// FunctionTypeT FunctionType;

// virtual 2D array and integer ( screen) coordinate
// Indexes of array starts from 0 not 1 
//unsigned int ix, iy; // var
static unsigned int ixMin = 0;	// Indexes of array starts from 0 not 1
static unsigned int ixMax;	//
static unsigned int iWidth;	// horizontal dimension of array

static unsigned int iyMin = 0;	// Indexes of array starts from 0 not 1
static unsigned int iyMax;	//

static unsigned int iHeight = 2000;	//  
// The size of array has to be a positive constant integer 
static unsigned long long int iSize;	// = iWidth*iHeight; 

// memmory 1D array 
unsigned char *data;
unsigned char *edge;
unsigned char *edge2;
//unsigned char *edge2;

// unsigned int i; // var = index of 1D array
//static unsigned int iMin = 0; // Indexes of array starts from 0 not 1
unsigned int iMax;	// = i2Dsize-1  = 
// The size of array has to be a positive constant integer 
// unsigned int i1Dsize ; // = i2Dsize  = (iMax -iMin + 1) =  ;  1D array with the same size as 2D array



// see SetPlane

double radius = 1.5; 
complex double plane_center = 0.0 ;
double  DisplayAspectRatio  = 1.0; // https://en.wikipedia.org/wiki/Aspect_ratio_(image)
// dx = dy compare setup : iWidth = iHeight;
double ZxMin; //= -1.3;	//-0.05;
double ZxMax;// = 1.3;	//0.75;
double ZyMin;// = -1.3;	//-0.1;
double ZyMax;// = 1.3;	//0.7;
double PixelWidth;	// =(ZxMax-ZxMin)/ixMax;
double PixelHeight;	// =(ZyMax-ZyMin)/iyMax;



double ratio; 

complex double trap_center;
complex double trap_center2;
double ER;
double ER2;			//= 1e60;
double AR; // bigger values do not works
double AR2;




int IterMax = 1000;
int IterMax_LSM = 1000;
//int IterMax_DEM = 10000000;

/* colors = shades of gray from 0 to 255 

   unsigned char colorArray[2][2]={{255,231},    {123,99}};
   color = 245;  exterior 
   
   here are two period 2 basins: basin1  and basin2
   each basin is a basin of attraction of period 2 cycle
   
   Each cycle has immediate basin of attraction which consist of 2 components ( and it's preimages)
   so we need 4 colors 
   
   also exterior is a component oof one basin , 
   it is not a basin of attraction to infiiniity
   
   
   
   
   
*/
unsigned char iColorOfBasin1 = 170;
unsigned char iColorOfInterior = 150;
unsigned char iColorOfInterior1 = 250;
unsigned char iColorOfExterior = 225;

// for parabolic chessboards
unsigned char colorArray[2][2]={{255,231},
				   {123,99}}; /* shades of gray used in image */

unsigned char iColorOfBoundary = 0;
unsigned char iColorOfUnknown = 5;





// pixel counters
unsigned long long int uUnknown = 0;
unsigned long long int uInterior = 0;
unsigned long long int uExterior = 0;



/* critical point */

const int period = 1;
complex double zcr = -0.5; //

complex double c ;

// 

complex double zp0;  // a nondegenerated 1-parabolic fixed point


// for MBD
static double TwoPi=2.0*M_PI; // texture
double t0 ; // manually tuned t for MBD
// see https://www.youtube.com/watch?v=JttLtB0Gkdk&t=894s
// 


/* ------------------------------------------ functions -------------------------------------------------------------*/

/* 

*/


// complex function b3 
complex double f(const complex double z0) {

  double complex z = z0;
  // ((z+1/2)/(z/2+1))^3
  complex double n = z +0.5;
  complex double d = 1.0 +(z/2.0); 
  z = cpow(n/d, 3.0);
  return  z;
}
	


double c_arg(complex double z)
{
 double arg;
 arg = carg(z);
 if (arg<0.0) arg+= TwoPi ; 
 return arg; 
}

double c_turn(complex double z)
{
 double arg;
 arg = c_arg(z);
 return arg/TwoPi; 
}





int is_z_outside(complex double z){

  if (creal(z) >ZxMax ||
      creal(z) <ZxMin ||
      cimag(z) >ZyMax ||
      cimag(z) <ZyMin)
    {return 1; } // is outside = true
      
    
  return 0; // is inside = false



}









// from screen to world coordinate ; linear mapping
// uses global cons
double GiveZx (int ix)
{
  return (ZxMin + ix * PixelWidth);
}

// uses globaal cons
double GiveZy (int iy)
{
  return (ZyMax - iy * PixelHeight);
}				// reverse y axis


complex double GiveZ (int ix, int iy)
{
  double Zx = GiveZx (ix);
  double Zy = GiveZy (iy);

  return Zx + Zy * I;




}







//------------------complex numbers -----------------------------------------------------

double cabs2(complex double z){

  return creal(z)*creal(z)+cimag(z)*cimag(z);


}






/* -----------  array functions = drawing -------------- */

/* gives position of 2D point (ix,iy) in 1D array  ; uses also global variable iWidth */
unsigned int Give_i (unsigned int ix, unsigned int iy)
{
  return ix + iy * iWidth;
}












/* 

   is it possible to adjust AR so that level curves in interior have figure 8?

   find such AR for internal LCM/J and LSM that level curves croses critical point and it's preimages
   for attracting ( also weakly attracting = parabolic) dynamics

   it may fail 
   * if one iteration is bigger then smallest distance between periodic point zp0 and Julia set
   * if critical point is attracted by another cycye ( then change periodic point zp0)

   Made with help of Claude Heiland-Allen


   attracting radius of circle around finite attractor
   there are 2 basins so  
  
  
   It would have to be done separately in each basin.

   A suggested method:

   For each critical point, forward iterate to find an attractor and then thin out the critical point set to only one per basin by removing all but one that converge to a common attractor, for each attractor.
   For each pixel, calculate a smoothed iteration value (e.g. using the methods in my GVC coloring ucl) and note which basin it is in.
   For each critical point in the reduced set, calculate a smoothed iteration value using the same method as in step 2.
   For each pixel, subtract from its smoothed iteration value the one found in step 3 for the critical point that shares its basin. Note that the critical point itself, if inside the image rectangle and in a pixel center, will end up with zero and some points may end up with negative values.
   The level set boundaries you want will now be the boundaries where the sign or the integer part of the modified smoothed iteration value changes. In particular, the -0.something to +0.something transition will pass through the critical point, the n.something to (n+1).something transitions for nonnegative n will pass through its images, and the same for negative n will pass through its preimages.

   pauldebrot 
   https://fractalforums.org/programming/11/crtical-points-and-level-curves/4323/msg29514#new
  


*/
double GiveTunedAR(const double iter_Max){

  fprintf(stdout, " Tuned AR = \n");

  complex double z = zcr; // initial point z0 = criical point 
  double iter;
  double r ;//= 10 * PixelWidth; // initial value 
  //double rMin = 30 * PixelWidth;
  // double t;
  
  // iterate critical point
  for (iter=0; iter< iter_Max; iter+=1.0 ){
	  
    
    //if ( r<rMin) {break;}
		
    z = f(z); // forward iteration
  	
  }
  // check distance between zn = f^n(zcr) and periodic point zp0
  
 
  
  r = cabs(z - zp0)/2.0;
  // use it as a AR
  return r;
	
	
}








// ****************** DYNAMICS = trap tests ( target sets) ****************************


// ???????

int IsInsideTrap(int ix, int iy){


  complex double z = GiveZ(ix, iy);
  if (  cabs2(trap_center -z) < AR2 )
    {return 1;}
  return 0;



}







/*
  1 basin  = not works here, because whole plane / sphere/ rectanlge is the same , the only one basin
  - unknown ( possibly empty set ) 

*/

unsigned char ComputeColorOfFatouBasins (complex double z)
{



	
	
  int i;			// number of iteration
  for (i = 0; i < IterMax; ++i)
    {


		
      // infinity is superattracting here !!!!!	
       if ( cabs2(trap_center2-z) < AR2 ){ return iColorOfInterior1;}
      // 1 Attraction basins 
      if ( cabs2(trap_center-z) < AR2 ){ return iColorOfInterior;}
		 
      			
     
      z = f(z);		//  iteration: z(n+1) = f(zn)
	

    }

  
  return iColorOfUnknown;


}


/*
  2 basins 
  
  - - basin 1 
  - - basin 2
  - unknown ( possibly empty set ) 

*/

unsigned char ComputeColorOfFatouComponents (complex double z)
{



	
	
  


  int i;			// number of iteration
  for (i = 0; i < IterMax; ++i)
    {


      // infinity is superattracting here !!!!!	
       if ( cabs2(z) > ER2 ){ return iColorOfExterior;}
	
  
      //1 Attraction basins 
      if ( cabs2(trap_center-z) < AR2 ){ return iColorOfBasin1 - (i % period)*20;}
	 
		
     
      z = f(z);		//  iteration: z(n+1) = f(zn)
	

    }

  
  return iColorOfUnknown;


}




/*
attracting petals ( gray curves)
take 2 points: last point of critical orbit and fixed point.
draw circle which is passing thru above 2 points and with diameter equal to distance between such 2 points. Such circle is the smallest ( here not in general) attracting petal

*/

unsigned char ComputeColorOfLSM (complex double z)
{



	
	
  //double cabs2z;
  //double cabs2zAR;


  int i;			// number of iteration
  for (i = 0; i < IterMax_LSM; ++i)
    {


	//cabs2z = cabs2(z);	
	//cabs2zAR = cabs2(trap_center - z);

      // infinity is superattracting here !!!!!	
       
       if ( cabs2(trap_center2-z) < AR2 || ( cabs2(trap_center-z) < AR2 ))
		 
       		{ return (15*i) % 255;} // cabs2(zp0-z) = cabs2(z) because zp0 = zcr = 0
	 
     
	
      z = f(z);	

    }

  return iColorOfUnknown;


}

/*
attracting petals ( gray curves)
take 2 points: last point of critical orbit and fixed point.
draw circle which is passing thru above 2 points and with diameter equal to distance between such 2 points. Such circle is the smallest ( here not in general) attracting petal

*/

unsigned char ComputeColorOfLS2M (complex double z)
{



	
	
  double cabs2z;
  double cabs2zAR;


  int i;			// number of iteration
  for (i = 0; i < IterMax_LSM; ++i)
    {


	//cabs2z = cabs2(z);	
	//cabs2zAR = cabs2(trap_center - z);

     
       if (  cabs2(trap_center2-z) < AR2 || ( cabs2(trap_center-z) < AR2 ))
       		{ 
       			if (i %2) 
       				{return 255;}
       				else {return 230;}}
       		//return (15*i) % 255;} // cabs2(zp0-z) = cabs2(z) because zp0 = zcr = 0
	 
     
	
      z = f(z);	

    }

  return iColorOfUnknown;


}


unsigned char ComputeColorOfBD (complex double z)
{



	
	
  double cabs2z;
  double cabs2zAR;

  int i;			// number of iteration
  for (i = 0; i < IterMax_LSM; ++i)
    {


	//cabs2z = cabs2(z); // numerical speed up : cabs2(zp0-z) = cabs2(z) because zp0 = zcr = 0	
	//cabs2zAR = cabs2(trap_center - z);
      // 
       if (  cabs2(trap_center2-z) < AR2 || ( cabs2(trap_center-z) < AR2 ) ) // if z is inside target set ( orbit trap) 
       		{ 
       			if (creal(z) > 0) // binary decomposition of target set
       				{  return 0;}
       				else {return 255; }
     
      
      		}
	 
     
	
      z = f(z);	

    }

  return iColorOfUnknown;


}


// Modified BD
unsigned char ComputeColorOfMBD (complex double z)
{



	
	
  double cabs2z;
  double cabs2zAR;
  double turn; 

  int i;			// number of iteration
  for (i = 0; i < IterMax_LSM; ++i)
    {


	cabs2z = cabs2(z); // numerical speed up : cabs2(zp0-z) = cabs2(z) because zp0 = zcr = 0
	cabs2zAR = cabs2(trap_center - z);	

      //  if z is inside target set ( orbit trap) = exterior of circle with radius ER 
       if ( cabs2z > ER2  ) // exterior
       		{ 
       			if (creal(z) > 0) // binary decomposition of target set
       				{  return 150;}
       				else {return 255; }
     
      
      		}
      		
      	if ( cabs2zAR  < AR2 ) // if z is inside target set ( orbit trap) = interior of cirlce with radius AR
      		{
      			turn = c_turn(z);
      			if (turn < t0 || turn > t0+0.5) // modified binary decomposition of target set
      				{  return 150;}
       				else {return 255; }
     
      		
      		}
	 
     
	
      z = f(z);	

    }

  return iColorOfUnknown;


}


// Modified BD
unsigned char ComputeColorOfIBD (complex double z)
{



	
	
  double cabs2z;
  double cabs2zAR;
  double turn; 

  int i;			// number of iteration
  for (i = 0; i < IterMax_LSM; ++i)
    {


	cabs2z = cabs2(z); // numerical speed up : cabs2(zp0-z) = cabs2(z) because zp0 = zcr = 0
	cabs2zAR = cabs2(trap_center - z);	

      //  if z is inside target set ( orbit trap) = exterior of circle with radius ER 
       if ( cabs2z > ER2  ) // exterior
       		{ 
       			return iColorOfExterior;
     
      
      		}
      		
      	if ( cabs2zAR  < AR2 ) // if z is inside target set ( orbit trap) = interior of cirlce with radius AR
      		{
      			turn = c_turn(z);
      			if (turn < t0 || turn > t0+0.5) // modified binary decomposition of target set
      				{  return 150;}
       				else {return 255; }
     
      		
      		}
	 
     
	
      z = f(z);	

    }

  return iColorOfUnknown;


}




 



// 
unsigned char ComputeColorOfParabolicCheckerboard (complex double z)
{



	
	
  double cabs2z;
  double cabs2zAR;
  //double turn; 
  int m;
  int n;

  int i;			// number of iteration
  for (i = 0; i < IterMax_LSM; ++i)
    {


	cabs2z = cabs2(z); // numerical speed up : cabs2(zp0-z) = cabs2(z) because zp0 = zcr = 0
	cabs2zAR = cabs2(trap_center - z);	

      //  if z is inside target set ( orbit trap) = exterior of circle with radius ER 
       if ( cabs2z > ER2  ) // exterior
       		{ 
       			return iColorOfExterior;
     
      
      		}
      		
      	if ( cabs2zAR  < AR2 ) // if z is inside target set ( orbit trap) = interior of cirlce with radius AR
      		{
      			m = (cimag(z) > 0 ? 0 : 1); // petal part
   			n = (i % 2); // attraction time 
   			return colorArray[m][n]; //iColor
     
      		
      		}
	 
     
	
      z = f(z);	

    }

  return iColorOfUnknown;


}



// 
unsigned char ComputeColorOfParabolicCheckerboard2 (complex double z)
{



	
	
  double cabs2z;
  double cabs2zAR;
  double angle; 
  int pMax = 6; // ? child period  ?

  int i;			// number of iteration
  for (i = 0; i < IterMax_LSM; ++i)
    {


	cabs2z = cabs2(z); // numerical speed up : cabs2(zp0-z) = cabs2(z) because zp0 = zcr = 0
	cabs2zAR = cabs2(trap_center - z);	

      //  if z is inside target set ( orbit trap) = exterior of circle with radius ER 
       if ( cabs2z > ER2  ) // exterior
       		{ 
       			return iColorOfExterior;
     
      
      		}
      		
      		
      	for (int p=0; p < pMax; p++){ 	
      	if ( cabs2zAR  < AR2 ) // if z is inside target set ( orbit trap) = interior of cirlce with radius AR
      		{
      			angle = c_turn(z - trap_center); // now in (0,1) range
      		  
      		  
      		  angle = angle*200.0; // repeated gradient
      		  
      		  return angle* 255; // now in (0,255) range
      		       
      		
      		}
      		
	 
     
	
      z = f(z);	
      }

    }

  return iColorOfUnknown;


}



/* ==================================================================================================
   ============================= Draw functions ===============================================================
   =====================================================================================================
*/ 
unsigned char ComputeColor(FunctionTypeT FunctionType, complex double z){

  unsigned char iColor;
	
	
	
  switch(FunctionType){
  
  case FatouBasins :{iColor = ComputeColorOfFatouBasins(z); break;}
  	
  case FatouComponents :{iColor = ComputeColorOfFatouComponents(z); break;}
  
 
  
  case LSM :{iColor = ComputeColorOfLSM(z); break;}
  
  case LS2M : {iColor = ComputeColorOfLS2M(z); break;  }	
  
  
  // case DEM : {iColor = ComputeColorOfDEMJ(z); break;}
	
   	
  	//case Unknown : {iColor = ComputeColorOfUnknown(z); break;}
		
  	case BD : {iColor = ComputeColorOfBD(z); break;}
		
  	case MBD : {iColor = ComputeColorOfMBD(z); break;}
  	
  	case IBD : {iColor = ComputeColorOfIBD(z); break;}
  	
  	case ParabolicCheckerboard: {iColor = ComputeColorOfParabolicCheckerboard (z); break;}
  	
  	case ParabolicCheckerboard2: {iColor = ComputeColorOfParabolicCheckerboard2 (z); break;}
		
  	//case SAC : {iColor = ComputeColorOfSAC(z); break;}
  
  	//case DLD : {iColor = ComputeColorOfDLD(z); break;}
		
  	//case ND : {iColor = ComputeColorOfND(z); break;}
		
  	//case NP : {iColor = ComputeColorOfNP(z); break;}
		
  	//case POT : {iColor = ComputeColorOfPOT(z); break;}
		
  	//case Blend : {iColor = ComputeColorOfBlend(z); break;}
   	
  	
  
  	
  	
	
  default: {}
	
	
  }
	
  return iColor;



}


// plots raster point (ix,iy) 
int DrawPoint ( unsigned char A[], FunctionTypeT FunctionType, int ix, int iy)
{
  int i;			/* index of 1D array */
  unsigned char iColor;
  complex double z;


  i = Give_i (ix, iy);		/* compute index of 1D array from indices of 2D array */
  if(i<0 && i> iMax)
    { return 1;}
  
  z = GiveZ(ix,iy);
  iColor = ComputeColor(FunctionType, z);
  A[i] = iColor ;		// 
  		
  	  
  return 0;
}




int DrawImage ( unsigned char A[], FunctionTypeT FunctionType)
{
  unsigned int ix, iy;		// pixel coordinate 

  fprintf (stderr, "compute image %d \n", FunctionType);
  // for all pixels of image 
#pragma omp parallel for schedule(dynamic) private(ix,iy) shared(A, ixMax , iyMax, uUnknown, uInterior, uExterior)
  for (iy = iyMin; iy <= iyMax; ++iy)
    {
      fprintf (stderr, " %d from %d \r", iy, iyMax);	//info 
      for (ix = ixMin; ix <= ixMax; ++ix)
	DrawPoint(A, FunctionType, ix, iy);	//  
    }
  fprintf (stderr, "\n");	//info 
  return 0;
}







int PlotPoint(const complex double z, unsigned char A[]){

	
  unsigned int ix = (creal(z)-ZxMin)/PixelWidth;
  unsigned int iy = (ZyMax - cimag(z))/PixelHeight;
  unsigned int i = Give_i(ix,iy); /* index of _data array */
	
  	
	
  if(i>-1 && i< iMax)
    {A[i]= 0; // 255-A[i];
    }
	
	
  return 0;
	
}






int IsInsideCircle (int x, int y, int xcenter, int ycenter, int r){

	
  double dx = x- xcenter;
  double dy = y - ycenter;
  double d = sqrt(dx*dx+dy*dy);
  if (d<r) {    return 1;}
  return 0;
	  

} 

// Big point = disk 
int PlotBigPoint(const complex double z, unsigned char iColor,  double p_size, unsigned char A[]){

	
  unsigned int ix_seed = (creal(z)-ZxMin)/PixelWidth;
  unsigned int iy_seed = (ZyMax - cimag(z))/PixelHeight;
  unsigned int i;
	
	
  if (  is_z_outside(z)) 
    {fprintf (stdout,"PlotBigPoint :  z= %.16f %+.16f*I is outside\n", creal(z), cimag(z)); return 1;} // do not plot	
	
  /* mark seed point by big pixel */
  int iSide =p_size*iWidth/2000.0 ; /* half of width or height of big pixel */
  int iY;
  int iX;
  for(iY=iy_seed-iSide;iY<=iy_seed+iSide;++iY){ 
    for(iX=ix_seed-iSide;iX<=ix_seed+iSide;++iX){ 
      if (IsInsideCircle(iX, iY, ix_seed, iy_seed, iSide)) {
	i= Give_i(iX,iY); /* index of _data array */
	//if(i>-1 && i< iMax)
	A[i]= iColor; // 255; ;//  -A[i];
      }
      // else {printf(" bad point \n");}
	
    }}
	
	
  return 0;
	
}



int PlotAllPoints(const complex double zz[], int kMax, unsigned char iColor, double p_size,unsigned char A[]){

  int k;
	
	
  //printf("kMax = %d \n",kMax);
	

  for (k = 0; k < kMax; ++k)
    {
      //fprintf(stderr, "z = %+f %+f \n", creal(zz[k]),cimag(zz[k]));
      PlotBigPoint(zz[k], iColor, p_size, A);}
  return 0;





}




int DrawForwardOrbit(const complex double z0, const unsigned long long int i_Max,unsigned char iColor, double p_size, unsigned char A[]){
 

  
  unsigned long long int i; /* nr of point of critical orbit */
  complex double z = z0;
  printf("draw forward orbit \n");
 
  PlotBigPoint(z, iColor, p_size, A);
  
  /* forward orbit of critical point  */
  for (i=1;i<i_Max ; ++i)
    {
      z  = f(z);
      //if (cabs2(z - z2a) > 2.0) {return 1;} // escaping
      PlotBigPoint(z, iColor, p_size/2 , A);
    }
  
  fprintf (stdout,"first point of the orbit z0= %.16f %+.16f*I \n", creal(z0), cimag(z0));
  fprintf (stdout,"last point of the orbit z= %.16f %+.16f*I \n", creal(z), cimag(z));
   
  return 0;
 
}



// ***********************************************************************************************
// ********************** draw line segment ***************************************
// ***************************************************************************************************




// plots raster point (ix,iy) 
int iDrawPoint(unsigned int ix, unsigned int iy, unsigned char iColor, unsigned char A[])
{ 

  /* i =  Give_i(ix,iy) compute index of 1D array from indices of 2D array */
  if (ix >=ixMin && ix<=ixMax && iy >=iyMin && iy<=iyMax )
    {A[Give_i(ix,iy)] = iColor;}
  else {fprintf (stdout,"iDrawPoint :   (%d; %d) is outside\n", ix,iy); }

  return 0;
}



/*
  http://rosettacode.org/wiki/Bitmap/Bresenham%27s_line_algorithm
  Instead of swaps in the initialisation use error calculation for both directions x and y simultaneously:
*/
void iDrawLine( int x0, int y0, int x1, int y1, unsigned char iColor, unsigned char A[]) 
{
  int x=x0; int y=y0;
  int dx = abs(x1-x0), sx = x0<x1 ? 1 : -1;
  int dy = abs(y1-y0), sy = y0<y1 ? 1 : -1; 
  int err = (dx>dy ? dx : -dy)/2, e2;

  for(;;){
    iDrawPoint(x, y, iColor, A);
    if (x==x1 && y==y1) break;
    e2 = err;
    if (e2 >-dx) { err -= dy; x += sx; }
    if (e2 < dy) { err += dx; y += sy; }
  }
}




int dDrawLineSegment(double complex Z0, double complex Z1, int color, unsigned char *array) 
{

  double Zx0 = creal(Z0);
  double Zy0 = cimag(Z0);
  double Zx1 = creal(Z1);
  double Zy1 = cimag(Z1);
  unsigned int ix0, iy0; // screen coordinate = indices of virtual 2D array 
  unsigned int ix1, iy1; // screen coordinate = indices of virtual 2D array

  // first step of clipping
  //if (  Zx0 < ZxMax &&  Zx0 > ZxMin && Zy0 > ZyMin && Zy0 <ZyMax 
  // && Zx1 < ZxMax &&  Zx1 > ZxMin && Zy1 > ZyMin && Zy1 <ZyMax )
   	
  ix0= (Zx0- ZxMin)/PixelWidth; 
  iy0 = (ZyMax - Zy0)/PixelHeight; // inverse Y axis 
  ix1= (Zx1- ZxMin)/PixelWidth; 
  iy1= (ZyMax - Zy1)/PixelHeight; // inverse Y axis 
   	
  // second step of clipping
  if (ix0 >=ixMin && ix0<=ixMax && ix0 >=ixMin && ix0<=ixMax && iy0 >=iyMin && iy0<=iyMax && iy1 >=iyMin && iy1<=iyMax )
    iDrawLine(ix0,iy0,ix1,iy1,color, array) ;

  return 0;
}




int DrawAttractors(const complex double zpp[],  double p_size, unsigned char A[]){

	unsigned char color = 0;
	
	
  	
  	// join points by lin to create closed curve
  	for (int i=0;i<period ; ++i){
  		dDrawLineSegment(zpp[i], zpp[i+1],color,A);}
  	dDrawLineSegment(zpp[period-1], zpp[0],color,A);
  	// 
  	PlotAllPoints(zpp, period, color, p_size,A);
	
	
	

  return 0;

}




int MarkTraps(unsigned char A[]){

  unsigned int ix, iy;		// pixel coordinate 
  unsigned int i;


  fprintf (stderr, "Mark traps \n");
  // for all pixels of image 
#pragma omp parallel for schedule(dynamic) private(ix,iy) shared(A, ixMax , iyMax)
  for (iy = iyMin; iy <= iyMax; ++iy)
    {
      fprintf (stderr, " %d from %d \r", iy, iyMax);	//info 
      for (ix = ixMin; ix <= ixMax; ++ix){
	if (IsInsideTrap(ix, iy)) {
	  	i= Give_i(ix,iy); /* index of _data array */
	  	A[i]= 255;// 255-A[i]; // inverse color
	}}}
  return 0;
}









// ***********************************************************************************************
// ********************** mark immediate basin of attracting cycle***************************************
// ***************************************************************************************************


int FillContour(complex double seed,  unsigned char color, unsigned char _data[])
{ 
  /* 
     fills contour with black border ( color = iColorOfBoundary)  using seed point inside contour 
     and horizontal lines 
     it starts from seed point, saves max right( iXmaxLocal) and max left ( iXminLocal) interior points of horizontal line,
     in new line ( iY+1 or iY-1) it computes new interior point  : iXmidLocal=iXminLocal + (iXmaxLocal-iXminLocal)/2;
     result is stored in _data array : 1D array of 1-bit colors ( shades of gray)
     it does not check if index of _data array is good  so memory error is possible 
     
     it need array with components boundaries mrked by iColorOfBoundary
     
  */
	double dXseed = creal(seed);
	double dYseed = cimag(seed);
	// from 
  	int iXseed = (int)((dXseed - ZxMin)/PixelWidth);
  	int iYseed = (int)((ZyMax - dYseed )/PixelHeight); // reversed Y axis
  	
  	
  	
  	
  	int iX; /* seed integer coordinate */
    	int iY = iYseed;
    	/* most interior point of line iY */
    	int iXmidLocal=iXseed;
    	/* min and max of interior points of horizontal line iY */
    	int iXminLocal;
    	int iXmaxLocal; 
  	int i ; /* index of _data array */;


	//fprintf (stderr, "FillContour seed = %.16f %+.16f = %d %+d\n",creal(seed), cimag(seed), iXseed,iYseed);
  
  	/* ---------  move up --------------- */ 
  do{
    iX=iXmidLocal;
    i =Give_i(iX,iY); /* index of _data array */;
  
    /* move to right */
    while (_data[i] != iColorOfBoundary) 
      { _data[i]=color;
	iX+=1; 
	i=Give_i(iX,iY);  
      }
    iXmaxLocal=iX-1;

    /* move to left */
    iX=iXmidLocal-1; 
    i=Give_i(iX,iY);
    while (_data[i] != iColorOfBoundary) 
      { _data[i]=color;
	iX-=1; 
	i=Give_i(iX,iY); 
      }
    iXminLocal=iX+1; 

    iY+=1; /* move up */
    iXmidLocal=iXminLocal + (iXmaxLocal-iXminLocal)/2; /* new iX inside contour */
    i=Give_i(iXmidLocal,iY); /* index of _data array */;
    if ( _data[i] == iColorOfBoundary)  break; /*  it should not cross the border */
 
  } while  (iY<iyMax); 
  
  
  /* ------  move down ----------------- */
  iXmidLocal=iXseed;
  iY=iYseed-1;
  
  
  do{
    iX=iXmidLocal;
    i =Give_i(iX,iY); /* index of _data array */;
  
    /* move to right */
    while (_data[i] != iColorOfBoundary) /*  */
      { _data[i]=color;
	iX+=1;
	i=Give_i(iX,iY);  
      }
    iXmaxLocal=iX-1;

    /* move to left */
    iX=iXmidLocal-1; 
    i=Give_i(iX,iY);
    while (_data[i] != iColorOfBoundary) /*  */
      { _data[i]=color;
	iX-=1; /* move to right */
	i=Give_i(iX,iY);  
      }
    iXminLocal=iX+1; 
  
    iY-=1; /* move down */
    iXmidLocal=iXminLocal + (iXmaxLocal-iXminLocal)/2; /* new iX inside contour */
    i=Give_i(iXmidLocal,iY); /* index of _data array */;
    if ( _data[i]== iColorOfBoundary)  break; /*  it should not cross the border */
    
  } while  (0<iY); 

	//fprintf (stderr, "FillContour done \n");
  return 0;
}



// needs zpp and period global var
int  MarkImmediateBasin( unsigned char A[]){
	fprintf (stderr, "mark immediate basin of attracting cycle \n");
	
	//printf("  \n");
	unsigned char iColor = 100;
	
	
	if (period==1)
		{ FillContour(zp0, iColor , A);}
	//for (int i=0;i<period ; ++i){
  		// FillContour(zpp[i], iColor , A);
  		 
  		//}
 	return 0;
 	}

















// ***********************************************************************************************
// ********************** edge detection usung Sobel filter ***************************************
// ***************************************************************************************************

// from Source to Destination
int ComputeBoundaries(unsigned char S[], unsigned char D[])
{
 
  unsigned int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  unsigned int i; /* index of 1D array  */
  /* sobel filter */
  unsigned char G, Gh, Gv; 
  // boundaries are in D  array ( global var )
 
  // clear D array
  memset(D, iColorOfBasin1, iSize*sizeof(*D)); // for heap-allocated arrays, where N is the number of elements = FillArrayWithColor(D , iColorOfBasin1);
 
  // printf(" find boundaries in S array using  Sobel filter\n");   
#pragma omp parallel for schedule(dynamic) private(i,iY,iX,Gv,Gh,G) shared(iyMax,ixMax)
  for(iY=1;iY<iyMax-1;++iY){ 
    for(iX=1;iX<ixMax-1;++iX){ 
      Gv= S[Give_i(iX-1,iY+1)] + 2*S[Give_i(iX,iY+1)] + S[Give_i(iX-1,iY+1)] - S[Give_i(iX-1,iY-1)] - 2*S[Give_i(iX-1,iY)] - S[Give_i(iX+1,iY-1)];
      Gh= S[Give_i(iX+1,iY+1)] + 2*S[Give_i(iX+1,iY)] + S[Give_i(iX-1,iY-1)] - S[Give_i(iX+1,iY-1)] - 2*S[Give_i(iX-1,iY)] - S[Give_i(iX-1,iY-1)];
      G = sqrt(Gh*Gh + Gv*Gv);
      i= Give_i(iX,iY); /* compute index of 1D array from indices of 2D array */
      if (G==0) {D[i]=255;} /* background */
      else {D[i]=iColorOfBoundary;}  /* boundary */
    }
  }
 
   
 
  return 0;
}



// copy from Source to Destination
int CopyBoundaries(unsigned char S[],  unsigned char D[])
{
 
  unsigned int iX,iY; /* indices of 2D virtual array (image) = integer coordinate */
  unsigned int i; /* index of 1D array  */
 
 
  //printf("copy boundaries from S array to D array \n");
  for(iY=1;iY<iyMax-1;++iY)
    for(iX=1;iX<ixMax-1;++iX)
      {i= Give_i(iX,iY); if (S[i]==iColorOfBoundary) D[i]=iColorOfBoundary;}
 
 
 
  return 0;
}







// FillAllArrayWithColor
//memset (data, 255, sizeof (unsigned char ) * iSize);








// *******************************************************************************************
// ********************************** save A array to pgm file ****************************
// *********************************************************************************************

int SaveArray2PGMFile (unsigned char A[],  char * n, char *comment)
{

  FILE *fp;
  const unsigned int MaxColorComponentValue = 255;	/* color component is coded from 0 to 255 ;  it is 8 bit color file */
  char name[100];		/* name of file */
  snprintf (name, sizeof name, "%s", n);	/* radius and iHeght are global variables */
  char *filename = strcat (name, ".pgm");
  char long_comment[200];
  sprintf (long_comment, "Numerical approximation of Julia set for f(z)= z^2 + c %s", comment);





  // save image array to the pgm file 
  fp = fopen (filename, "wb");	// create new file,give it a name and open it in binary mode 
  fprintf (fp, "P5\n # %s\n %u %u\n %u\n", long_comment, iWidth, iHeight, MaxColorComponentValue);	// write header to the file
  size_t rSize = fwrite (A, sizeof(A[0]), iSize, fp);	// write whole array with image data bytes to the file in one step 
  fclose (fp);

  // info 
  if ( rSize == iSize) 
    {
      printf ("File %s saved ", filename);
      if (long_comment == NULL || strlen (long_comment) == 0)
	printf ("\n");
      else { printf (". Comment = %s \n", long_comment); }
    }
  else {printf("wrote %zu elements out of %llu requested\n", rSize,  iSize);}

  return 0;
}




int PrintCInfo ()
{

  printf ("gcc version: %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);	// https://stackoverflow.com/questions/20389193/how-do-i-check-my-gcc-c-compiler-version-for-my-eclipse
  // OpenMP version is displayed in the console : export  OMP_DISPLAY_ENV="TRUE"

  printf ("__STDC__ = %d\n", __STDC__);
  printf ("__STDC_VERSION__ = %ld\n", __STDC_VERSION__);
  printf ("c dialect = ");
  switch (__STDC_VERSION__)
    {				// the format YYYYMM 
    case 199409L:
      printf ("C94\n");
      break;
    case 199901L:
      printf ("C99\n");
      break;
    case 201112L:
      printf ("C11\n");
      break;
    case 201710L:
      printf ("C18\n");
      break;
      //default : /* Optional */

    }

  return 0;
}


int
PrintProgramInfo ()
{


 	// display info messages
  	fprintf (stdout, "Numerical approximation of dynamic plane with Julia set for f(z)= z^3 + c  \n");
  	fprintf (stdout, "c =  %.16f %+.16f*i  is the parabolic parameter with internal angle 0\n", creal (c), cimag (c) );
  	fprintf(stdout, "critical point zcr =  %.16f %+.16f*i \n", creal (zcr), cimag (zcr));

  	fprintf (stdout, "Image Width = %f in world coordinate\n", ZxMax - ZxMin);
  	fprintf (stdout, "PixelWidth = %.16f \n", PixelWidth);
  	


  	fprintf (stdout, "plane description \n");
  	fprintf (stdout, "plane center z =  %.16f %+.16f*i  and radius = %.16f \n", creal (plane_center), cimag (plane_center), radius);
  	// center and radius
  	// center and zoom
  	// GradientRepetition
  	fprintf (stdout, "Maximal number of iterations = iterMax = %d \n", IterMax);
  	fprintf (stdout, "ratio of image  = %f ; it should be 1.000 ...\n", ratio);
  
   	fprintf (stdout, "sizes of traps around attractors \n");
  	fprintf (stdout, "Escaping Radius = ER = %.16f = %f *PixelWidth = %f %% of ImageWidth \n", ER, ER / PixelWidth, ER / ZxMax - ZxMin);
	fprintf (stdout, "trap center z =  %.16f %+.16f*i  and radius = %.16f \n", creal (trap_center), cimag (trap_center), AR);
	fprintf (stdout, "Atracting Radius = AR = %.16f = %f *PixelWidth = %f %% of ImageWidth \n", AR, AR / PixelWidth, AR / ZxMax - ZxMin);
   	fprintf(stdout, "periodic cycle = parabolic fixed point z =  %.16f %+.16f*i \n", creal (zp0), cimag (zp0));
   	
  		
  	//




  return 0;
}



int SetPlane(complex double plane_center, double radius, double a_ratio){

  ZxMin = creal(plane_center) - radius*a_ratio;	
  ZxMax = creal(plane_center) + radius*a_ratio;	//0.75;
  ZyMin = cimag(plane_center) - radius;	// inv
  ZyMax = cimag(plane_center) + radius;	//0.7;
  return 0;

}



// Check Orientation of z-plane image : mark first quadrant of complex plane 
// it should be in the upper right position
// uses global var :  ...
int CheckZPlaneOrientation(unsigned char A[] )
{
 
  double Zx, Zy; //  Z= Zx+ZY*i;
  unsigned i; /* index of 1D array */
  unsigned int ix, iy;		// pixel coordinate 
	
  fprintf(stderr, "compute image CheckOrientation\n");
  // for all pixels of image 
#pragma omp parallel for schedule(dynamic) private(ix,iy, i, Zx, Zy) shared(A, ixMax , iyMax) 
  for (iy = iyMin; iy <= iyMax; ++iy){
    //fprintf (stderr, " %d from %d \r", iy, iyMax);	//info 
    for (ix = ixMin; ix <= ixMax; ++ix){
      // from screen to world coordinate 
      Zy = GiveZy(iy);
      Zx = GiveZx(ix);
      i = Give_i(ix, iy); /* compute index of 1D array from indices of 2D array */
      if (Zx>0 && Zy>0) A[i]=255-A[i];   // check the orientation of Z-plane by marking first quadrant */
    }
  }
   
   
  return 0;
}







// *****************************************************************************
//;;;;;;;;;;;;;;;;;;;;;;  setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
// **************************************************************************************

int setup ()
{

  fprintf (stderr, "setup start\n");



	c = 2.0/(3.0*sqrt(3.0)); // https://arxiv.org/pdf/1510.00043v2.pdf


  /* 2D array ranges */

  iWidth = iHeight* DisplayAspectRatio ;
  iSize = iWidth * iHeight;	// size = number of points in array 
  // iy
  iyMax = iHeight - 1;		// Indexes of array starts from 0 not 1 so the highest elements of an array is = array_name[size-1].
  //ix

  ixMax = iWidth - 1;

  /* 1D array ranges */
  // i1Dsize = i2Dsize; // 1D array with the same size as 2D array
  iMax = iSize - 1;		// Indexes of array starts from 0 not 1 so the highest elements of an array is = array_name[size-1].

  
  SetPlane( plane_center, radius,  DisplayAspectRatio );	
  /* Pixel sizes */
  PixelWidth = (ZxMax - ZxMin) / ixMax;	//  ixMax = (iWidth-1)  step between pixels in world coordinate 
  PixelHeight = (ZyMax - ZyMin) / iyMax;
  ratio = ((ZxMax - ZxMin) / (ZyMax - ZyMin)) / ((double) iWidth / (double) iHeight);	// it should be 1.000 ...
  
   // compute fixed point 
    
   zp0 = 1.0; // 1.0/ sqrt(3.0); //  
  
  // LSM  
  // Escape Radius ( of circle around infinity 
   ER = 200.0; // 
  ER2 = ER*ER;
  // Attraction Radius
  AR = GiveTunedAR(100); // quality of the image ,compute  after setting zp0 value !!!
  AR2 = AR * AR;
  trap_center =  zp0 - AR;
  trap_center2 = zp0 + AR;
  
  
  
  
  
 
  
  
  
  // for MBD
  t0 = 0.0; // Is it iternal angle from inetrnal adress  ???
  
  
  // DEM
 // BoundaryWidth = 0.5*iWidth/2000.0  ; //  measured in pixels ( when iWidth = 2000) 
  //distanceMax = BoundaryWidth*PixelWidth;



  /* create dynamic 1D arrays for colors ( shades of gray ) */
  	data = malloc (iSize * sizeof (unsigned char));

  	edge = malloc (iSize * sizeof (unsigned char));
 	edge2 = malloc (iSize * sizeof (unsigned char));
  if (data == NULL || edge == NULL || edge2 == NULL)
    {
      fprintf (stderr, " Could not allocate memory");
      return 1;
    }
  




 


  fprintf (stderr, " end of setup \n");

  return 0;

}				// ;;;;;;;;;;;;;;;;;;;;;;;;; end of the setup ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




int end ()
{


  fprintf (stderr, " allways free memory (deallocate )  to avoid memory leaks \n");	// https://en.wikipedia.org/wiki/C_dynamic_memory_allocation
  free (data);
  free(edge);
 	free(edge2);

  PrintProgramInfo ();
  PrintCInfo ();
  return 0;

}

// ********************************************************************************************************************
/* -----------------------------------------  main   -------------------------------------------------------------*/
// ********************************************************************************************************************

int main ()
{
  setup ();
  
  
   
  
  
    DrawImage (data, FatouBasins);	 
    SaveArray2PGMFile (data,  "FatouBasins" , "FatouBasins ");
    
    ComputeBoundaries(data,edge);
    SaveArray2PGMFile (edge,  "FatouBasins_LCM" , "FatouBasins_LCM ");
    
    CopyBoundaries(edge, data);
   SaveArray2PGMFile (data,  "FatouBasins_LSCM" , "FatouBasins_LSCM");
  
   MarkTraps(data);
   PlotBigPoint(zp0, 255, 10.0, data);
  SaveArray2PGMFile (data,  "FatouBasins_LSCM_trap" , "FatouBasins_LSCM_trap");
  
  
  //DrawImage (data, FatouComponents);	 
  //SaveArray2PGMFile (data,  "FatouComponents" , "FatouComponents ");
    
    
  //ComputeBoundaries(data,edge);
  //SaveArray2PGMFile (edge,  "FatouComponents_LCM" , "FatouComponents_LCM ");
    
  //CopyBoundaries(edge, data);
  //SaveArray2PGMFile (data,  "FatouComponents_LSCM" , "FatouComponents_LSCM");
    
    
 // MarkTraps(data);
  //SaveArray2PGMFile (data,  "FatouComponents_LSCM_trap" , "FatouComponents_LSCM_trap");
   
  
 //MarkImmediateBasin( edge);
  //DrawAttractors(zpp, 10, edge);
  //SaveArray2PGMFile (edge,  "imm" , "imm"); 
 
 
 
 //PlotBigPoint(trap_center, 0, 10.0, data);
 // DrawAttractors(zpp, 10, data);
 //SaveArray2PGMFile (data,  "FatouBasins_LSCM_zp" , "FatouComponents_LSCM_zp");
     
    
    
 
  DrawImage (data, LSM);	
 SaveArray2PGMFile (data,  "LSM" , "LSM");
  
 ComputeBoundaries(data,edge2);
  SaveArray2PGMFile (edge2,  "LCM" , "LCM ");
    
  CopyBoundaries(edge2, data);
 SaveArray2PGMFile (data,  "LSCM" , "LSCM");
 
 
 
    
   DrawImage (data, LS2M);	
 SaveArray2PGMFile (data,  "LS2M" , "LS2M");
  
    
  CopyBoundaries(edge2, data);
 SaveArray2PGMFile (data,  "LS2CM" , "LS2CM");
 
 
 //int DrawForwardOrbit(const complex double z0, const unsigned long long int i_Max,unsigned char iColor, double p_size, unsigned char A[]){
 DrawForwardOrbit(zcr, IterMax, 0, 20.0, data);
     SaveArray2PGMFile (data,  "LS2CM_cr" , "LS2CM_cr");
    
      DrawImage (data, MBD);	
 SaveArray2PGMFile (data,  "MBD" , "MBD");
 
 ComputeBoundaries(data,edge);
  SaveArray2PGMFile (edge,  "MBD_LCM" , "MBD_LCM ");
  
  
  
  
    
  CopyBoundaries(edge, data);
  SaveArray2PGMFile (data,  "MBD_LSCM" , "MBD_LSCM");
  
  CopyBoundaries(edge2, edge);
  SaveArray2PGMFile (edge,  "MBD_LSM_LCM" , "MBD_LSM_LCM ");
    
    CopyBoundaries(edge, data);
  SaveArray2PGMFile (data,  "MBD_LSM_LSCM" , "MBD_LSM_LSCM");
  
  
       DrawImage (data, IBD);	
 SaveArray2PGMFile (data,  "IBD" , "IBD");
 
 ComputeBoundaries(data,edge);
  SaveArray2PGMFile (edge,  "IBD_LCM" , "IBD_LCM ");
  
  
  
  
    
  CopyBoundaries(edge, data);
  SaveArray2PGMFile (data,  "IBD_LSCM" , "IBD_LSCM");
  
  CopyBoundaries(edge2, edge);
  SaveArray2PGMFile (edge,  "IBD_LSM_LCM" , "IBD_LSM_LCM ");
    
    CopyBoundaries(edge, data);
  SaveArray2PGMFile (data,  "IBD_LSM_LSCM" , "IBD_LSM_LSCM");
  
  
  
   DrawImage (data, ParabolicCheckerboard);	
 SaveArray2PGMFile (data,  "ParabolicCheckerboard_LSM" , "ParabolicCheckerboard_LSM");
  
 //ComputeBoundaries(data,edge2);
  //SaveArray2PGMFile (edge2,  "ParabolicCheckerboard_LCM" , "ParabolicCheckerboard_LCM ");
    
  CopyBoundaries(edge2, data);
 SaveArray2PGMFile (data,  "ParabolicCheckerboard_LSCM" , "ParabolicCheckerboard_LSCM");
 
  
   DrawImage (data, ParabolicCheckerboard2);	
 SaveArray2PGMFile (data,  "ParabolicCheckerboard2_LSM" , "ParabolicCheckerboard2_LSM");
 
 
 ComputeBoundaries(data,edge);
  SaveArray2PGMFile (edge,  "ParabolicCheckerboard2_LCM" , "ParabolicCheckerboard2_LCM ");
    
  CopyBoundaries(edge2, data);
 SaveArray2PGMFile (data,  "ParabolicCheckerboard2_LSCM2" , "ParabolicCheckerboard2_LSCM2");
  
 
   
  end ();

  return 0;
}
