/* estimate_variability.cpp as part of variability_estimator_lsd0.9
this code is used for determining multi-color structure functions of QSOs
*/

#include <ctime>  
//#include "fitsio.h"
#include <iostream>
#include <time.h>
#include <set>
#include <iostream>
#include <algorithm>
#include <vector>
#include <iterator>

#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_sort.h> 
#include <gsl/gsl_permute.h> 
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/param.h>

#include <iostream>
#include <cmath>

#include <cstdlib>      // std::rand, std::srand
#include <ctime>
#include <fstream>

#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <iomanip> 
#include <string>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h> 
#include <time.h>
#include <math.h>
#include <cmath>
/* GSL - GNU Scientific Library  */
/* #define GSL_CHECK_RANGE_OFF   */
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#include <gsl/gsl_cdf.h>
/* ----------------------------- */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
/* ----------------------------- */
#include <gsl/gsl_statistics_double.h>

#include <algorithm>    // std::random_shuffle



using namespace std;



//    rad - convert degrees to rad
double rad(double x)
{
 double pi=3.141592654;
 double radvalue=x*2.0*pi/360.0;

 return radvalue;
}

//     deg - convert rad to degrees
double deg (double x)
{
  double pi=3.141592654;
  double degvalue = x*360.0/(2.0*pi);

  return degvalue;
}

//     getb - get Galactic Latitude
double getb(double dgc, double rgc,double d, double r)
{
  
 return asin(sin(dgc)*sin(d)+cos(dgc)*cos(d)*(cos(r-rgc)));
}

//    getl - get Galactic Longitude
double getl(double dgc, double rgc,double lcp, double d, double r, double b)
{
  double pi=3.141592654;

  double sin_t,cos_t,t,l;
  double getl;

  //     Calculate sin(theta) and cos(theta)

  sin_t=(cos(d)*sin(r-rgc)/cos(b));
  cos_t=((cos(dgc)*sin(d)-sin(dgc)*cos(d)*cos(r-rgc))/cos(b));

  //     Use signs of sin_t and cos_t to determine which quadrant t lies in

  if((sin_t>=0) && (cos_t>0))
  {
     t=asin(sin_t);
  }
  
  if((sin_t>=0) && (cos_t<0))
  {
     t=pi-asin(sin_t);
  }
  
  if((sin_t<0) && (cos_t<0))
  {
     t=pi-asin(sin_t);
  }
  if((sin_t<0) && (cos_t>=0))
  {
     t=2.0*pi-asin(sin_t);
  }
  
  //     Convert t to l

  l=lcp-t;

  //     Ensure l lies in the range 0-2pi 

  if(l<0) 
  {
    l=l+(2.0*pi);
  }
  if(l>(2.0*pi)) 
  {
    l=l-(2.0*pi);
  }
  
  
  return l;
}


void radecTolb(double ra, double dec, double& longitude_l, double& latitude_b)
{
  
    double ragc=192.85948; /* RA of North Gal Pole, 2000 */
    double decgc=27.12825; /* Dec of North Gal Pole, 2000 */
    double lcp=122.932; /* longitude of the North celestial pole */
  
    double pi=3.141592654;

 double decr=0;
 double rar=0;
 double decgcr=0;
 double ragcr=0;
 double lcpr=0;
 double b=0;
 double l=0;
 longitude_l=0;
 latitude_b=0;

	/*
11.46<RA<13.7 , -1<DEC<1.24
	 */

      //calculate l and b
      
      
      
      	    
  decr=rad(dec);
    rar=rad(ra);
    decgcr=rad(decgc);
    ragcr=rad(ragc);
    lcpr=rad(lcp);

     //   CALCULATE GALACTIC COORDINATES 

     b=getb(decgcr,ragcr,decr,rar); //in decimal where 360 degree = 2 PI
     l=getl(decgcr,ragcr,lcpr,decr,rar,b); //in decimal where 360 degree = 2 PI
    
	    longitude_l =deg(l);
	    latitude_b = deg(b);
	    
	   // std::cout<<"b "<<b<<" l "<<l<<" longitude_l "<<longitude_l<<" latitude_b "<<latitude_b<<std::endl;
      }


string intToString(int value){

  stringstream ss; 
  ss << value;
  return (ss.str());
}


string doubleToString(double d){
  
  stringstream ss; 
  ss << d;
  return (ss.str());
}

double stringToDouble(string s){

  std::istringstream stream (s);
  double t;
  stream >> t;
  return t;
}
    
int stringToInt(string s){

  std::istringstream stream (s);
  int t;
  stream >> t;
  return t;
}
    
    

//solve v=Ax for x
int lusolve(gsl_matrix *A,gsl_vector *v, gsl_vector *&x)
{
  double ldet;
  int signum;
  gsl_permutation *p = gsl_permutation_alloc(A->size1);
  
  gsl_matrix *tmpA = gsl_matrix_alloc(A->size1,A->size2);
  gsl_matrix_memcpy(tmpA , A);
 
  gsl_linalg_LU_decomp(tmpA , p , &signum);
  gsl_linalg_LU_solve (tmpA,  p, v, x);
  
  gsl_permutation_free(p);
  gsl_matrix_free(tmpA);

  return 0;
} 


int inv(gsl_matrix *&A) {

  double ldet;
  int signum;
  gsl_permutation *p = gsl_permutation_alloc(A->size1);
  
  gsl_matrix *tmpA = gsl_matrix_alloc(A->size1, A->size2);
  gsl_matrix_memcpy(tmpA , A);
 
  gsl_linalg_LU_decomp(tmpA , p , &signum);
  gsl_linalg_LU_invert (tmpA, p, A);
 
  gsl_permutation_free(p);
  gsl_matrix_free(tmpA);
   
  return 0;
}


 double det(gsl_matrix *A) {

  double determinant;
  int signum;
  gsl_permutation *p = gsl_permutation_alloc(A->size1);
  
  gsl_matrix *tmpA = gsl_matrix_alloc(A->size1, A->size2);
  gsl_matrix_memcpy(tmpA , A);
 
  gsl_linalg_LU_decomp(tmpA , p , &signum);
  determinant = gsl_linalg_LU_det(tmpA , signum);
  gsl_permutation_free(p);
  gsl_matrix_free(tmpA);
   
  return determinant;
}



int singular (gsl_matrix *A)
{
  int i;
  int n = A->size1;
  
  for (i = 0; i < n; i++)
  {
     double u = gsl_matrix_get (A, i, i);
     if (u == 0) return 1;
  }
  
  if(det(A)==0)
  {return 1;}
  else
  {
	 
 return 0;
  }
}




 /* 
This function computes the logarithm of the absolute value of the determinant of a matrix A,
 $\ln|\det(A)|$, from its LU decomposition, LU. This function may be useful if the direct computation
of the determinant would overflow or underflow.
*/
double logdet(gsl_matrix *A)
{
  double ldet;
  int signum;
  gsl_permutation *p = gsl_permutation_alloc(A->size1);
  
  gsl_matrix *tmpA = gsl_matrix_alloc(A->size1, A->size2);
  gsl_matrix_memcpy(tmpA , A);
 
  gsl_linalg_LU_decomp(tmpA , p , &signum);   
  ldet = gsl_linalg_LU_lndet (tmpA);
    
  gsl_permutation_free(p);
  gsl_matrix_free(tmpA);
   
  return ldet;
}
 


int goodObservationFlag(int flag/*, string &badflag*/)
{
  
  
  /*
   badflags = (ps.flags['FITFAIL'] | ps.flags['POORFIT'] |
            ps.flags['SATSTAR'] | ps.flags['BLEND'] |
            ps.flags['BADPSF'] | ps.flags['DEFECT'] | ps.flags['SATURATED'] |
            ps.flags['CR_LIMIT'] | ps.flags['EXT_LIMIT'] |
            ps.flags['MOMENTS_FAILURE'] | ps.flags['SKY_FAILURE'] |
            ps.flags['SKYVAR_FAILURE'] | ps.flags['BIG_RADIUS'] |
            ps.flags['SIZE_SKIPPED'] | ps.flags['ON_SPIKE'] |
            ps.flags['ON_GHOST'] | ps.flags['OFF_CHIP'])
   */
  
  // Bit flags to distinguish analysis results
	// When adding to or subtracting from this list, please also modify pmSourceMaskHeader
	typedef enum {
	    PM_SOURCE_MODE_DEFAULT                = 0x00000000, ///< Initial value: resets all bits
	    PM_SOURCE_MODE_PSFMODEL               = 0x00000001, ///< Source fitted with a psf model (linear or non-linear)
	    PM_SOURCE_MODE_EXTMODEL               = 0x00000002, ///< Source fitted with an extended-source model
	    PM_SOURCE_MODE_FITTED                 = 0x00000004, ///< Source fitted with non-linear model (PSF or EXT; good or bad)
	    PM_SOURCE_MODE_FAIL                   = 0x00000008, ///< Fit (non-linear) failed (non-converge, off-edge, run to zero)
	    PM_SOURCE_MODE_POOR                   = 0x00000010, ///< Fit succeeds, but low-SN, high-Chisq, or large (for PSF -- drop?)
	    PM_SOURCE_MODE_PAIR                   = 0x00000020, ///< Source fitted with a double psf
	    PM_SOURCE_MODE_PSFSTAR                = 0x00000040, ///< Source used to define PSF model
	    PM_SOURCE_MODE_SATSTAR                = 0x00000080, ///< Source model peak is above saturation
	    PM_SOURCE_MODE_BLEND                  = 0x00000100, ///< Source is a blend with other sources
	    PM_SOURCE_MODE_EXTERNAL               = 0x00000200, ///< Source based on supplied input position
	    PM_SOURCE_MODE_BADPSF                 = 0x00000400, ///< Failed to get good estimate of object's PSF
	    PM_SOURCE_MODE_DEFECT                 = 0x00000800, ///< Source is thought to be a defect
	    PM_SOURCE_MODE_SATURATED              = 0x00001000, ///< Source is thought to be saturated pixels (bleed trail)
	    PM_SOURCE_MODE_CR_LIMIT               = 0x00002000, ///< Source has crNsigma above limit
	//    PM_SOURCE_MODE_EXT_LIMIT              = 0x00004000, ///< Source has extNsigma above limit
	    PM_SOURCE_MODE_MOMENTS_FAILURE        = 0x00008000, ///< could not measure the moments
	    PM_SOURCE_MODE_SKY_FAILURE            = 0x00010000, ///< could not measure the local sky
	    PM_SOURCE_MODE_SKYVAR_FAILURE         = 0x00020000, ///< could not measure the local sky variance
	    PM_SOURCE_MODE_BELOW_MOMENTS_SN       = 0x00040000, ///< moments not measured due to low S/N
	    PM_SOURCE_MODE_BIG_RADIUS             = 0x00100000, ///< poor moments for small radius, try large radius
	    PM_SOURCE_MODE_AP_MAGS                = 0x00200000, ///< source has an aperture magnitude
	    PM_SOURCE_MODE_BLEND_FIT              = 0x00400000, ///< source was fitted as a blend
	    PM_SOURCE_MODE_EXTENDED_FIT           = 0x00800000, ///< full extended fit was used
	    PM_SOURCE_MODE_EXTENDED_STATS         = 0x01000000, ///< extended aperture stats calculated
	    PM_SOURCE_MODE_LINEAR_FIT             = 0x02000000, ///< source fitted with the linear fit
	    PM_SOURCE_MODE_NONLINEAR_FIT          = 0x04000000, ///< source fitted with the non-linear fit
	    PM_SOURCE_MODE_RADIAL_FLUX            = 0x08000000, ///< radial flux measurements calculated
	    PM_SOURCE_MODE_SIZE_SKIPPED           = 0x10000000, ///< size could not be determined
	    PM_SOURCE_MODE_ON_SPIKE               = 0x20000000, ///< peak lands on diffraction spike
	    PM_SOURCE_MODE_ON_GHOST               = 0x40000000, ///< peak lands on ghost or glint
	    PM_SOURCE_MODE_OFF_CHIP               = 0x80000000, ///< peak lands off edge of chip
	} pmSourceMode;

 
	//add them to badflag
 
int counter=0;
  //one has to be set
  if ((flag & (PM_SOURCE_MODE_FAIL))
         == (PM_SOURCE_MODE_FAIL ))
  {
    //badflag +=";PM_SOURCE_MODE_FAIL";
    counter++;
  }
  if ((flag & (PM_SOURCE_MODE_POOR))
         == (PM_SOURCE_MODE_POOR ))
  {
    //badflag+=";PM_SOURCE_MODE_POOR";
    counter++;
  }
  if ((flag & (PM_SOURCE_MODE_SATSTAR))
         == (PM_SOURCE_MODE_SATSTAR ))
  {
    //badflag +=";PM_SOURCE_MODE_SATSTAR";
    counter++;
  }
  if ((flag & (PM_SOURCE_MODE_BLEND))
         == (PM_SOURCE_MODE_BLEND ))
  {
    //badflag+=";PM_SOURCE_MODE_BLEND";
    counter++;
  }
  if ((flag & (PM_SOURCE_MODE_BADPSF))
         == (PM_SOURCE_MODE_BADPSF ))
  {
    //badflag+=";PM_SOURCE_MODE_BADPSF";
    counter++;
  }
  if ((flag & (PM_SOURCE_MODE_DEFECT))
         == (PM_SOURCE_MODE_DEFECT ))
  {
    //badflag+=";PM_SOURCE_MODE_DEFECT";
    counter++;
  }
  if ((flag & (PM_SOURCE_MODE_SATURATED))
         == (PM_SOURCE_MODE_SATURATED ))
  {
    //badflag+=";PM_SOURCE_MODE_SATURATED";
     counter++;
  }
  if ((flag & (PM_SOURCE_MODE_CR_LIMIT))
         == (PM_SOURCE_MODE_CR_LIMIT ))
  {
   // badflag+=";PM_SOURCE_MODE_CR_LIMIT";
    counter++;
  }
  /*if ((flag & (PM_SOURCE_MODE_EXT_LIMIT))
         == (PM_SOURCE_MODE_EXT_LIMIT ))
  {
    badflag+=";PM_SOURCE_MODE_EXT_LIMIT";
    counter++;
  }*/
 if ((flag & (PM_SOURCE_MODE_MOMENTS_FAILURE))
         == (PM_SOURCE_MODE_MOMENTS_FAILURE ))
  {
    //badflag+=";PM_SOURCE_MODE_MOMENTS_FAILURE";
     counter++;
  }
  if ((flag & (PM_SOURCE_MODE_SKY_FAILURE))
         == (PM_SOURCE_MODE_SKY_FAILURE ))
  {
    //badflag+=";PM_SOURCE_MODE_SKY_FAILURE";
     counter++;
  }
  if ((flag & (PM_SOURCE_MODE_SKYVAR_FAILURE))
         == (PM_SOURCE_MODE_SKYVAR_FAILURE ))
  {
    //badflag+=";PM_SOURCE_MODE_SKYVAR_FAILURE";
     counter++;
  }
  if ((flag & (PM_SOURCE_MODE_BIG_RADIUS))
         == (PM_SOURCE_MODE_BIG_RADIUS ))
  {
    //badflag+=";PM_SOURCE_MODE_BIG_RADIUS";
    counter++;
  }
  if ((flag & (PM_SOURCE_MODE_SIZE_SKIPPED))
         == (PM_SOURCE_MODE_SIZE_SKIPPED ))
  {
    //badflag+=";PM_SOURCE_MODE_SIZE_SKIPPED";
     counter++;
  }
  if ((flag & (PM_SOURCE_MODE_ON_SPIKE))
         == (PM_SOURCE_MODE_ON_SPIKE ))
  {
    //badflag+=";PM_SOURCE_MODE_ON_SPIKE";
    counter++;
  }
  if ((flag & (PM_SOURCE_MODE_ON_GHOST))
         == (PM_SOURCE_MODE_ON_GHOST ))
  {
    //badflag+=";PM_SOURCE_MODE_ON_GHOST";
    counter++; 
  }
  if ((flag & (PM_SOURCE_MODE_OFF_CHIP))
         == (PM_SOURCE_MODE_OFF_CHIP ))
  {
    //badflag+=";PM_SOURCE_MODE_OFF_CHIP";
     counter++;
  }/*
  if(badflag.size()>0)
  {
    if (badflag.substr(0,1)==";")
    { 
      badflag = badflag.substr (1, badflag.size()-1);
    }
  }*/
  
  
  if(counter>0)
  {return -1;}
  else
  {return 0;
  }
   
  
   
}




void applyColorcut(gsl_vector *gband_cleaned, gsl_vector *rband_cleaned, gsl_vector *iband_cleaned, gsl_vector *zband_cleaned, gsl_vector *yband_cleaned,
		       gsl_vector *gbanderr_cleaned, gsl_vector *rbanderr_cleaned, gsl_vector *ibanderr_cleaned, gsl_vector *zbanderr_cleaned, gsl_vector *ybanderr_cleaned,
                       double &errorweightedmean_g, double &errorweightedmean_r, double &errorweightedmean_i, double &errorweightedmean_z, 
		       double &errorweightedmean_y, int &basic_colorcut_passed, int &stellar_locus)
{
	  //and then calculate the error-weighted means for applying color-mag cut
	  
	    
	  double sum1=0;
	  double sum2=0;
	  double factor1=0;
	 

	  double dim_gband_cleaned=0;
	  double dim_rband_cleaned=0;
	  double dim_iband_cleaned=0;
	  double dim_zband_cleaned=0;
	  double dim_yband_cleaned=0;
	  
	  double chi2g=0;
	  double chi2r=0;
	  double chi2i=0;
	  double chi2z=0;
	  double chi2y=0;  
	
	  double value=0;
	  
	  //g
	  if(gband_cleaned!=0)
	  {
	  
	    sum1=0;
	    sum2=0;
	    factor1=0;
	    
	    for(int i=0;i<gbanderr_cleaned->size;i++)
	    {
	      sum1 += 1/ (pow(gsl_vector_get(gbanderr_cleaned,i),2));
	    }
	    factor1=1/sum1;
	    
	    
	    for(int i=0;i<gbanderr_cleaned->size;i++)
	    {
	      sum2 += gsl_vector_get(gband_cleaned,i) / (pow(gsl_vector_get(gbanderr_cleaned,i),2));
	    }
	    
	    errorweightedmean_g=factor1*sum2;
	    dim_gband_cleaned =gband_cleaned->size;
	  }
	  
	  //r
	  if(rband_cleaned!=0)
	  {
	  
	  sum1=0;
	  sum2=0;
	  factor1=0;
	  
	  for(int i=0;i<rbanderr_cleaned->size;i++)
	  {
	    sum1 += 1/ (pow(gsl_vector_get(rbanderr_cleaned,i),2));
	  }
	  factor1=1/sum1;
	  
	  
	  for(int i=0;i<rbanderr_cleaned->size;i++)
	  {
	    sum2 += gsl_vector_get(rband_cleaned,i) / (pow(gsl_vector_get(rbanderr_cleaned,i),2));
	  }
	  
	  errorweightedmean_r=factor1*sum2;
	  dim_rband_cleaned =rband_cleaned->size;
	}
	  //i
	  if(iband_cleaned!=0)
	  {
	  
	  sum1=0;
	  sum2=0;
	  factor1=0;
	  
	  for(int i=0;i<ibanderr_cleaned->size;i++)
	  {
	    sum1 += 1/ (pow(gsl_vector_get(ibanderr_cleaned,i),2));
	  }
	  factor1=1/sum1;
	  
	  
	  for(int i=0;i<ibanderr_cleaned->size;i++)
	  {
	    sum2 += gsl_vector_get(iband_cleaned,i) / (pow(gsl_vector_get(ibanderr_cleaned,i),2));
	  }
	  
	  errorweightedmean_i=factor1*sum2;
	  dim_iband_cleaned =iband_cleaned->size;
	}
	  
	  //z
	  if(zband_cleaned!=0)
	  {
	  
	    sum1=0;
	    sum2=0;
	    factor1=0;
	    
	    for(int i=0;i<zbanderr_cleaned->size;i++)
	    {
	      sum1 += 1/ (pow(gsl_vector_get(zbanderr_cleaned,i),2));
	    }
	    factor1=1/sum1;
	    
	    
	    for(int i=0;i<zbanderr_cleaned->size;i++)
	    {
	      sum2 += gsl_vector_get(zband_cleaned,i) / (pow(gsl_vector_get(zbanderr_cleaned,i),2));
	    }
	    
	    errorweightedmean_z=factor1*sum2;
	    dim_zband_cleaned =zband_cleaned->size;
	  }
	  
	  //y
	  
	  if(yband_cleaned!=0)
	  {
	    sum1=0;
	    sum2=0;
	    factor1=0;
	    
	    for(int i=0;i<ybanderr_cleaned->size;i++)
	    {
	      sum1 += 1/ (pow(gsl_vector_get(ybanderr_cleaned,i),2));
	    }
	    factor1=1/sum1;
	    
	    
	    for(int i=0;i<ybanderr_cleaned->size;i++)
	    {
	      sum2 += gsl_vector_get(yband_cleaned,i) / (pow(gsl_vector_get(ybanderr_cleaned,i),2));
	    }
	    
	    errorweightedmean_y=factor1*sum2;
	    dim_yband_cleaned =yband_cleaned->size;
	  }
	  
    // calculate \chi^ 2_{\lambda}
    
  
  
    
    //general color-mag cut for all objects, removes possibly saturated or to weak objects
    
    value=0;
    //|(g-r) - (g-i)| = |r-i| > 0.2 && g-r<1 && g-i<1
    
  basic_colorcut_passed=0;
  stellar_locus=0;
  
	/*      
	      if( ( 14 < errorweightedmean_g ) && (errorweightedmean_g< 21) 
		&&( -1.0 < errorweightedmean_g-errorweightedmean_r) && (errorweightedmean_g-errorweightedmean_r < 1.0))
	      {basic_colorcut_passed=1;}
	  */    

  if( (errorweightedmean_g - errorweightedmean_r < -0.6) 
     || (errorweightedmean_r - errorweightedmean_i < -0.6)
     || (errorweightedmean_i - errorweightedmean_z < -0.6)
     || (errorweightedmean_z - errorweightedmean_y < -0.6)   
     || (errorweightedmean_g - errorweightedmean_r > 1.1) 
     || (errorweightedmean_r - errorweightedmean_i > 4)
     || (errorweightedmean_i - errorweightedmean_z > 4)
     || (errorweightedmean_z - errorweightedmean_y > 4)
     || (errorweightedmean_g - errorweightedmean_i > 1.0) 
  )
  {
    //color cut not passed
    basic_colorcut_passed = 0;
  }
  else
  {
     if(   (errorweightedmean_g < 21)  && (errorweightedmean_g > 15) 
        && (errorweightedmean_r > 15) && (errorweightedmean_i > 15) 
        && (errorweightedmean_z > 15 ) && (errorweightedmean_y >15) 
        && (errorweightedmean_g - errorweightedmean_r < 1.3) 
        && (errorweightedmean_g - errorweightedmean_i < 1.3)  )
        {
	    //color cut passed!   
            basic_colorcut_passed=1;
        }
  }


       if( (   ( (errorweightedmean_g-errorweightedmean_i)<1.4*(errorweightedmean_g-errorweightedmean_r)+0.05  )
	 &&  ( (errorweightedmean_g-errorweightedmean_i)>1.4*(errorweightedmean_g-errorweightedmean_r)-0.05   ) 
     ) )
       {stellar_locus=1;}
    
    


}







void calculate_q(gsl_vector *gband_cleaned, gsl_vector *rband_cleaned, gsl_vector *iband_cleaned, gsl_vector *zband_cleaned, gsl_vector *yband_cleaned,
		       gsl_vector *gbanderr_cleaned, gsl_vector *rbanderr_cleaned, gsl_vector *ibanderr_cleaned, gsl_vector *zbanderr_cleaned, gsl_vector *ybanderr_cleaned,
                       double &q, double &errorweightedmean_g, double &errorweightedmean_r, double &errorweightedmean_i, 
		       double &errorweightedmean_z, double &errorweightedmean_y,
		       int &basic_colorcut_passed, int &stellar_locus)
{
     
  
  //get error-weighted means
  
	errorweightedmean_g = 0;
	 errorweightedmean_r = 0;
	errorweightedmean_i = 0;
	  errorweightedmean_z = 0;
	  errorweightedmean_y = 0;
	   basic_colorcut_passed=0;
	   stellar_locus=0;
	  
  applyColorcut(gband_cleaned, rband_cleaned, iband_cleaned, zband_cleaned, yband_cleaned,
		       gbanderr_cleaned, rbanderr_cleaned, ibanderr_cleaned, zbanderr_cleaned, ybanderr_cleaned,
                       errorweightedmean_g, errorweightedmean_r, errorweightedmean_i, errorweightedmean_z, errorweightedmean_y,
		       basic_colorcut_passed, stellar_locus);

  
    
	    //make a chi2gcomponents_vec, ... c++ vector
	     
     std::vector<double> chi2gcomponents_vec;
     std::vector<double> chi2rcomponents_vec;
     std::vector<double> chi2icomponents_vec;
     std::vector<double> chi2zcomponents_vec;
     std::vector<double> chi2ycomponents_vec;
          std::vector<double> maxchiperband_vec;
     
     
    std::vector<double> chi2gcomponentssquared_vec;
    std::vector<double> chi2rcomponentssquared_vec;
    std::vector<double> chi2icomponentssquared_vec;
    std::vector<double> chi2zcomponentssquared_vec;
    std::vector<double> chi2ycomponentssquared_vec;
    
        
    double chi2g=0;
    double chi2r=0;
    double chi2i=0;
    double chi2z=0;
    double chi2y=0;   
    double chi2source=0;
      int bandsused=0;
    int numberdatapoints=0;
    double ndof=0;
    double value=0;
    
    if(gband_cleaned!=0)
    {
      for(int i=0;i<gband_cleaned->size;i++)
      {
	value =(gsl_vector_get(gband_cleaned,i) - errorweightedmean_g) / sqrt( pow(gsl_vector_get(gbanderr_cleaned,i),2)+pow(0.015, 2));
	chi2gcomponents_vec.push_back(value);
	chi2gcomponentssquared_vec.push_back(pow(value,2));
      }
    }
    
    if(rband_cleaned!=0)
    {
      for(int i=0;i<rband_cleaned->size;i++)
      {
	value =(gsl_vector_get(rband_cleaned,i) - errorweightedmean_r) / sqrt( pow(gsl_vector_get(rbanderr_cleaned,i),2)+pow(0.015, 2));
	chi2rcomponents_vec.push_back(value);
	chi2rcomponentssquared_vec.push_back(pow(value,2));
      }
    }
    
    if(iband_cleaned!=0)
    {
      for(int i=0;i<iband_cleaned->size;i++)
      {
	value =(gsl_vector_get(iband_cleaned,i) - errorweightedmean_i) / sqrt( pow(gsl_vector_get(ibanderr_cleaned,i),2)+pow(0.015, 2));
	chi2icomponents_vec.push_back(value);
	chi2icomponentssquared_vec.push_back(pow(value,2));
      }
    }
    
    if(zband_cleaned!=0)
    {
      for(int i=0;i<zband_cleaned->size;i++)
      {
	value =(gsl_vector_get(zband_cleaned,i) - errorweightedmean_z) / sqrt( pow(gsl_vector_get(zbanderr_cleaned,i),2)+pow(0.015, 2));
	chi2zcomponents_vec.push_back(value);
	chi2zcomponentssquared_vec.push_back(pow(value,2));
      }
    }
    
    if(yband_cleaned!=0)
    {
      for(int i=0;i<yband_cleaned->size;i++)
      {
	value =(gsl_vector_get(yband_cleaned,i) - errorweightedmean_y) / sqrt( pow(gsl_vector_get(ybanderr_cleaned,i),2)+pow(0.015, 2));
	chi2ycomponents_vec.push_back(value);
	chi2ycomponentssquared_vec.push_back(pow(value,2));
      }
    }    
  
    

//sum up vector components to make chi2g, ... chi2y
    //also output them for histogram of chi^2 contributions
    
       
    chi2g=0;
    chi2r=0;
    chi2i=0;
    chi2z=0;
    chi2y=0;   
    
for(int i=0;i<chi2gcomponents_vec.size();i++)
{
  chi2g+=chi2gcomponentssquared_vec[i];
}

for(int i=0;i<chi2rcomponents_vec.size();i++)
{
   chi2r+=chi2rcomponentssquared_vec[i];
}


for(int i=0;i<chi2icomponents_vec.size();i++)
{
   chi2i+=chi2icomponentssquared_vec[i];
}


for(int i=0;i<chi2zcomponents_vec.size();i++)
{
   chi2z+=chi2zcomponentssquared_vec[i];
}

for(int i=0;i<chi2ycomponents_vec.size();i++)
{
   chi2y+=chi2ycomponentssquared_vec[i];
}
    
    
    
//     chi^2_dof = sum chi^2_lamba / (N - n_bands)
     
     chi2source=0;
    
     bandsused=0;
     numberdatapoints=0;
    if(gband_cleaned!=0)
    {
      chi2source +=  chi2g;
      bandsused++;
      numberdatapoints+=gband_cleaned->size;
    }
    
    if(rband_cleaned!=0)
    {
      chi2source += chi2r;
      bandsused++;
      numberdatapoints+=rband_cleaned->size;
    }
    
    if(iband_cleaned!=0)
    {
      chi2source += chi2i;
      bandsused++;
      numberdatapoints+=iband_cleaned->size;
    }
    
    if(zband_cleaned!=0)
    {
      chi2source += chi2z;
      bandsused++;
      numberdatapoints+=zband_cleaned->size;
    }
    
    if(yband_cleaned!=0)
    {
      chi2source += chi2y;
      bandsused++;
      numberdatapoints+=yband_cleaned->size;
    }
    
    
    ndof = numberdatapoints-bandsused;
  //for chi2/dof  
    q=chi2source;
    chi2source/=ndof;

 //for Poisson
    q-=ndof;
    q/=sqrt(2*ndof);

  
	if(chi2gcomponents_vec.size()>0)
	{
	  chi2gcomponents_vec.clear();
	}
	
	if(chi2rcomponents_vec.size()>0)
	{
	  chi2rcomponents_vec.clear();
	}
	
	if(chi2icomponents_vec.size()>0)
	{
	  chi2icomponents_vec.clear();
	}
	
	if(chi2zcomponents_vec.size()>0)
	{
	  chi2zcomponents_vec.clear();
	}
	
	if(chi2ycomponents_vec.size()>0)
	{
	  chi2ycomponents_vec.clear();
	}
	
	if(maxchiperband_vec.size()>0)
	{
	  maxchiperband_vec.clear();
	}
	
	if(chi2gcomponentssquared_vec.size()>0)
	{
	  chi2gcomponentssquared_vec.clear();
	}
	
	if(chi2rcomponentssquared_vec.size()>0)
	{
	  chi2rcomponentssquared_vec.clear();
	}
	
	if(chi2icomponentssquared_vec.size()>0)
	{
	  chi2icomponentssquared_vec.clear();
	}
	
	if(chi2zcomponentssquared_vec.size()>0)
	{
	  chi2zcomponentssquared_vec.clear();
	}
	
	if(chi2ycomponentssquared_vec.size()>0)
	{
	  chi2ycomponentssquared_vec.clear();
	}
	
    
    
    
    
    
}

	  /*
	   # psf_inst_mag psf_inst_mag_sig ra dec detra detdec filterid exptime mjd_obs x_psf y_psf airmass ap_mag psf_qf psf_qf_perfect 
	   median_fwhm_geo zp ucalmag flags flags2
	   -7.69377 0.169709 0.0554833 -1.24875 0.0554889 -1.24883 r 60 55098.5 3012.01 432.575 1.28 -7.6218 0.909623 0.909623 7.66449 0 0 3012.01 0
-8.03874 0.124927 0.0554833 -1.24875 0.0554231 -1.24887 r 60 55098.5 4645.4 1375.07 1.297 -8.38809 0.999132 0.999132 7.29549 0 0 4645.4 0
-7.70149 0.127478 0.0554833 -1.24875 0.0553861 -1.24883 r 60 55102.5 478.239 2461.89 1.122 -7.07466 0.999139 0.999139 5.17646 29.1267 21.4252 478.239 0

	   
	   ('psf_inst_mag_sig', '>f4'), # looks like it's uncertainty in the
magnitude, but in fact it is 1/(S/N), i.e., you don't need the factor
of ln(10)/2.5 to convert this to S/N.  
--> delta m = 2.5 log10(1+n/s) = 2.5 log10(1+psf_inst_mag_sig) 


	   */
	  
int cleanPS1file(gsl_vector *gband_orig, gsl_vector *rband_orig, gsl_vector *iband_orig, gsl_vector *zband_orig, gsl_vector *yband_orig,
		  gsl_vector *gbanderr_orig, gsl_vector *rbanderr_orig, gsl_vector *ibanderr_orig, gsl_vector *zbanderr_orig, 
		  gsl_vector *ybanderr_orig,
		  gsl_vector *gbandtime_orig, gsl_vector *rbandtime_orig, gsl_vector *ibandtime_orig, gsl_vector *zbandtime_orig, 
		  gsl_vector *ybandtime_orig, 
		  vector<int> gband_flags_vec_orig, vector<int> rband_flags_vec_orig, vector<int> iband_flags_vec_orig, 
		  vector<int> zband_flags_vec_orig, vector<int> yband_flags_vec_orig, 		 
		  gsl_vector *gband_psf_qf_perfect_vec,gsl_vector *rband_psf_qf_perfect_vec,gsl_vector *iband_psf_qf_perfect_vec,
		  gsl_vector *zband_psf_qf_perfect_vec,gsl_vector *yband_psf_qf_perfect_vec,
		  gsl_vector *gband_ap_mag_vec,gsl_vector *rband_ap_mag_vec,gsl_vector *iband_ap_mag_vec, gsl_vector *zband_ap_mag_vec,
		  gsl_vector *yband_ap_mag_vec,
		  gsl_vector *gband_psf_inst_mag_vec, gsl_vector *rband_psf_inst_mag_vec, gsl_vector *iband_psf_inst_mag_vec, 
		  gsl_vector *zband_psf_inst_mag_vec, gsl_vector *yband_psf_inst_mag_vec, 
		  gsl_vector *&gband_cleaned, gsl_vector *&rband_cleaned, gsl_vector *&iband_cleaned, gsl_vector *&zband_cleaned,gsl_vector *&yband_cleaned,
		  gsl_vector *&gbanderr_cleaned, gsl_vector *&rbanderr_cleaned, gsl_vector *&ibanderr_cleaned, gsl_vector *&zbanderr_cleaned,gsl_vector *&ybanderr_cleaned,
		  gsl_vector *&gbandtime_cleaned, gsl_vector *&rbandtime_cleaned, gsl_vector *&ibandtime_cleaned, 
		  gsl_vector *&zbandtime_cleaned,gsl_vector *&ybandtime_cleaned, 
		 /*int &number_gband_badflag,  int &number_rband_badflag, 
		  int &number_iband_badflag,  int &number_zband_badflag,  int &number_yband_badflag, */
		 /* string &g_badflags, string &r_badflags, string &i_badflags, string &z_badflags, string &y_badflags,*/
		  int &number_gband_out_as_gband_psf_qf_perfect_less0p95,
		  int &number_rband_out_as_rband_psf_qf_perfect_less0p95,
		  int &number_iband_out_as_iband_psf_qf_perfect_less0p95,
		  int &number_zband_out_as_zband_psf_qf_perfect_less0p95,
		  int &number_yband_out_as_yband_psf_qf_perfect_less0p95,
		  int &number_gband_out_as_gband_ap_mag,int &number_rband_out_as_rband_ap_mag,int &number_iband_out_as_iband_ap_mag,
		  int &number_zband_out_as_zband_ap_mag,int &number_yband_out_as_yband_ap_mag, int &total_beforechi2,
		  int &number_gband_outchi2, int &number_rband_outchi2, int &number_iband_outchi2, int &number_zband_outchi2, int &number_yband_outchi2,
   		  int only_by_flags )	  
{
   
  
  /*
i)   compute the individual $\chi$s
ii)  compute the median of the $\chi$s
iii) compute the interquartile range, this is a number computed by IQR = Q3 −  Q1, where Q3 and Q1 are the upper and lower quartiles
iv)  throw away the points having a $| \chi | > \chi_{\mathrm{\median}}+  3 \times IQR \times 1.349/2 $, at maximum 10 \% of all light curve points
   */
  
  gband_cleaned=0;
  gbanderr_cleaned=0;
  gbandtime_cleaned=0;
  
  rband_cleaned=0;
  rbanderr_cleaned=0;
  rbandtime_cleaned=0;
 
  iband_cleaned=0;
  ibanderr_cleaned=0;
  ibandtime_cleaned=0;
  
  zband_cleaned=0;
  zbanderr_cleaned=0;
  zbandtime_cleaned=0;
  
  yband_cleaned=0;
  ybanderr_cleaned=0;
  ybandtime_cleaned=0;
  
  gsl_vector *gband=0;
  gsl_vector *gbanderr=0;
  gsl_vector *gbandtime=0;
  gsl_vector *rband=0;
  gsl_vector *rbanderr=0;
  gsl_vector *rbandtime=0;
  gsl_vector *iband=0;
  gsl_vector *ibanderr=0;
  gsl_vector *ibandtime=0;
  gsl_vector *zband=0;
  gsl_vector *zbanderr=0;
  gsl_vector *zbandtime=0;
  gsl_vector *yband=0;
  gsl_vector *ybanderr=0;
  gsl_vector *ybandtime=0;
  
  gsl_vector *gband_cleaned_tmp=0;
  gsl_vector *gbanderr_cleaned_tmp=0;
  gsl_vector *gbandtime_cleaned_tmp=0;
  gsl_vector *rband_cleaned_tmp=0;
  gsl_vector *rbanderr_cleaned_tmp=0;
  gsl_vector *rbandtime_cleaned_tmp=0;
  gsl_vector *iband_cleaned_tmp=0;
  gsl_vector *ibanderr_cleaned_tmp=0;
  gsl_vector *ibandtime_cleaned_tmp=0;
  gsl_vector *zband_cleaned_tmp=0;
  gsl_vector *zbanderr_cleaned_tmp=0;
  gsl_vector *zbandtime_cleaned_tmp=0;
  gsl_vector *yband_cleaned_tmp=0;
  gsl_vector *ybanderr_cleaned_tmp=0;
  gsl_vector *ybandtime_cleaned_tmp=0;
 
  
 //0) remove all entries having bad flags
 
 /*
  g_badflags="";
  r_badflags="";
  i_badflags="";
  z_badflags="";
  y_badflags="";
  */
  
 /*
  string gband_badflags="";
  string rband_badflags="";
  string iband_badflags="";
  string zband_badflags="";
  string yband_badflags="";*/
  
  int goodobs;
  
  int dim_gband_orig=0;
  int dim_rband_orig=0;
  int dim_iband_orig=0;
  int dim_zband_orig=0;
  int dim_yband_orig=0;
    
  int dim_gband_goodobs=0;
  int dim_rband_goodobs=0;
  int dim_iband_goodobs=0;
  int dim_zband_goodobs=0;
  int dim_yband_goodobs=0;
  

  
  /*
  It might be that we're cutting too aggressively now.  I guess you are
  losing all of the y band points because they have large errors, and	
  therefore tend to have |ap-psf| > 0.1 frequently.  Depending on how
  your tests go, you could also consider doing something like:
  |ap-psf| > max(3*uncertainty, 0.1) or something, to keep in low S/N
  points.  I'm not sure how valuable those y band points are going to be
  in general, though.

   */
  
  //use only detections with good flags, psf_qf_perfect>0.95,|ap_mag - psf_inst_mag| < max(4*uncertainty, 0.1)
  
  int errfactor = 4;
  /*
  int number_gband_badflag=0;
  int number_rband_badflag=0;
  int number_iband_badflag=0;
  int number_zband_badflag=0;
  int number_yband_badflag=0;*/
  
  number_gband_out_as_gband_psf_qf_perfect_less0p95=0;
  number_rband_out_as_rband_psf_qf_perfect_less0p95=0;
  number_iband_out_as_iband_psf_qf_perfect_less0p95=0;
  number_zband_out_as_zband_psf_qf_perfect_less0p95=0;
  number_yband_out_as_yband_psf_qf_perfect_less0p95=0;
  number_gband_out_as_gband_ap_mag=0;
  number_rband_out_as_rband_ap_mag=0;
  number_iband_out_as_iband_ap_mag=0;
  number_zband_out_as_zband_ap_mag=0;
  number_yband_out_as_yband_ap_mag=0;
  
  
   number_gband_outchi2=0;
   number_rband_outchi2=0;
   number_iband_outchi2=0;
   number_zband_outchi2=0;
   number_yband_outchi2=0;
  
  

  if(gband_orig!=0)
  {
    dim_gband_orig = gband_orig->size;

      
    for(int i=0;i<gband_orig->size;i++)
    {
      goodobs = goodObservationFlag( gband_flags_vec_orig[i]/*,g_badflags*/);
   
     /* if(goodobs!=0)
      {number_gband_badflag++;
	
      }
      */
      if(!(gsl_vector_get(gband_psf_qf_perfect_vec,i)>0.95))
      {number_gband_out_as_gband_psf_qf_perfect_less0p95++;}
      
      if(!(abs(     gsl_vector_get(gband_ap_mag_vec,i) - gsl_vector_get(gband_psf_inst_mag_vec,i) ) < max(errfactor*gsl_vector_get(gbanderr_orig,i), 0.1) ))
      {number_gband_out_as_gband_ap_mag++;}
      
      if( (goodobs==0) && (gsl_vector_get(gband_psf_qf_perfect_vec,i)>0.95) 
	  && (abs(     gsl_vector_get(gband_ap_mag_vec,i) - gsl_vector_get(gband_psf_inst_mag_vec,i) ) < max(errfactor*gsl_vector_get(gbanderr_orig,i), 0.1) )   )
      {
	
	dim_gband_goodobs++;
      }
   
      
    }
  }

  if(rband_orig!=0)
  {
     dim_rband_orig = rband_orig->size;

    for(int i=0;i<rband_orig->size;i++)
    {
      
      goodobs = goodObservationFlag(  rband_flags_vec_orig[i]/*, r_badflags*/);
      
      /*if(goodobs!=0)
      {number_rband_badflag++;
      }*/
      
      if(!(gsl_vector_get(rband_psf_qf_perfect_vec,i)>0.95))
      {number_rband_out_as_rband_psf_qf_perfect_less0p95++;}
      
      if(!(abs(     gsl_vector_get(rband_ap_mag_vec,i) - gsl_vector_get(rband_psf_inst_mag_vec,i) ) < max(errfactor*gsl_vector_get(rbanderr_orig,i), 0.1) ))
      {number_rband_out_as_rband_ap_mag++;}
      
      if( (goodobs==0) && (gsl_vector_get(rband_psf_qf_perfect_vec,i)>0.95) 
	  && (abs(     gsl_vector_get(rband_ap_mag_vec,i) - gsl_vector_get(rband_psf_inst_mag_vec,i) ) < max(errfactor*gsl_vector_get(rbanderr_orig,i), 0.1) )   )
      {
	dim_rband_goodobs++;
      }
      
    }
  }
  
  if(iband_orig!=0)
  {
     dim_iband_orig = iband_orig->size;

    for(int i=0;i<iband_orig->size;i++)
    {
      goodobs = goodObservationFlag(  iband_flags_vec_orig[i]/*, i_badflags*/);
      
      /*if(goodobs!=0)
      {number_iband_badflag++;
      }*/
      
      if(!(gsl_vector_get(iband_psf_qf_perfect_vec,i)>0.95))
      {number_iband_out_as_iband_psf_qf_perfect_less0p95++;}
      
      if(!(abs(     gsl_vector_get(iband_ap_mag_vec,i) - gsl_vector_get(iband_psf_inst_mag_vec,i) ) < max(errfactor*gsl_vector_get(ibanderr_orig,i), 0.1) ))
      {number_iband_out_as_iband_ap_mag++;}
   
      if( (goodobs==0) && (gsl_vector_get(iband_psf_qf_perfect_vec,i)>0.95) 
	 && (abs(     gsl_vector_get(iband_ap_mag_vec,i) - gsl_vector_get(iband_psf_inst_mag_vec,i) ) < max(errfactor*gsl_vector_get(ibanderr_orig,i), 0.1) )   )
      {
	dim_iband_goodobs++;
      }
      
    }
  }
  
  if(zband_orig!=0)
  {

   dim_zband_orig = zband_orig->size;
    for(int i=0;i<zband_orig->size;i++)
    {
      goodobs = goodObservationFlag(  zband_flags_vec_orig[i]/*, z_badflags*/);
      
      /*if(goodobs!=0)
      {number_zband_badflag++;
	
      }*/
      
      if(!(gsl_vector_get(zband_psf_qf_perfect_vec,i)>0.95))
      {number_zband_out_as_zband_psf_qf_perfect_less0p95++;}
      
      if(!(abs(     gsl_vector_get(zband_ap_mag_vec,i) - gsl_vector_get(zband_psf_inst_mag_vec,i) ) < max(errfactor*gsl_vector_get(zbanderr_orig,i), 0.1) ))
      {number_zband_out_as_zband_ap_mag++;}
      
      
      if( (goodobs==0) && (gsl_vector_get(zband_psf_qf_perfect_vec,i)>0.95) 
	  && (abs(     gsl_vector_get(zband_ap_mag_vec,i) - gsl_vector_get(zband_psf_inst_mag_vec,i) ) < max(errfactor*gsl_vector_get(zbanderr_orig,i), 0.1) )   )
      {
	dim_zband_goodobs++;
      }

    }
  }
  
  if(yband_orig!=0)
  {
     dim_yband_orig = yband_orig->size;

    for(int i=0;i<yband_orig->size;i++)
    {
      goodobs = goodObservationFlag(  yband_flags_vec_orig[i]/*, y_badflags*/);
      
      /*if(goodobs!=0)
      {number_yband_badflag++;
	
      }*/
      
      if(!(gsl_vector_get(yband_psf_qf_perfect_vec,i)>0.95))
      {number_yband_out_as_yband_psf_qf_perfect_less0p95++;}
      
      if(!(abs(     gsl_vector_get(yband_ap_mag_vec,i) - gsl_vector_get(yband_psf_inst_mag_vec,i) ) < max(errfactor*gsl_vector_get(ybanderr_orig,i), 0.1) ))
      {number_yband_out_as_yband_ap_mag++;}
      
      if( (goodobs==0) && (gsl_vector_get(yband_psf_qf_perfect_vec,i)>0.95) 
	  && (abs(     gsl_vector_get(yband_ap_mag_vec,i) - gsl_vector_get(yband_psf_inst_mag_vec,i) ) < max(errfactor*gsl_vector_get(ybanderr_orig,i), 0.1) )   )
      {
	dim_yband_goodobs++;
      }

    }
  }
  
  /*
  std::cout<<"dim_gband_orig "<<dim_gband_orig<<std::endl;
  std::cout<<"dim_rband_orig "<<dim_rband_orig<<std::endl;
  std::cout<<"dim_iband_orig "<<dim_iband_orig<<std::endl;
  std::cout<<"dim_zband_orig "<<dim_zband_orig<<std::endl;
  std::cout<<"dim_yband_orig "<<dim_yband_orig<<std::endl;
  std::cout<<"total before cleaning: "<<dim_gband_orig+dim_rband_orig+dim_iband_orig+dim_zband_orig+dim_yband_orig<<std::endl;
  
  std::cout<<"......................................"<<std::endl;
  std::cout<<"good obs. (flags and yband_psf_qf_perfect_vec, yband_ap_mag_vec,yband_psf_inst_mag_vec  )"<<std::endl;
  std::cout<<"dim_gband_goodobs "<<dim_gband_goodobs<<std::endl;
  std::cout<<"dim_rband_goodobs "<<dim_rband_goodobs<<std::endl;
  std::cout<<"dim_iband_goodobs "<<dim_iband_goodobs<<std::endl;
  std::cout<<"dim_zband_goodobs "<<dim_zband_goodobs<<std::endl;
  std::cout<<"dim_yband_goodobs "<<dim_yband_goodobs<<std::endl;
  */
  
  /*
 
  gband_badflags="";
  rband_badflags="";
  iband_badflags="";
  zband_badflags="";
  yband_badflags="";*/
  
  
  if(dim_gband_goodobs!=0)
  {
    int j=0;
    gband = gsl_vector_calloc(dim_gband_goodobs);
    gbanderr = gsl_vector_calloc(dim_gband_goodobs);
    gbandtime = gsl_vector_calloc(dim_gband_goodobs);
    
    for(int i=0;i<gband_orig->size;i++)
    {
      goodobs = goodObservationFlag(  gband_flags_vec_orig[i]/*, gband_badflags*/);
   
        if( (goodobs==0) && (gsl_vector_get(gband_psf_qf_perfect_vec,i)>0.95) 
	  && (abs(     gsl_vector_get(gband_ap_mag_vec,i) - gsl_vector_get(gband_psf_inst_mag_vec,i) ) < max(errfactor*gsl_vector_get(gbanderr_orig,i), 0.1) )   )
      {
	
	gsl_vector_set(gband,j,gsl_vector_get(gband_orig,i));
	gsl_vector_set(gbanderr,j,gsl_vector_get(gbanderr_orig,i));
	gsl_vector_set(gbandtime,j,gsl_vector_get(gbandtime_orig,i));
	
	j++;
      }
    
    }
  
  }

  if(dim_rband_goodobs!=0)
  {
    int j=0;
    rband = gsl_vector_calloc(dim_rband_goodobs);
    rbanderr = gsl_vector_calloc(dim_rband_goodobs);
    rbandtime = gsl_vector_calloc(dim_rband_goodobs);
    
    for(int i=0;i<rband_orig->size;i++)
    {
      goodobs = goodObservationFlag(  rband_flags_vec_orig[i]/*, rband_badflags*/);
      if( (goodobs==0) && (gsl_vector_get(rband_psf_qf_perfect_vec,i)>0.95) 
	  && (abs(     gsl_vector_get(rband_ap_mag_vec,i) - gsl_vector_get(rband_psf_inst_mag_vec,i) ) < max(errfactor*gsl_vector_get(rbanderr_orig,i), 0.1) )   )
      {
	gsl_vector_set(rband,j,gsl_vector_get(rband_orig,i));
	gsl_vector_set(rbanderr,j,gsl_vector_get(rbanderr_orig,i));
	gsl_vector_set(rbandtime,j,gsl_vector_get(rbandtime_orig,i));
	
	j++;
      }
    }
  }
  
  
  if(dim_iband_goodobs!=0)
  {
    int j=0;
    iband = gsl_vector_calloc(dim_iband_goodobs);
    ibanderr = gsl_vector_calloc(dim_iband_goodobs);
    ibandtime = gsl_vector_calloc(dim_iband_goodobs);
    
    for(int i=0;i<iband_orig->size;i++)
    {
      goodobs = goodObservationFlag(  iband_flags_vec_orig[i]/*, iband_badflags*/);
      if( (goodobs==0) && (gsl_vector_get(iband_psf_qf_perfect_vec,i)>0.95) 
	  && (abs(     gsl_vector_get(iband_ap_mag_vec,i) - gsl_vector_get(iband_psf_inst_mag_vec,i) ) < max(errfactor*gsl_vector_get(ibanderr_orig,i), 0.1) )   )
      {
	gsl_vector_set(iband,j,gsl_vector_get(iband_orig,i));
	gsl_vector_set(ibanderr,j,gsl_vector_get(ibanderr_orig,i));
	gsl_vector_set(ibandtime,j,gsl_vector_get(ibandtime_orig,i));
	
	j++;
      }
    }
  }
   
  if(dim_zband_goodobs!=0)
  {
    int j=0;
    zband = gsl_vector_calloc(dim_zband_goodobs);
    zbanderr = gsl_vector_calloc(dim_zband_goodobs);
    zbandtime = gsl_vector_calloc(dim_zband_goodobs);
    
    for(int i=0;i<zband_orig->size;i++)
    {
      goodobs = goodObservationFlag(  zband_flags_vec_orig[i]/*, zband_badflags*/);
      if( (goodobs==0) && (gsl_vector_get(zband_psf_qf_perfect_vec,i)>0.95) 
	  && (abs(     gsl_vector_get(zband_ap_mag_vec,i) - gsl_vector_get(zband_psf_inst_mag_vec,i) ) < max(errfactor*gsl_vector_get(zbanderr_orig,i), 0.1) )   )
      {
	gsl_vector_set(zband,j,gsl_vector_get(zband_orig,i));
	gsl_vector_set(zbanderr,j,gsl_vector_get(zbanderr_orig,i));
	gsl_vector_set(zbandtime,j,gsl_vector_get(zbandtime_orig,i));
	
	j++;
      }
    }
    
  }
   
  if(dim_yband_goodobs!=0)
  {
    int j=0;
    yband = gsl_vector_calloc(dim_yband_goodobs);
    ybanderr = gsl_vector_calloc(dim_yband_goodobs);
    ybandtime = gsl_vector_calloc(dim_yband_goodobs);
    
    for(int i=0;i<yband_orig->size;i++)
    {
      goodobs = goodObservationFlag(  yband_flags_vec_orig[i]/*, yband_badflags*/);
     
      if( (goodobs==0) && (gsl_vector_get(yband_psf_qf_perfect_vec,i)>0.95) 
	  && (abs(     gsl_vector_get(yband_ap_mag_vec,i) - gsl_vector_get(yband_psf_inst_mag_vec,i) ) < max(errfactor*gsl_vector_get(ybanderr_orig,i), 0.1) )   )
      {
	gsl_vector_set(yband,j,gsl_vector_get(yband_orig,i));
	gsl_vector_set(ybanderr,j,gsl_vector_get(ybanderr_orig,i));
	gsl_vector_set(ybandtime,j,gsl_vector_get(ybandtime_orig,i));
	
	j++;
      }
    }
  }
	  
	  if(only_by_flags==1)
	    //cleaning is now done
	  {
	    
	     
  int allbands_dim = dim_gband_goodobs+dim_rband_goodobs+dim_iband_goodobs+dim_zband_goodobs+dim_yband_goodobs; 
 
 
  if(allbands_dim>0)
  {
	   // copy
	   if(dim_gband_goodobs>0)
	   {
	    gband_cleaned = gsl_vector_calloc(dim_gband_goodobs);
	    gbanderr_cleaned = gsl_vector_calloc(dim_gband_goodobs);
	    gbandtime_cleaned = gsl_vector_calloc(dim_gband_goodobs);
	     
	    gsl_vector_memcpy(gband_cleaned, gband);
	    gsl_vector_memcpy(gbanderr_cleaned, gbanderr);
	    gsl_vector_memcpy(gbandtime_cleaned, gbandtime);
	    
	   }
	   
	   if(dim_rband_goodobs>0)
	   {
	    rband_cleaned = gsl_vector_calloc(dim_rband_goodobs);
	    rbanderr_cleaned = gsl_vector_calloc(dim_rband_goodobs);
	    rbandtime_cleaned = gsl_vector_calloc(dim_rband_goodobs);
	    
	    gsl_vector_memcpy(rband_cleaned, rband);
	    gsl_vector_memcpy(rbanderr_cleaned, rbanderr);
	    gsl_vector_memcpy(rbandtime_cleaned, rbandtime);
	   
	   }
	   if (dim_iband_goodobs>0)
	   {  
	    iband_cleaned = gsl_vector_calloc(dim_iband_goodobs);
	    ibanderr_cleaned = gsl_vector_calloc(dim_iband_goodobs);
	    ibandtime_cleaned = gsl_vector_calloc(dim_iband_goodobs);
	    
	    gsl_vector_memcpy(iband_cleaned, iband);
	    gsl_vector_memcpy(ibanderr_cleaned, ibanderr);
	    gsl_vector_memcpy(ibandtime_cleaned, ibandtime);
	    
	   }
	   if(dim_zband_goodobs>0)
	   {
	    zband_cleaned = gsl_vector_calloc(dim_zband_goodobs);
	    zbanderr_cleaned = gsl_vector_calloc(dim_zband_goodobs);
	    zbandtime_cleaned = gsl_vector_calloc(dim_zband_goodobs);
	    
	    gsl_vector_memcpy(zband_cleaned, zband);
	    gsl_vector_memcpy(zbanderr_cleaned, zbanderr);
	    gsl_vector_memcpy(zbandtime_cleaned, zbandtime);
	   
	   }
	   if(dim_yband_goodobs>0)
	   {
	    yband_cleaned = gsl_vector_calloc(dim_yband_goodobs);
	    ybanderr_cleaned = gsl_vector_calloc(dim_yband_goodobs);
	    ybandtime_cleaned = gsl_vector_calloc(dim_yband_goodobs);
	
	    gsl_vector_memcpy(yband_cleaned, yband);
	    gsl_vector_memcpy(ybanderr_cleaned, ybanderr);
	    gsl_vector_memcpy(ybandtime_cleaned, ybandtime);
	 
	   }
	  
	   
	 
	
	// clean up:
	    
	    
  int status;
  if( (gband_cleaned==0) && (rband_cleaned==0) && (iband_cleaned==0) && (zband_cleaned==0) && (yband_cleaned==0))
  {
    status= -1;
  }
  else
  {status = 0;
  }
  
		
      if( gband!=0)
      {
	gsl_vector_free(gband);
      }
      
      if( gbanderr!=0)
      {
	gsl_vector_free(gbanderr);
      }

      if( gbandtime!=0)
      {
	gsl_vector_free(gbandtime);
      }
      
	
      if( rband!=0)
      {
	gsl_vector_free(rband);
      }
      
      if( rbanderr!=0)
      {
	gsl_vector_free(rbanderr);
      }

      if( rbandtime!=0)
      {
	gsl_vector_free(rbandtime);
      }
      
	
      if( iband!=0)
      {
	gsl_vector_free(iband);
      }
      
      if( ibanderr!=0)
      {
	gsl_vector_free(ibanderr);
      }

      if( ibandtime!=0)
      {
	gsl_vector_free(ibandtime);
      }
      
	
      if( zband!=0)
      {
	gsl_vector_free(zband);
      }
      
      if( zbanderr!=0)
      {
	gsl_vector_free(zbanderr);
      }

      if( zbandtime!=0)
      {
	gsl_vector_free(zbandtime);
      }
	    
		
      if( yband!=0)
      {
	gsl_vector_free(yband);
      }
      
      if( ybanderr!=0)
      {
	gsl_vector_free(ybanderr);
      }

      if( ybandtime!=0)
      {
	gsl_vector_free(ybandtime);
      }


  
  
 return status;
  }
  
  else
  {
    
    gband_cleaned=0;
    gbanderr_cleaned=0;
    gbandtime_cleaned=0;
    
    rband_cleaned=0;
    rbanderr_cleaned=0;
    rbandtime_cleaned=0;
  
    iband_cleaned=0;
    ibanderr_cleaned=0;
    ibandtime_cleaned=0;
    
    zband_cleaned=0;
    zbanderr_cleaned=0;
    zbandtime_cleaned=0;
    
    yband_cleaned=0;
    ybanderr_cleaned=0;
    ybandtime_cleaned=0;
  
    return -1;
    
  }
	
	  }
	 else{

  int dim_gband=0;
  if(gband!=0)
  {
    dim_gband = gband->size;
    gband_cleaned_tmp = gsl_vector_calloc(dim_gband);   
    gbanderr_cleaned_tmp = gsl_vector_calloc(dim_gband);   
    gbandtime_cleaned_tmp = gsl_vector_calloc(dim_gband);   
  }
  int dim_rband=0;
  if(rband!=0)
  {
    dim_rband = rband->size;
    rband_cleaned_tmp = gsl_vector_calloc(dim_rband);   
    rbanderr_cleaned_tmp = gsl_vector_calloc(dim_rband);   
    rbandtime_cleaned_tmp = gsl_vector_calloc(dim_rband);   
  }
  int dim_iband=0;
  if(iband!=0)
  {
    dim_iband = iband->size;
    iband_cleaned_tmp = gsl_vector_calloc(dim_iband);   
    ibanderr_cleaned_tmp = gsl_vector_calloc(dim_iband);   
    ibandtime_cleaned_tmp = gsl_vector_calloc(dim_iband);   
  }
  int dim_zband=0;
  if(zband!=0)
  {
    dim_zband = zband->size;
    zband_cleaned_tmp = gsl_vector_calloc(dim_zband);   
    zbanderr_cleaned_tmp = gsl_vector_calloc(dim_zband);
    zbandtime_cleaned_tmp = gsl_vector_calloc(dim_zband);
  }
  int dim_yband=0;
  if(yband!=0)
  {
    dim_yband = yband->size;
    yband_cleaned_tmp = gsl_vector_calloc(dim_yband);   
    ybanderr_cleaned_tmp = gsl_vector_calloc(dim_yband);   
    ybandtime_cleaned_tmp = gsl_vector_calloc(dim_yband);   
  }
  
  
  total_beforechi2=dim_gband+dim_rband+dim_iband+dim_zband+dim_yband;
 
	  //ib) compute chi terms

  
	std::vector<double> chiterms_g_vec;
	std::vector<double> chiterms_r_vec;
        std::vector<double> chiterms_i_vec;
	std::vector<double> chiterms_z_vec;
	std::vector<double> chiterms_y_vec;
    
    
	
	std::vector<double> chiterms_vec;
	
	
	
    double value;
	    
    
    double meanmag_g=0;
    double meanmag_r=0;
    double meanmag_i=0;
    double meanmag_z=0;
    double meanmag_y=0;
    
    
    double sum1=0;
    double sum2=0;
    double factor1=0;
    double factor2=0;
    
     if(dim_gband!=0)
	  {
	  
	    sum1=0;
	    sum2=0;
	    factor1=0;
	    
	    for(int i=0;i<gbanderr->size;i++)
	    {
	      sum1 += 1/ (pow(gsl_vector_get(gbanderr,i),2));
	    }
	    factor1=1/sum1;
	    
	    
	    for(int i=0;i<gbanderr->size;i++)
	    {
	      sum2 += gsl_vector_get(gband,i) / (pow(gsl_vector_get(gbanderr,i),2));
	    }
	    
	    meanmag_g=factor1*sum2;
	   
	    
	  }
    
    
     if(dim_rband!=0)
	  {
	  
	    sum1=0;
	    sum2=0;
	    factor1=0;
	    
	    for(int i=0;i<rbanderr->size;i++)
	    {
	      sum1 += 1/ (pow(gsl_vector_get(rbanderr,i),2));
	    }
	    factor1=1/sum1;
	    
	    
	    for(int i=0;i<rbanderr->size;i++)
	    {
	      sum2 += gsl_vector_get(rband,i) / (pow(gsl_vector_get(rbanderr,i),2));
	    }
	    
	    meanmag_r=factor1*sum2;
	   
	    
	  }
    
    
    
     if(dim_iband!=0)
	  {
	  
	    sum1=0;
	    sum2=0;
	    factor1=0;
	    
	    for(int i=0;i<ibanderr->size;i++)
	    {
	      sum1 += 1/ (pow(gsl_vector_get(ibanderr,i),2));
	    }
	    factor1=1/sum1;
	    
	    
	    for(int i=0;i<ibanderr->size;i++)
	    {
	      sum2 += gsl_vector_get(iband,i) / (pow(gsl_vector_get(ibanderr,i),2));
	    }
	    
	    meanmag_i=factor1*sum2;
	   
	    
	  }
    
    
     if(dim_zband!=0)
	  {
	  
	    sum1=0;
	    sum2=0;
	    factor1=0;
	    
	    for(int i=0;i<zbanderr->size;i++)
	    {
	      sum1 += 1/ (pow(gsl_vector_get(zbanderr,i),2));
	    }
	    factor1=1/sum1;
	    
	    
	    for(int i=0;i<zbanderr->size;i++)
	    {
	      sum2 += gsl_vector_get(zband,i) / (pow(gsl_vector_get(zbanderr,i),2));
	    }
	    
	    meanmag_z=factor1*sum2;
	   
	    
	  }
    
    
    
    
    
    
     if(dim_yband!=0)
	  {
	  
	    sum1=0;
	    sum2=0;
	    factor1=0;
	    
	    for(int i=0;i<ybanderr->size;i++)
	    {
	      sum1 += 1/ (pow(gsl_vector_get(ybanderr,i),2));
	    }
	    factor1=1/sum1;
	    
	    
	    for(int i=0;i<ybanderr->size;i++)
	    {
	      sum2 += gsl_vector_get(yband,i) / (pow(gsl_vector_get(ybanderr,i),2));
	    }
	    
	    meanmag_y=factor1*sum2;
	   
	    
	  }
    
    
   
    if(dim_gband!=0)
    {
      for(int i=0;i<gband->size;i++)
      {
	value =(gsl_vector_get(gband,i) - meanmag_g) / sqrt( pow(/*1.3**/gsl_vector_get(gbanderr,i),2)+pow(0.015, 2));
	chiterms_vec.push_back(value);
      }
    }
      
  if(dim_rband!=0)
    {
      for(int i=0;i<rband->size;i++)
      {
	value =(gsl_vector_get(rband,i) - meanmag_r) / sqrt( pow(/*1.3**/gsl_vector_get(rbanderr,i),2)+pow(0.015, 2));
	chiterms_vec.push_back(value);
      }
    }
    
    if(dim_iband!=0)
    {
      for(int i=0;i<iband->size;i++)
      {
	value =(gsl_vector_get(iband,i) - meanmag_i) / sqrt( pow(/*1.3**/gsl_vector_get(ibanderr,i),2)+pow(0.015, 2));
	chiterms_vec.push_back(value);
      }
    }
    
    if(dim_zband!=0)
    {
      for(int i=0;i<zband->size;i++)
      {
	value =(gsl_vector_get(zband,i) - meanmag_z) / sqrt( pow(/*1.3**/gsl_vector_get(zbanderr,i),2)+pow(0.015, 2));
	chiterms_vec.push_back(value);
      }
    }
    
    if(dim_yband!=0)
    {
      for(int i=0;i<yband->size;i++)
      {
	value =(gsl_vector_get(yband,i) - meanmag_y) / sqrt( pow(/*1.3**/gsl_vector_get(ybanderr,i),2)+pow(0.015, 2));
	chiterms_vec.push_back(value);
      }
    }
    
    
  
  double median_chiterms;
  
  double IQR_chiterms;
  
  // ii)  compute the median of the $\chi$s
  // iii) compute the interquartile range, this is a number computed by IQR = Q3 −  Q1, where Q3 and Q1 are the upper and lower quartiles
  

    double* chiterms_tmp=0;
    
    
  int allbands_dim = dim_gband+dim_rband+dim_iband+dim_zband+dim_yband;
   
    if(allbands_dim !=0)
    {
       chiterms_tmp=new double [allbands_dim];
       
       
      for(int i=0;i<allbands_dim;i++)
      {
	chiterms_tmp[i]=chiterms_vec[i];
      }
      
      gsl_sort (chiterms_tmp, 1, allbands_dim);
      median_chiterms = gsl_stats_median_from_sorted_data (chiterms_tmp, 1, allbands_dim);
      IQR_chiterms = gsl_stats_quantile_from_sorted_data (chiterms_tmp, 1, allbands_dim, 0.75) 
		      - gsl_stats_quantile_from_sorted_data (chiterms_tmp, 1,allbands_dim, 0.25) ;

    }
    
    
    

  /*iv)  throw away the points having a $| \chi - \chi_{\mathrm{\median}}| > |dist \times IQR \times 1.349/2 | $, 
   * at maximum 10 \% of all light curve points*/
  
   int dist = 5;//2;
   
   

 /*
    the 10-percent-rule makes it a bit more complicated, so 
    we need temporary arrays for sorting: 
    make one vector of concatenated chiterms_..._vecs, 
    make one vector of numbers 0 to dim_gband+dim_rband+dim_iband+dim_zband+dim_yband (this is p)
 */
    
  

  
  if(allbands_dim>0)
  {
    
    
  gsl_vector *chiterms_minus_max_chi_allowed_allbands_vec = gsl_vector_calloc(allbands_dim); 

   //copy into chiterms_allbands_vec .....
   
  for(int i=0;i<allbands_dim;i++)
    {      
      gsl_vector_set(chiterms_minus_max_chi_allowed_allbands_vec,i, 
		     /*abs*/(abs(chiterms_vec[i] - median_chiterms)-abs(dist*IQR_chiterms * 1.349/2)) );
      
    } 
    
    

  
  
  
  //sorted by ascending order
  //p=new int [allbands_dim];
  //gsl_sort_index((size_t*)p, chiterms_minus_max_chi_allowed_allbands_vec,1, allbands_dim);
  //gsl_permute((size_t*)p, chiterms_minus_max_chi_allowed_allbands_vec,1, allbands_dim);
  
  
  
    gsl_permutation * p = gsl_permutation_alloc (allbands_dim);
    gsl_permutation_init (p);
  
  gsl_sort_vector_index (p,chiterms_minus_max_chi_allowed_allbands_vec);
 
  gsl_permute_vector(p,chiterms_minus_max_chi_allowed_allbands_vec);
 
  //maximum number of points to throw away
  int max_away_points = allbands_dim / 10;
 
  //and now throw away as making cleaned vectors _tmp
  
  int gband_cleaned_idx=0;
  int rband_cleaned_idx=0;
  int iband_cleaned_idx=0;
  int zband_cleaned_idx=0;
  int yband_cleaned_idx=0;
  
  int counter=0;
  
  for(int i=0;i<allbands_dim;i++)
  {
      //sorted by ascending order, so:
      //if chiterms_minus_max_chi_allowed_allbands_vec entry > 0 AND i > allbands_dim - max_away_points - 1: do not copy into cleaned vectors
      if(!(      (gsl_vector_get(chiterms_minus_max_chi_allowed_allbands_vec,i) > 0) && (i > allbands_dim - max_away_points - 1)) )
      {
counter++;
	//copy into _cleaned_tmp
	int index = gsl_permutation_get(p,i);

	//std::cout<<"index "<<index<<std::endl;
	if(index < dim_gband)
	{
	  gsl_vector_set(gband_cleaned_tmp,gband_cleaned_idx,gsl_vector_get(gband,index));
	  gsl_vector_set(gbanderr_cleaned_tmp,gband_cleaned_idx,gsl_vector_get(gbanderr,index));
	  gsl_vector_set(gbandtime_cleaned_tmp,gband_cleaned_idx,gsl_vector_get(gbandtime,index));
	  gband_cleaned_idx++;
	//  std::cout<<"copy to g"<<std::endl;
	}
	
	if( (index < dim_gband+dim_rband) && (index >= dim_gband) )
	{
	  gsl_vector_set(rband_cleaned_tmp,rband_cleaned_idx,gsl_vector_get(rband,index-dim_gband));
	  gsl_vector_set(rbanderr_cleaned_tmp,rband_cleaned_idx,gsl_vector_get(rbanderr,index-dim_gband));
	  gsl_vector_set(rbandtime_cleaned_tmp,rband_cleaned_idx,gsl_vector_get(rbandtime,index-dim_gband));

	  rband_cleaned_idx++;
	}
	
	if( (index < dim_gband+dim_rband+dim_iband) && (index >= dim_gband+dim_rband) )
	{
	  gsl_vector_set(iband_cleaned_tmp,iband_cleaned_idx,gsl_vector_get(iband,index-dim_gband-dim_rband));
	  gsl_vector_set(ibanderr_cleaned_tmp,iband_cleaned_idx,gsl_vector_get(ibanderr,index-dim_gband-dim_rband));
	  gsl_vector_set(ibandtime_cleaned_tmp,iband_cleaned_idx,gsl_vector_get(ibandtime,index-dim_gband-dim_rband));

	  iband_cleaned_idx++;
	}
	
	if( (index < dim_gband+dim_rband+dim_iband+dim_zband) && (index >= dim_gband+dim_rband+dim_iband) )
	{
	  gsl_vector_set(zband_cleaned_tmp,zband_cleaned_idx,gsl_vector_get(zband,index-dim_gband-dim_rband-dim_iband));
	  gsl_vector_set(zbanderr_cleaned_tmp,zband_cleaned_idx,gsl_vector_get(zbanderr,index-dim_gband-dim_rband-dim_iband));
	  gsl_vector_set(zbandtime_cleaned_tmp,zband_cleaned_idx,gsl_vector_get(zbandtime,index-dim_gband-dim_rband-dim_iband));

	  zband_cleaned_idx++;
	}
	
	if( (index < dim_gband+dim_rband+dim_iband+dim_zband+dim_yband) && (index >= dim_gband+dim_rband+dim_iband+dim_zband) )
	{
	  gsl_vector_set(yband_cleaned_tmp,yband_cleaned_idx,gsl_vector_get(yband,index-dim_gband-dim_rband-dim_iband-dim_zband));
	  gsl_vector_set(ybanderr_cleaned_tmp,yband_cleaned_idx,gsl_vector_get(ybanderr,index-dim_gband-dim_rband-dim_iband-dim_zband));
	  gsl_vector_set(ybandtime_cleaned_tmp,yband_cleaned_idx,gsl_vector_get(ybandtime,index-dim_gband-dim_rband-dim_iband-dim_zband));

	  yband_cleaned_idx++;
	}
	
      }
      
      
  }
  
  
  gsl_permutation * p_g;
  gsl_permutation * p_r;
  gsl_permutation * p_i;
  gsl_permutation * p_z;
  gsl_permutation * p_y;
  
  //sort vectors
  
  if(dim_gband>0)
  {
    p_g = gsl_permutation_alloc (dim_gband);
    gsl_permutation_init (p_g);
  
    gsl_sort_vector_index (p_g,gbandtime_cleaned_tmp);
 
    gsl_permute_vector(p_g,gbandtime_cleaned_tmp);
    gsl_permute_vector(p_g,gband_cleaned_tmp);
    gsl_permute_vector(p_g,gbanderr_cleaned_tmp);
  }
  
  if(dim_rband>0)
  {
    p_r = gsl_permutation_alloc (dim_rband);
    gsl_permutation_init (p_r);
  
    gsl_sort_vector_index (p_r,rbandtime_cleaned_tmp);
 
    gsl_permute_vector(p_r,rbandtime_cleaned_tmp);
    gsl_permute_vector(p_r,rband_cleaned_tmp);
    gsl_permute_vector(p_r,rbanderr_cleaned_tmp);
  }
  
  if(dim_iband>0)
  {
    p_i = gsl_permutation_alloc (dim_iband);
    gsl_permutation_init (p_i);
  
    gsl_sort_vector_index (p_i,ibandtime_cleaned_tmp);
 
    gsl_permute_vector(p_i,ibandtime_cleaned_tmp);
    gsl_permute_vector(p_i,iband_cleaned_tmp);
    gsl_permute_vector(p_i,ibanderr_cleaned_tmp);
  }
  
  if(dim_zband>0)
  {
    p_z = gsl_permutation_alloc (dim_zband);
    gsl_permutation_init (p_z);
  
    gsl_sort_vector_index (p_z,zbandtime_cleaned_tmp);
 
    gsl_permute_vector(p_z,zbandtime_cleaned_tmp);
    gsl_permute_vector(p_z,zband_cleaned_tmp);
    gsl_permute_vector(p_z,zbanderr_cleaned_tmp);
  }
  
  if(dim_yband>0)
  {
    p_y = gsl_permutation_alloc (dim_yband);
    gsl_permutation_init (p_y);
    
    gsl_sort_vector_index (p_y,ybandtime_cleaned_tmp);
  
    gsl_permute_vector(p_y,ybandtime_cleaned_tmp);
    gsl_permute_vector(p_y,yband_cleaned_tmp);
    gsl_permute_vector(p_y,ybanderr_cleaned_tmp);
  }
 
  //make cleaned vectors 
  
  int dim_gband_cleaned=0;
  int dim_rband_cleaned=0;
  int dim_iband_cleaned=0;
  int dim_zband_cleaned=0;
  int dim_yband_cleaned=0;
  if(dim_gband>0)
  {
    for(int i=0;i<gband_cleaned_tmp->size;i++)
    {
      if(gsl_vector_get(gband_cleaned_tmp,i)!=0)
      {
	dim_gband_cleaned++;
      }
    }
  }
  
  if(dim_rband>0)
  {
    for(int i=0;i<rband_cleaned_tmp->size;i++)
    {
      if(gsl_vector_get(rband_cleaned_tmp,i)!=0)
      {
	dim_rband_cleaned++;
      }
    }
  }
  
  if(dim_iband>0)
  {
    for(int i=0;i<iband_cleaned_tmp->size;i++)
    {
      if(gsl_vector_get(iband_cleaned_tmp,i)!=0)
      {
	dim_iband_cleaned++;
      }
    }
  }
  
  if(dim_zband>0)
  {
    for(int i=0;i<zband_cleaned_tmp->size;i++)
    {
      if(gsl_vector_get(zband_cleaned_tmp,i)!=0)
      {
	dim_zband_cleaned++;
      }
    }
  }
    
  if(dim_yband>0)
  {
    for(int i=0;i<yband_cleaned_tmp->size;i++)
    {
    
      if(gsl_vector_get(yband_cleaned_tmp,i)!=0)
      {
	dim_yband_cleaned++;
      }
    }
  }
  
  //calloc
  if(dim_gband_cleaned > 0)
  {
    gband_cleaned = gsl_vector_calloc (dim_gband_cleaned);
    gbanderr_cleaned = gsl_vector_calloc (dim_gband_cleaned);
    gbandtime_cleaned = gsl_vector_calloc (dim_gband_cleaned);
  }
  
  if(dim_rband_cleaned > 0)
  {
    rband_cleaned = gsl_vector_calloc (dim_rband_cleaned);
    rbanderr_cleaned = gsl_vector_calloc (dim_rband_cleaned);
    rbandtime_cleaned = gsl_vector_calloc (dim_rband_cleaned);
  }
  
  if(dim_iband_cleaned > 0)
  {
    iband_cleaned = gsl_vector_calloc (dim_iband_cleaned);
    ibanderr_cleaned = gsl_vector_calloc (dim_iband_cleaned);
    ibandtime_cleaned = gsl_vector_calloc (dim_iband_cleaned);
  }
  
  if(dim_zband_cleaned > 0)
  {
    zband_cleaned = gsl_vector_calloc (dim_zband_cleaned);
    zbanderr_cleaned = gsl_vector_calloc (dim_zband_cleaned);
    zbandtime_cleaned = gsl_vector_calloc (dim_zband_cleaned);
  }
  
  if(dim_yband_cleaned > 0)
  {
    yband_cleaned = gsl_vector_calloc (dim_yband_cleaned);
    ybanderr_cleaned = gsl_vector_calloc (dim_yband_cleaned);
    ybandtime_cleaned = gsl_vector_calloc (dim_yband_cleaned);
  }
  
  /*
 std::cout<<"dim_gband_cleaned "<<dim_gband_cleaned<<std::endl; 
 std::cout<<"dim_rband_cleaned "<<dim_rband_cleaned<<std::endl;
 std::cout<<"dim_iband_cleaned "<<dim_iband_cleaned<<std::endl;
 std::cout<<"dim_zband_cleaned "<<dim_zband_cleaned<<std::endl;
 std::cout<<"dim_yband_cleaned "<<dim_yband_cleaned<<std::endl;
  */
  
  int j=0;
 
  
  if(dim_gband_cleaned > 0)
  {
    for(int i=0;i<gband_cleaned_tmp->size;i++)
    {
      if(gsl_vector_get(gband_cleaned_tmp,i)!=0)
      {
	gsl_vector_set(gband_cleaned, j, gsl_vector_get(gband_cleaned_tmp,i));
	gsl_vector_set(gbanderr_cleaned, j, gsl_vector_get(gbanderr_cleaned_tmp,i));
	gsl_vector_set(gbandtime_cleaned, j, gsl_vector_get(gbandtime_cleaned_tmp,i));
	  
	j++;
      }
    }
  }
  j=0;
 
  if(dim_rband_cleaned > 0)
  {
    for(int i=0;i<rband_cleaned_tmp->size;i++)
    {
      if(gsl_vector_get(rband_cleaned_tmp,i)!=0)
      {
	gsl_vector_set(rband_cleaned, j, gsl_vector_get(rband_cleaned_tmp,i));
	gsl_vector_set(rbanderr_cleaned, j, gsl_vector_get(rbanderr_cleaned_tmp,i));
	gsl_vector_set(rbandtime_cleaned, j, gsl_vector_get(rbandtime_cleaned_tmp,i));
	  
	j++;
      }
    }
  
  }
  
  j=0;
  if(dim_iband_cleaned > 0)
  {
    for(int i=0;i<iband_cleaned_tmp->size;i++)
    {
      if(gsl_vector_get(iband_cleaned_tmp,i)!=0)
      {
	gsl_vector_set(iband_cleaned, j, gsl_vector_get(iband_cleaned_tmp,i));
	gsl_vector_set(ibanderr_cleaned, j, gsl_vector_get(ibanderr_cleaned_tmp,i));
	gsl_vector_set(ibandtime_cleaned, j, gsl_vector_get(ibandtime_cleaned_tmp,i));
	  
	j++;
      }
    }
  }
  j=0;
  
  if(dim_zband_cleaned > 0)
  {
    for(int i=0;i<zband_cleaned_tmp->size;i++)
    {
      if(gsl_vector_get(zband_cleaned_tmp,i)!=0)
      {
	gsl_vector_set(zband_cleaned, j, gsl_vector_get(zband_cleaned_tmp,i));
	gsl_vector_set(zbanderr_cleaned, j, gsl_vector_get(zbanderr_cleaned_tmp,i));
	gsl_vector_set(zbandtime_cleaned, j, gsl_vector_get(zbandtime_cleaned_tmp,i));
	  
	j++;
      }
    }
  }
  
  
  j=0;
  if(dim_yband_cleaned > 0)
  {
      
    for(int i=0;i<yband_cleaned_tmp->size;i++)
    {
      if(gsl_vector_get(yband_cleaned_tmp,i)!=0)
      {
	gsl_vector_set(yband_cleaned, j, gsl_vector_get(yband_cleaned_tmp,i));
	gsl_vector_set(ybanderr_cleaned, j, gsl_vector_get(ybanderr_cleaned_tmp,i));
	gsl_vector_set(ybandtime_cleaned, j, gsl_vector_get(ybandtime_cleaned_tmp,i));
	  
	j++;
      }
    }
  }
  
  number_gband_outchi2=dim_gband-dim_gband_cleaned;
  number_rband_outchi2=dim_rband-dim_rband_cleaned;
  number_iband_outchi2=dim_iband-dim_iband_cleaned;
  number_zband_outchi2=dim_zband-dim_zband_cleaned;
  number_yband_outchi2=dim_yband-dim_yband_cleaned;
  
  
  
  //free
  
  if( gband_cleaned_tmp!=0)
  {
    gsl_vector_free(gband_cleaned_tmp);
  }
  if( gbanderr_cleaned_tmp!=0)
  {
    gsl_vector_free(gbanderr_cleaned_tmp);
  }
  if( gbandtime_cleaned_tmp!=0)
  {
    gsl_vector_free(gbandtime_cleaned_tmp);
  }
  if( rband_cleaned_tmp!=0)
  {
    gsl_vector_free(rband_cleaned_tmp);
  }
  if( rbanderr_cleaned_tmp!=0)
  {
    gsl_vector_free(rbanderr_cleaned_tmp);
  }
  if( rbandtime_cleaned_tmp!=0)
  {
    gsl_vector_free(rbandtime_cleaned_tmp);
  }
  if( iband_cleaned_tmp!=0)
  {
    gsl_vector_free(iband_cleaned_tmp);
  }
  if( ibanderr_cleaned_tmp!=0)
  {
    gsl_vector_free(ibanderr_cleaned_tmp);
  }
  if( ibandtime_cleaned_tmp!=0)
  {
    gsl_vector_free(ibandtime_cleaned_tmp);
  }
  if( zband_cleaned_tmp!=0)
  {
    gsl_vector_free(zband_cleaned_tmp);
  }
  if( zbanderr_cleaned_tmp!=0)
  {
    gsl_vector_free(zbanderr_cleaned_tmp);
  }
  if( zbandtime_cleaned_tmp!=0)
  {
    gsl_vector_free(zbandtime_cleaned_tmp);
  }
  
  if( yband_cleaned_tmp!=0)
  {
    gsl_vector_free(yband_cleaned_tmp);
  }
  
  if( ybanderr_cleaned_tmp!=0)
  {
    gsl_vector_free(ybanderr_cleaned_tmp);
  }
  
  if( ybandtime_cleaned_tmp!=0)
  {
    gsl_vector_free(ybandtime_cleaned_tmp);
  }
  
  gsl_permutation_free(p);    //gsl permutation
  
  
  
  return 0;
  }
  
  else
  {
    
    gband_cleaned=0;
    gbanderr_cleaned=0;
    gbandtime_cleaned=0;
    
    rband_cleaned=0;
    rbanderr_cleaned=0;
    rbandtime_cleaned=0;
  
    iband_cleaned=0;
    ibanderr_cleaned=0;
    ibandtime_cleaned=0;
    
    zband_cleaned=0;
    zbanderr_cleaned=0;
    zbandtime_cleaned=0;
    
    yband_cleaned=0;
    ybanderr_cleaned=0;
    ybandtime_cleaned=0;
    return -1;
    
  }
  
}  
  
  
    
}
	  

	  
	  
	  
double loglikelihood_strucparam_multicolor(gsl_vector *contband, gsl_vector *bandindex_vec, gsl_matrix *m_mat, 
					   gsl_vector *meanmag_per_band_vec,
					   gsl_vector *time, gsl_vector *contbanderr, gsl_vector *strucfunc_params, 
					   gsl_vector* wavelengths_vec,
					   int fitmodel, gsl_vector *missingbands )
{

  int dim =contband->size;
  int thebandindex=0;
  gsl_matrix *cmatrix = gsl_matrix_alloc(dim,dim);//rows, columns
 
  gsl_matrix *cinvmatrix = gsl_matrix_alloc(dim,dim);//rows, columns
  double logdet_c;

  int used_number_of_bands = m_mat->size2;
 // std::cout<<"2 used_number_of_bands "<<used_number_of_bands<<std::endl;

  int number_of_bands = 5;
  double loglikelihood=0;
  double chisquare=0;
  double logdet_Vmu=0;
  
  double a_r = gsl_vector_get(strucfunc_params, 0);
  double tau = gsl_vector_get(strucfunc_params, 1);
  double alpha = gsl_vector_get(strucfunc_params, 2);
   
  double A=0;
 
  double wavelength_base = gsl_vector_get(wavelengths_vec, 1); //r band
  gsl_vector *a_vec = gsl_vector_alloc(number_of_bands);
 
  double val;
  for(int i=0;i<a_vec->size;i++)
  {
    val = a_r * pow(gsl_vector_get(wavelengths_vec, i)/wavelength_base, alpha);
    gsl_vector_set(a_vec, i, val);
    
  //  std::cout<<val<<std::endl;
  }
  //std::cout<<"---"<<std::endl;

 double baseline = 5;
  if(fitmodel==1) //power law model
  {
 
    for (int i = 0; i < dim; i++)
    {    
      for (int j = 0; j < dim; j++)
      {   
	A=sqrt(gsl_vector_get(a_vec, gsl_vector_get(bandindex_vec,i))*gsl_vector_get(a_vec, gsl_vector_get(bandindex_vec,j)));

	gsl_matrix_set (cmatrix, i, j, pow(A,2)*(pow(baseline,tau) - 0.5*pow(abs(gsl_vector_get(time,i)-gsl_vector_get(time,j))/365,tau))  );  //power law
      }
      A=sqrt(gsl_vector_get(a_vec, gsl_vector_get(bandindex_vec,i))*gsl_vector_get(a_vec, gsl_vector_get(bandindex_vec,i)));
     gsl_matrix_set(cmatrix,i,i,pow(A,2)*pow(baseline,tau)+ pow(gsl_vector_get(contbanderr,i),2) );
    }
  }
  
  
  if(fitmodel==2) //drw
  {

   
    for (int i = 0; i < dim; i++)
    {    
      for (int j = 0; j < dim; j++)
      {   
	A=sqrt(gsl_vector_get(a_vec, gsl_vector_get(bandindex_vec,i))*gsl_vector_get(a_vec, gsl_vector_get(bandindex_vec,j)));
 
	
	 gsl_matrix_set (cmatrix, i, j,  A*A*exp( -(abs(gsl_vector_get(time,i)-gsl_vector_get(time,j)))/tau)  ); //drw
     
    }
   
     
    A=sqrt(gsl_vector_get(a_vec, gsl_vector_get(bandindex_vec,i))*gsl_vector_get(a_vec, gsl_vector_get(bandindex_vec,i)));
    gsl_matrix_set (cmatrix, i, i, A*A+ pow(gsl_vector_get(contbanderr,i),2) );  
    
    
  
    }
  }
 
   
  int nlin = used_number_of_bands;
  int npt = dim;
  double logdet_tlmat2=0;
   
  double det_cpn, det_tlmat2;  
   
  gsl_matrix *ctempmatrix = gsl_matrix_alloc(dim,dim);//rows, columns
 
  
  gsl_matrix *cinvperpmatrix = gsl_matrix_alloc(dim,dim);//rows, columns

  gsl_matrix *tlmat1 = gsl_matrix_alloc(nlin, npt);    
  gsl_matrix *tlmat3 = gsl_matrix_alloc(nlin, npt);
  gsl_matrix *tlmat5 = gsl_matrix_alloc(npt,nlin);
  gsl_matrix *lmat = gsl_matrix_alloc(npt,nlin);
  gsl_matrix *tlmat2 = gsl_matrix_alloc(nlin, nlin);
  gsl_matrix *tlmat2inv = gsl_matrix_alloc(nlin, nlin);
   
  gsl_vector *linpar = gsl_vector_alloc(nlin);
   
  gsl_vector *lpar = gsl_vector_alloc(nlin);
 
  gsl_matrix_memcpy (lmat, m_mat);
   

  
  //here: set lmat slightly different for "removing trend"
  
  
  /*
   up to now: lmat contains 1 in the column of corresponding band
   */
  
  
  
   //calculate the inverse of cmatrix 

  
 if (singular (cmatrix)) 
 {
   //std::cout<<"cmatrix singular"<<std::endl;
   loglikelihood=-1e32;
 }
 else{
   
  gsl_matrix_memcpy (cinvmatrix, cmatrix);

  
  inv(cinvmatrix);
 //std::cout<<"inverted cinvmatrix"<<std::endl;
  logdet_c = logdet(cmatrix);
  
  
  //calculate tlmat1 = lmat^T cinvmatrix 

  gsl_blas_dgemm (CblasTrans, CblasNoTrans,
                       1.0, m_mat, cinvmatrix,
                       0.0, tlmat1);
  
   

  gsl_blas_dgemv(CblasNoTrans,1.0,tlmat1,contband,0.0,linpar);    
 
  /*            tlmat2 =        L^T invC  L        = tlmat1 L
  calculate tlmat2 = lmag^T cinvmatrix lmat = tlmat1 lmat*/

  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                       1.0, tlmat1, m_mat,
                       0.0, tlmat2);
 
  
 if (singular (tlmat2)) 
 {
   //std::cout<<"tlmat2 singular"<<std::endl;
   loglikelihood=-1e32;
 }
else{
  
  
  logdet_tlmat2 =  abs( logdet(tlmat2)   ); 

  //calculate the inverse of tlmat2
      
  gsl_matrix_memcpy (tlmat2inv, tlmat2);
  inv(tlmat2inv);
 // std::cout<<"inverted tlmat2inv"<<std::endl;
   

  /*
  solve the equation to get lpar q:
              q = (L^T invC L)^(-1) L^T invC y  
  */

  gsl_blas_dgemv(CblasNoTrans,1.0,tlmat2inv,linpar,0.0,lpar);
 
  
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,
                       1.0, tlmat2inv, tlmat1,
                       0.0, tlmat3);


  gsl_blas_dgemm (CblasTrans, CblasNoTrans,
                       1.0, tlmat1, tlmat3,
                       0.0, ctempmatrix);

//std::cout<<"gsl_blas_dgemm done"<<std::endl;
  for(int i=0; i<npt; i++)
  {
    for(int j=0; j<npt; j++)
    {
  
      gsl_matrix_set(cinvperpmatrix, i,j, gsl_matrix_get(cinvmatrix,i,j)   -  gsl_matrix_get(ctempmatrix,i,j)   );
  
    }
  
  }
  
  //calculate chi^2
  //        chi = y^T[invC - invC L (L^T invC L)^(-1) L^T invC] y

  chisquare = 0.0;
  
 
  for(int i=0; i<npt; i++)
    {
      for(int j=0; j<npt; j++)
      {
	//chisquare+= (gsl_vector_get(contband,i))*(gsl_vector_get(contband,j))*   gsl_matrix_get(cinvperpmatrix,i,j);
	chisquare+= (gsl_vector_get(contband,i))*(gsl_vector_get(contband,j))* (gsl_matrix_get(cinvmatrix,i,j)-gsl_matrix_get(ctempmatrix,i,j));
      }
    }
    
    //the determinants were already taken log_e in the subroutine logdet

    loglikelihood = - 0.5*(logdet_c+logdet_tlmat2) - 0.5*chisquare;

  
    }
 }
   
  
  
  for(int i=0;i<lpar->size;i++)
  {
    
    gsl_vector_set(meanmag_per_band_vec,i,gsl_vector_get(lpar,i));
    
  }
  

  gsl_matrix_free(ctempmatrix);
  gsl_matrix_free(cmatrix);
  gsl_matrix_free(cinvmatrix);
  gsl_matrix_free(cinvperpmatrix);
  gsl_matrix_free(tlmat1);
  gsl_matrix_free(tlmat3);
  gsl_matrix_free(tlmat5);
  gsl_matrix_free(lmat);
  gsl_matrix_free(tlmat2);
  gsl_matrix_free(tlmat2inv);
  gsl_vector_free(linpar);
  gsl_vector_free(lpar);
   
//std::cout<<"loglikelihood "<<loglikelihood<<std::endl;
  return loglikelihood;
}

  
void calculate_classification_params(string objid, gsl_vector *&gband, gsl_vector *&rband, gsl_vector *&iband, gsl_vector *&zband, gsl_vector *&yband, 
				     gsl_vector *&gbanderr, gsl_vector *&rbanderr, gsl_vector *&ibanderr, gsl_vector *&zbanderr, gsl_vector *&ybanderr,
				     gsl_vector *&gbandtime_vec, gsl_vector *&rbandtime_vec, gsl_vector *&ibandtime_vec, gsl_vector *&zbandtime_vec, 
				     gsl_vector *&ybandtime_vec, 
				     vector<int>& gbandflag_vec,  vector<int>& rbandflag_vec,
				     vector<int>& ibandflag_vec, vector<int>& zbandflag_vec, vector<int>& ybandflag_vec,
				     gsl_vector *&gband_psf_qf_perfect_vec, gsl_vector *&rband_psf_qf_perfect_vec, gsl_vector *&iband_psf_qf_perfect_vec, gsl_vector *&zband_psf_qf_perfect_vec,
				     gsl_vector *&yband_psf_qf_perfect_vec,
				     gsl_vector *&gband_ap_mag, gsl_vector *&rband_ap_mag, gsl_vector *&iband_ap_mag, gsl_vector *&zband_ap_mag, 
				     gsl_vector *&yband_ap_mag,
				     gsl_vector *&gband_psf_inst_mag, gsl_vector *&rband_psf_inst_mag, gsl_vector *&iband_psf_inst_mag, 
				     gsl_vector *&zband_psf_inst_mag, gsl_vector *&yband_psf_inst_mag, 
				     double ra, double dec,
				     int &processed,
				     double &ra_output, double &dec_output,
				     double &maxloglik_output, double &q_output, int &stellar_locus_output, double &omega_best_output,
				     double &tau_best_output, double &errorweightedmeanmag_g_output, double &errorweightedmeanmag_r_output,
				     double &errorweightedmeanmag_i_output, double &errorweightedmeanmag_z_output, double &errorweightedmeanmag_y_output,
				     double &strucfuncmeanmag_g_output, double &strucfuncmeanmag_r_output, double &strucfuncmeanmag_i_output,
				     double &strucfuncmeanmag_z_output, double &strucfuncmeanmag_y_output,
				     int &gbandepochs_cleaned_output, int &rbandepochs_cleaned_output,
				     int &ibandepochs_cleaned_output,  int &zbandepochs_cleaned_output,
				     int &ybandepochs_cleaned_output)
{
  
processed=0;
  //process all light curves with:
  //(after outlier cleaning, q>5 off stellar locus, q>30 within stellar locus)
  //after outlier cleaning, at least 10 epochs
  
  //output: objid ra dec longitude_l latitude_b maxloglik q stellar_locus omega_best tau_best meanmag_g meanmag_r meanmag_i meanmag_z meanmag_y gbandepochs rbandepochs ibandepochs zbandepochs ybandepochs SDSStype

 
  int fitmodel = 2;   //1: power law, 2: drw

  int returnread=0;
  
  int dim_gband=0;
  int dim_rband=0;
  int dim_iband=0;
  int dim_zband=0;
  int dim_yband=0;
  int dim_gband_cleaned=0;
  int dim_rband_cleaned=0;
  int dim_iband_cleaned=0;
  int dim_zband_cleaned=0;
  int dim_yband_cleaned=0;

  gsl_vector *gband_cleaned=0;
  gsl_vector *rband_cleaned=0;    
  gsl_vector *iband_cleaned=0;
  gsl_vector *zband_cleaned=0;
  gsl_vector *yband_cleaned=0;
  
  gsl_vector *gbanderr_cleaned=0;  
  gsl_vector *rbanderr_cleaned=0;
  gsl_vector *ibanderr_cleaned=0;
  gsl_vector *zbanderr_cleaned=0;
  gsl_vector *ybanderr_cleaned=0;  
  
  gsl_vector *gbandtime_cleaned=0;
  gsl_vector *rbandtime_cleaned=0;
  gsl_vector *ibandtime_cleaned=0;
  gsl_vector *zbandtime_cleaned=0;
  gsl_vector *ybandtime_cleaned=0;
 
  gsl_vector *gband_flags_cleaned=0;
  gsl_vector *rband_flags_cleaned=0;
  gsl_vector *iband_flags_cleaned=0;
  gsl_vector *zband_flags_cleaned=0;
  gsl_vector *yband_flags_cleaned=0;
  /*
  
  vector<int> gband_flags_vec;
  vector<int> rband_flags_vec;
  vector<int> iband_flags_vec;
  vector<int> zband_flags_vec;
  vector<int> yband_flags_vec;
	*/			     
  
  
  int number_gband_out_as_gband_psf_qf_perfect_less0p95,number_rband_out_as_rband_psf_qf_perfect_less0p95;
  int number_iband_out_as_iband_psf_qf_perfect_less0p95,number_zband_out_as_zband_psf_qf_perfect_less0p95;
  int number_yband_out_as_yband_psf_qf_perfect_less0p95,number_gband_out_as_gband_ap_mag,number_rband_out_as_rband_ap_mag;
  int number_iband_out_as_iband_ap_mag;
  int number_zband_out_as_zband_ap_mag,number_yband_out_as_yband_ap_mag,total_beforechi2;
  int number_gband_outchi2,number_rband_outchi2,number_iband_outchi2,number_zband_outchi2,number_yband_outchi2;

  /*
  int number_gband_badflag=0;
  int number_rband_badflag=0;
  int number_iband_badflag=0;
  int number_zband_badflag=0;
  int number_yband_badflag=0;

   */
				
  gsl_vector *contband_vec=0;
  gsl_vector *time_vec=0;
  gsl_vector *m_err_vec=0;
  gsl_vector *bandindex_vec=0;
  gsl_vector *wavelengths_vec=0;
 
 int number_of_bands=5;
   wavelengths_vec=gsl_vector_calloc(5);//grizy 
   
 //PS1
 
   gsl_vector_set(wavelengths_vec,0, 4810); //g
   gsl_vector_set(wavelengths_vec,1, 6170); //r
   gsl_vector_set(wavelengths_vec,2, 7520); //i
   gsl_vector_set(wavelengths_vec,3, 8660); //z
   gsl_vector_set(wavelengths_vec,4, 9620); //y
   
  
  gsl_matrix *m_mat=0;
  
  
  gsl_vector *mean_per_band;
  
  
  gsl_vector *strucfunc_params_vec;
  strucfunc_params_vec = gsl_vector_calloc(3); //omega_r, tau, alpha
  gsl_vector *strucfunc_params_bestfit_vec;
  
//  gsl_vector *strucfunc_par_CI068_min;
//  strucfunc_par_CI068_min = gsl_vector_calloc(3);
  
//  gsl_vector *strucfunc_par_CI095_min;
//  strucfunc_par_CI095_min = gsl_vector_calloc(3);
    

  double logposterior;
  double q;
  int counter=0;
  gsl_vector *missingbands=0;
  int used_number_of_bands=0;
  gsl_vector *meanmag_per_band_vec=0;

  gsl_vector* strucfunc_params = gsl_vector_calloc(3);
  
  size_t maxindex_omegar=0;
  size_t maxindex_tau=0;
  //10x10 grid, uniform in omega, logarithmic in tau

  double omega_r_min;
  double omega_r_max;
  
  double tau_min;
  double tau_max;
  
  if(fitmodel==1) //power law
  {
     omega_r_min=0.01;
     omega_r_max = 1.0;
  
     tau_min = 0.01;
     tau_max = 1.0;
 
  }
if(fitmodel==2)
{
  
  omega_r_min=0.01;
  omega_r_max = 3.16;
  
  tau_min = 0.04;
  tau_max = 5000;   //730;
 
}
  
  int dim_omegardirection=20;//20;
  int dim_taudirection=30;//30;
  
  double stepsize_logomegar=(log10(omega_r_max)-log10(omega_r_min))/(dim_omegardirection-1);
  double stepsize_logtau =(log10(tau_max)-log10(tau_min))/(dim_taudirection-1);

  //double omega_r;
  //double tau;
  double omega_best=0;
  double tau_best=0;
 
  
  double maxloglik;
  
  double longitude_l=0;
  double latitude_b=0;
  int stellar_locus=0;
  int basic_colorcut_passed=0;
  
  double errorweightedmean_g=0;
  double errorweightedmean_r=0;
  double errorweightedmean_i=0;
  double errorweightedmean_z=0; 
  double errorweightedmean_y=0;
  
  
  int epochs_original =0;
  
  bool good_ap_mag = false;
  bool good_psf_qf_perfect = false;

  int epochs=0;
  
  //make structure function parameter grid (fill with omega_min to omega_max, tau_min to tau_max) par_grid
  
  gsl_matrix *par_grid_omegar = gsl_matrix_alloc(dim_omegardirection,dim_taudirection);
  gsl_matrix *par_grid_tau = gsl_matrix_alloc(dim_omegardirection,dim_taudirection);    
  
  for(int i=0;i<dim_omegardirection;i++)
  {
    for(int j=0;j<dim_taudirection;j++)
    {
      gsl_matrix_set(par_grid_omegar,i,j,pow(10,log10(omega_r_min)+i*stepsize_logomegar));
      gsl_matrix_set(par_grid_tau,i,j, pow(10,log10(tau_min)+j*stepsize_logtau));
    }
  }
   
  gsl_matrix *loglik_grid = gsl_matrix_alloc(dim_omegardirection,dim_taudirection);
  
  
 
      
    
    //int dim = 0;
    
    
//std::cout<<"initialized"<<std::endl;
  //  int flag=0;

    
     
  longitude_l=0;
  latitude_b=0;

  errorweightedmean_g=0;
  errorweightedmean_r=0;
  errorweightedmean_i=0;
  errorweightedmean_z=0;
  errorweightedmean_y=0;
  q=0;
	    
     
  basic_colorcut_passed=0;
  stellar_locus=0;

  
     //apply chi-based outlier cleaning
    
	   
	   returnread=cleanPS1file(gband, rband, iband, zband,yband, gbanderr, rbanderr, ibanderr, zbanderr, ybanderr,
				  gbandtime_vec, rbandtime_vec, ibandtime_vec, zbandtime_vec, ybandtime_vec,
				  gbandflag_vec, rbandflag_vec, ibandflag_vec, zbandflag_vec, ybandflag_vec,
				  gband_psf_qf_perfect_vec, rband_psf_qf_perfect_vec, iband_psf_qf_perfect_vec, zband_psf_qf_perfect_vec,
				  yband_psf_qf_perfect_vec,
				  gband_ap_mag, rband_ap_mag, iband_ap_mag, zband_ap_mag, yband_ap_mag,
				  gband_psf_inst_mag, rband_psf_inst_mag, iband_psf_inst_mag, zband_psf_inst_mag, yband_psf_inst_mag, 
				  gband_cleaned, rband_cleaned, iband_cleaned, zband_cleaned, yband_cleaned,
				  gbanderr_cleaned, rbanderr_cleaned, ibanderr_cleaned, zbanderr_cleaned, ybanderr_cleaned,
				  gbandtime_cleaned, rbandtime_cleaned, ibandtime_cleaned, zbandtime_cleaned, 
				  ybandtime_cleaned,
				 /* number_gband_badflag, number_rband_badflag, number_iband_badflag,  number_zband_badflag,  number_yband_badflag,*/
				  /*gband_badflags, rband_badflags, iband_badflags, zband_badflags, yband_badflags, */
				  number_gband_out_as_gband_psf_qf_perfect_less0p95,
				  number_rband_out_as_rband_psf_qf_perfect_less0p95,
				  number_iband_out_as_iband_psf_qf_perfect_less0p95,
				  number_zband_out_as_zband_psf_qf_perfect_less0p95,
				  number_yband_out_as_yband_psf_qf_perfect_less0p95,
				  number_gband_out_as_gband_ap_mag,number_rband_out_as_rband_ap_mag,number_iband_out_as_iband_ap_mag,
				  number_zband_out_as_zband_ap_mag,number_yband_out_as_yband_ap_mag,total_beforechi2,
			          number_gband_outchi2, number_rband_outchi2, number_iband_outchi2, number_zband_outchi2, 
				  number_yband_outchi2,0);     //1: flags only, 0: also chi^2
	 
		      if(returnread==0)
		      {
			
   
			good_ap_mag = false;
			good_psf_qf_perfect = false;
			
			epochs_original = 0;
 
			if(gbandtime_vec!=0)
			{	
			  epochs_original+=gbandtime_vec->size;
			}
			if(rbandtime_vec!=0)
			{	
			  epochs_original+=rbandtime_vec->size;
			}
			if(ibandtime_vec!=0)
			{	
			  epochs_original+=ibandtime_vec->size;
			}
			if(zbandtime_vec!=0)
			{	
			  epochs_original+=zbandtime_vec->size;
			}
			if(ybandtime_vec!=0)
			{	
			  epochs_original+=ybandtime_vec->size;
			}
			
   //out as too many bad by ap_mag?
   //out as too many bad by psf_qf_perfect_less0p95?
   
   if(number_gband_out_as_gband_psf_qf_perfect_less0p95 + number_rband_out_as_rband_psf_qf_perfect_less0p95 + 
     number_iband_out_as_iband_psf_qf_perfect_less0p95 + number_zband_out_as_zband_psf_qf_perfect_less0p95 + number_yband_out_as_yband_psf_qf_perfect_less0p95
     < (epochs_original/4))
   {good_psf_qf_perfect = true;};
   
   if(number_gband_out_as_gband_ap_mag + number_rband_out_as_rband_ap_mag + 
     number_iband_out_as_iband_ap_mag + number_zband_out_as_zband_ap_mag + number_yband_out_as_yband_ap_mag
     < (epochs_original/4))
   {good_ap_mag = true;};
		  
	if(good_ap_mag && good_psf_qf_perfect)
	{
			
			
 calculate_q(gband_cleaned, rband_cleaned, iband_cleaned, zband_cleaned, yband_cleaned,
		       gbanderr_cleaned, rbanderr_cleaned, ibanderr_cleaned, zbanderr_cleaned, ybanderr_cleaned,
                       q, errorweightedmean_g, errorweightedmean_r, errorweightedmean_i, 
		       errorweightedmean_z, errorweightedmean_y,
		       basic_colorcut_passed, stellar_locus);
 
  //after outlier cleaning
if(q<1e32)
{ 
     if(              ( (errorweightedmean_g > 15) && (errorweightedmean_g < 21.5) )
                    ||( (errorweightedmean_r > 15) && (errorweightedmean_r < 21.5) )
		    ||( (errorweightedmean_i > 15) && (errorweightedmean_i < 21.5) )     )
     
{
  
      			 dim_gband=0;	  
			 dim_rband=0;
			 dim_iband=0;
			 dim_zband=0;
			 dim_yband=0;
						
			 dim_gband_cleaned=0;	  
			 dim_rband_cleaned=0;
			 dim_iband_cleaned=0;
			 dim_zband_cleaned=0;
			 dim_yband_cleaned=0;
			
			missingbands = gsl_vector_calloc(number_of_bands);
			 
			used_number_of_bands=0;
			epochs=0;	  
			if(gbandtime_cleaned!=0)
			{	
			  epochs+=gbandtime_cleaned->size;
			  used_number_of_bands++;
			  dim_gband_cleaned = gbandtime_cleaned->size;
			}
			else
			{
			  gsl_vector_set(missingbands,0,1);
			}
			if(rbandtime_cleaned!=0)
			{	
			  epochs+=rbandtime_cleaned->size;
			  used_number_of_bands++;
			  dim_rband_cleaned = rbandtime_cleaned->size;
			}
			else
			{
			  gsl_vector_set(missingbands,1,1);
			}
			if(ibandtime_cleaned!=0)
			{	
			  epochs+=ibandtime_cleaned->size;
			  used_number_of_bands++;
			  dim_iband_cleaned = ibandtime_cleaned->size;
			}
			else
			{
			  gsl_vector_set(missingbands,2,1);
			}
			if(zbandtime_cleaned!=0)
			{	
			  epochs+=zbandtime_cleaned->size;
			  used_number_of_bands++;
			  dim_zband_cleaned = zbandtime_cleaned->size;
			}
			else
			{
			  gsl_vector_set(missingbands,3,1);
			}
			if(ybandtime_cleaned!=0)
			{	
			  epochs+=ybandtime_cleaned->size;
			  used_number_of_bands++;
			  dim_yband_cleaned = ybandtime_cleaned->size;
			}
			else
			{
			  gsl_vector_set(missingbands,4,1);
			}
					
//after outlier cleaning, at least 10 epochs
 if(epochs>=10)
 {
 		
     
 //  i) calculate loglik on parameter grid
   //DRW
   
	    //fill contband_vec, time_vec, bandindex_vec
	    
	    contband_vec = gsl_vector_calloc(epochs);
	    m_err_vec = gsl_vector_calloc(epochs);
	    time_vec = gsl_vector_calloc(epochs);
	     bandindex_vec = gsl_vector_calloc(epochs);
	     
	  
	    counter=0;
	    
	    if(gbandtime_cleaned!=0)
	    {
	      for(int i=0;i<gbandtime_cleaned->size;i++)
	      {
		gsl_vector_set(contband_vec,counter, gsl_vector_get(gband_cleaned,i));
		gsl_vector_set(m_err_vec,counter, gsl_vector_get(gbanderr_cleaned,i));
		gsl_vector_set(time_vec, counter, gsl_vector_get(gbandtime_cleaned,i));
		gsl_vector_set(bandindex_vec,i,0);
		counter++;
	      }
	    }
	    if(rbandtime_cleaned!=0)
	    {
	      for(int i=0;i<rbandtime_cleaned->size;i++)
	      {
		gsl_vector_set(contband_vec,counter, gsl_vector_get(rband_cleaned,i));
		gsl_vector_set(m_err_vec,counter, gsl_vector_get(rbanderr_cleaned,i));
		gsl_vector_set(time_vec, counter, gsl_vector_get(rbandtime_cleaned,i));
		gsl_vector_set(bandindex_vec,counter,1);
		counter++;
	      }
	    }
	   
	   if(ibandtime_cleaned!=0)
	   {
	    for(int i=0;i<ibandtime_cleaned->size;i++)
	    {
	      gsl_vector_set(contband_vec,counter, gsl_vector_get(iband_cleaned,i));
	      gsl_vector_set(m_err_vec,counter, gsl_vector_get(ibanderr_cleaned,i));
	      gsl_vector_set(time_vec, counter, gsl_vector_get(ibandtime_cleaned,i));
	      gsl_vector_set(bandindex_vec,counter,2);
	      counter++;
	    }
	   }
	   if(zbandtime_cleaned!=0)
	   {
	    for(int i=0;i<zbandtime_cleaned->size;i++)
	    {
	      gsl_vector_set(contband_vec,counter, gsl_vector_get(zband_cleaned,i));
	      gsl_vector_set(m_err_vec,counter, gsl_vector_get(zbanderr_cleaned,i));
	      gsl_vector_set(time_vec, counter, gsl_vector_get(zbandtime_cleaned,i));
	      gsl_vector_set(bandindex_vec,counter,3);
	      counter++;
	    }
	   }
	   if(ybandtime_cleaned!=0)
	   {
	    for(int i=0;i<ybandtime_cleaned->size;i++)
	    {
	      gsl_vector_set(contband_vec,counter, gsl_vector_get(yband_cleaned,i));
	      gsl_vector_set(m_err_vec,counter, gsl_vector_get(ybanderr_cleaned,i));
	      gsl_vector_set(time_vec, counter, gsl_vector_get(ybandtime_cleaned,i));
	      gsl_vector_set(bandindex_vec,counter,4);
	      counter++;
	    }
	   }

	
	    m_mat = gsl_matrix_alloc(epochs,used_number_of_bands);//rows, columns
  
	    
	    //set m_mat M
	    
	    gsl_matrix_set_zero(m_mat);

	    if(used_number_of_bands== number_of_bands)
	    {
	      for(int i=0;i<epochs;i++)
	      {
		gsl_matrix_set(m_mat,i,gsl_vector_get(bandindex_vec,i),1);
		
	      }
	    
	    }
	    else
	    {
	    int thebandindex;
	      
	    for(int i=0;i<epochs;i++)
	    {
	     	  
		thebandindex=gsl_vector_get(bandindex_vec,i);
	       
		 
		 
	       if(used_number_of_bands!=number_of_bands)
	      {
		if(gsl_vector_get(missingbands,0)==1) //no g band available
		{
		  thebandindex--;
		}
		if(gsl_vector_get(missingbands,1)==1) //no r band available
		{
		  if(gsl_vector_get(bandindex_vec,i)>1)
		  {
		    thebandindex--;
		  }
		}
		if(gsl_vector_get(missingbands,2)==1) //no i band available
		{
		  if(gsl_vector_get(bandindex_vec,i)>2)
		  {
		    thebandindex--;
		  }
		}
		if(gsl_vector_get(missingbands,3)==1) //no z band available
		{
		  if(gsl_vector_get(bandindex_vec,i)>3)
		  {
		    thebandindex--;
		  }
		}
		if(gsl_vector_get(missingbands,4)==1) //no y band available
		{
		  if(gsl_vector_get(bandindex_vec,i)>4)
		  {
		    thebandindex--;
		  }
		}
		
	      }
	       
	      gsl_matrix_set(m_mat,i,thebandindex,1);
	      
	    }
	    }
 
	    
	    
	     meanmag_per_band_vec=gsl_vector_calloc(used_number_of_bands);
 
	
	     
	//std::cout<<"dim_parametergrid "<<dim_parametergrid<<std::endl;
	
	     
	     
	for(int i=0; i<dim_omegardirection;i++)
	{
	  for(int j=0;j<dim_taudirection;j++)
	  {
//pass parameter   
	    
//	std::cout<<i<<" "<<dim_parametergrid<<std::endl;    	
	      
	   
	  gsl_vector_set(strucfunc_params, 0, gsl_matrix_get(par_grid_omegar,i,j));
	  gsl_vector_set(strucfunc_params, 1, gsl_matrix_get(par_grid_tau,i,j));
	  gsl_vector_set(strucfunc_params, 2, -0.65);    //alphavalue
 
	  logposterior = loglikelihood_strucparam_multicolor(contband_vec, bandindex_vec, 
								  m_mat, meanmag_per_band_vec,
								  time_vec, m_err_vec, strucfunc_params, wavelengths_vec, fitmodel, missingbands);
	 //here, we use flat prior, so don't add any logprior
	
	   //write into loglik_grid
	   gsl_matrix_set(loglik_grid,i,j,logposterior);
	}
	}
	
	
	//get the one with highest loglik_grid entry
	//std::cout<<"0"<<std::endl;
	gsl_matrix_max_index(loglik_grid,&maxindex_omegar,&maxindex_tau);
	
	maxloglik = gsl_matrix_get(loglik_grid,maxindex_omegar,maxindex_tau);
	omega_best = gsl_matrix_get(par_grid_omegar,maxindex_omegar,maxindex_tau);
	tau_best = gsl_matrix_get(par_grid_tau,maxindex_omegar,maxindex_tau);
  

  
	
	processed=1;
	ra_output=ra;
	dec_output=dec;

	maxloglik_output=maxloglik;
	q_output=q;
	stellar_locus_output=stellar_locus;
	omega_best_output=omega_best;
	tau_best_output=tau_best;
	errorweightedmeanmag_g_output=errorweightedmean_g;
	errorweightedmeanmag_r_output=errorweightedmean_r;
	errorweightedmeanmag_i_output=errorweightedmean_i;
	errorweightedmeanmag_z_output=errorweightedmean_z;
	errorweightedmeanmag_y_output=errorweightedmean_y;
	
	int missingbandindex=0;
	int meanmagindex=0;
	
	

	  if (gsl_vector_get(missingbands,0)==0) // g band is available
	  {
	    strucfuncmeanmag_g_output = gsl_vector_get(meanmag_per_band_vec,meanmagindex);
	    meanmagindex++;
	  }
	  else
	  {strucfuncmeanmag_g_output=0;}
	  

	  if (gsl_vector_get(missingbands,1)==0) // r band is available
	  {
	    strucfuncmeanmag_r_output = gsl_vector_get(meanmag_per_band_vec,meanmagindex);
	    meanmagindex++;
	  }
	  else
	  {strucfuncmeanmag_r_output=0;}
	  
	  
	  if (gsl_vector_get(missingbands,2)==0) // i band is available
	  {
	    strucfuncmeanmag_i_output = gsl_vector_get(meanmag_per_band_vec,meanmagindex);
	    meanmagindex++;
	  }
	  else
	  {strucfuncmeanmag_i_output=0;}
	  
	  
	  if (gsl_vector_get(missingbands,3)==0) // z band is available
	  {
	    strucfuncmeanmag_z_output = gsl_vector_get(meanmag_per_band_vec,meanmagindex);
	    meanmagindex++;
	  }
	  else
	  {strucfuncmeanmag_z_output=0;}
	  
	  if (gsl_vector_get(missingbands,4)==0) // y band is available
	  {
	    strucfuncmeanmag_y_output = gsl_vector_get(meanmag_per_band_vec,meanmagindex);
	    meanmagindex++;
	  }
	  else
	  {strucfuncmeanmag_y_output=0;}
	  
	
	
	
	gbandepochs_cleaned_output=dim_gband_cleaned;
	rbandepochs_cleaned_output=dim_rband_cleaned;
	ibandepochs_cleaned_output=dim_iband_cleaned;
	zbandepochs_cleaned_output=dim_zband_cleaned;
	ybandepochs_cleaned_output=dim_yband_cleaned;
				      
				      
				      
	
	     }
	

}}}
			//clean




  	if(gband_cleaned!=0)
	{
	  gsl_vector_free(gband_cleaned);
	}
	if(rband_cleaned!=0)
	{
	  gsl_vector_free(rband_cleaned);
	}
	if(iband_cleaned!=0)
	{
	  gsl_vector_free(iband_cleaned);
	}
	if(zband_cleaned!=0)
	{
	  gsl_vector_free(zband_cleaned);
	}
	if(yband_cleaned!=0)
	{
	  gsl_vector_free(yband_cleaned);
	}
	
	
	
  
 
	if(gbanderr_cleaned!=0)
	{
	  gsl_vector_free(gbanderr_cleaned);
	}
	if(rbanderr_cleaned!=0)
	{
	  gsl_vector_free(rbanderr_cleaned);
	}
	if(ibanderr_cleaned!=0)
	{
	  gsl_vector_free(ibanderr_cleaned);
	}
	if(zbanderr_cleaned!=0)
	{
	  gsl_vector_free(zbanderr_cleaned);
	}
	if(ybanderr_cleaned!=0)
	{
	  gsl_vector_free(ybanderr_cleaned);
	}

  
	
	if(gbandtime_cleaned!=0)
	{
	  gsl_vector_free(gbandtime_cleaned);
	}
	if(rbandtime_cleaned!=0)
	{
	  gsl_vector_free(rbandtime_cleaned);
	}
	if(ibandtime_cleaned!=0)
	{
	  gsl_vector_free(ibandtime_cleaned);
	}
	if(zbandtime_cleaned!=0)
	{
	  gsl_vector_free(zbandtime_cleaned);
	}
	if(ybandtime_cleaned!=0)
	{
	  gsl_vector_free(ybandtime_cleaned);
	}
	if(gband_flags_cleaned!=0)
	{
	  gsl_vector_free(gband_flags_cleaned);
	}

	
	
	
	if(rband_flags_cleaned!=0)
	{
	  gsl_vector_free(rband_flags_cleaned);
	}
	if(iband_flags_cleaned!=0)
	{
	  gsl_vector_free(iband_flags_cleaned);
	}
	if(zband_flags_cleaned!=0)
	{
	  gsl_vector_free(zband_flags_cleaned);
	}
	if(yband_flags_cleaned!=0)
	{
	  gsl_vector_free(yband_flags_cleaned);
	}
	
	
  
	if(contband_vec!=0)
	{
	  gsl_vector_free(contband_vec);
	}
	if(time_vec!=0)
	{
	  gsl_vector_free(time_vec);
	}
	if(m_err_vec!=0)
	{
	  gsl_vector_free(m_err_vec);
	}
	if(mean_per_band!=0)
	{
	  gsl_vector_free(mean_per_band);
	}
	if(wavelengths_vec!=0)
	{
	  gsl_vector_free(wavelengths_vec);
	}
	if(bandindex_vec!=0)
	{
	  gsl_vector_free(bandindex_vec);
	}
	
	
	
	
	if(strucfunc_params_vec!=0)
	{
	  gsl_vector_free(strucfunc_params_vec);
	}
	if(strucfunc_params_bestfit_vec!=0)
	{
	  gsl_vector_free(strucfunc_params_bestfit_vec);
	}
	if(missingbands!=0)
	{
	  gsl_vector_free(missingbands);
	}
	if(meanmag_per_band_vec!=0)
	{
	  gsl_vector_free(meanmag_per_band_vec);
	}
	if(strucfunc_params!=0)
	{
	  gsl_vector_free(strucfunc_params);
	}
	
			
	
  		
			
  
	if(m_mat!=0)
	{
	  gsl_matrix_free(m_mat);
	}
	
	if(par_grid_omegar!=0)
	{
	  gsl_matrix_free(par_grid_omegar);
	}
	
	if(par_grid_tau!=0)
	{
	  gsl_matrix_free(par_grid_tau);
	}
	
	if(loglik_grid!=0)
	{
	  gsl_matrix_free(loglik_grid);
	}
	
 
gband_cleaned=0;
rband_cleaned=0;    
iband_cleaned=0;
zband_cleaned=0;
yband_cleaned=0;
  
gbanderr_cleaned=0;  
rbanderr_cleaned=0;
ibanderr_cleaned=0;
zbanderr_cleaned=0;
ybanderr_cleaned=0;  
  
gbandtime_cleaned=0;
rbandtime_cleaned=0;
ibandtime_cleaned=0;
zbandtime_cleaned=0;
ybandtime_cleaned=0;
 
gband_flags_cleaned=0;
rband_flags_cleaned=0;
iband_flags_cleaned=0;
zband_flags_cleaned=0;
yband_flags_cleaned=0;
	
contband_vec=0;
time_vec=0;
m_err_vec=0;
bandindex_vec=0;
wavelengths_vec=0;
m_mat=0;
  
  
  
strucfunc_params_vec=0;
strucfunc_params_bestfit_vec=0;
missingbands=0;
meanmag_per_band_vec=0;

strucfunc_params=0;
par_grid_omegar=0;
par_grid_tau=0;
loglik_grid = 0;
  
			
		      }




  
				      
				      

}


  