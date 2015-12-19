#!/usr/bin/env python
import os


import numpy as numpy
from scipy.weave import inline


import pyfits
import ubercal_flat_bigflat_chmask

import esutil
from esutil.coords import eq2gal
restart =  True

from matplotlib.mlab import rec_append_fields

flags = {
'DEFAULT'          	  : 0x00000000, # Initial value: resets all bits
'PSFMODEL'         	  : 0x00000001, # Source fitted with a psf model (linear or non-linear)
'EXTMODEL'         	  : 0x00000002, # Source fitted with an extended-source model
'FITTED'           	  : 0x00000004, # Source fitted with non-linear model (PSF or EXT; good or bad)
'FITFAIL'             	  : 0x00000008, # Fit (non-linear) failed (non-converge, off-edge, run to zero)
'POORFIT'             	  : 0x00000010, # Fit succeeds, but low-SN, high-Chisq, or large (for PSF -- drop?)
'PAIR'             	  : 0x00000020, # Source fitted with a double psf
'PSFSTAR'          	  : 0x00000040, # Source used to define PSF model
'SATSTAR'          	  : 0x00000080, # Source model peak is above saturation
'BLEND'            	  : 0x00000100, # Source is a blend with other sources
'EXTERNALPOS'         	  : 0x00000200, # Source based on supplied input position
'BADPSF'           	  : 0x00000400, # Failed to get good estimate of object's PSF
'DEFECT'           	  : 0x00000800, # Source is thought to be a defect
'SATURATED'        	  : 0x00001000, # Source is thought to be saturated pixels (bleed trail)
'CR_LIMIT'         	  : 0x00002000, # Source has crNsigma above limit
'EXT_LIMIT'        	  : 0x00004000, # Source has extNsigma above limit
'MOMENTS_FAILURE'  	  : 0x00008000, # could not measure the moments
'SKY_FAILURE'      	  : 0x00010000, # could not measure the local sky
'SKYVAR_FAILURE'   	  : 0x00020000, # could not measure the local sky variance
'BELOW_MOMENTS_SN' 	  : 0x00040000, # moments not measured due to low S/N
'BIG_RADIUS'       	  : 0x00100000, # poor moments for small radius, try large radius
'AP_MAGS'          	  : 0x00200000, # source has an aperture magnitude
'BLEND_FIT'        	  : 0x00400000, # source was fitted as a blend
'EXTENDED_FIT'     	  : 0x00800000, # full extended fit was used
'EXTENDED_STATS'   	  : 0x01000000, # extended aperture stats calculated
'LINEAR_FIT'       	  : 0x02000000, # source fitted with the linear fit
'NONLINEAR_FIT'    	  : 0x04000000, # source fitted with the non-linear fit
'RADIAL_FLUX'      	  : 0x08000000, # radial flux measurements calculated
'SIZE_SKIPPED'     	  : 0x10000000, # size could not be determined
'ON_SPIKE'         	  : 0x20000000, # peak lands on diffraction spike
'ON_GHOST'         	  : 0x40000000, # peak lands on ghost or glint
'OFF_CHIP'         	  : 0x80000000, # peak lands off edge of chip
}
for x,y in flags.items():
    flags[y] = x

# definition of flags2 (see http://svn.pan-starrs.ifa.hawaii.edu/trac/ipp/wiki/CMF_PS1_V3)
flags2 = {
'DEFAULT'          	  : 0x00000000, # Initial value: resets all bits
'DIFF_WITH_SINGLE' 	  : 0x00000001, # diff source matched to a single positive detection
'DIFF_WITH_DOUBLE' 	  : 0x00000002, # diff source matched to positive detections in both images
'MATCHED'          	  : 0x00000004, # diff source matched to positive detections in both images
'ON_SPIKE'         	  : 0x00000008, # > 25% of (PSF-weighted) pixels land on diffraction spike
'ON_STARCORE'      	  : 0x00000010, # > 25% of (PSF-weighted) pixels land on starcore
'ON_BURNTOOL'      	  : 0x00000020, # > 25% of (PSF-weighted) pixels land on burntool
'ON_CONVPOOR'      	  : 0x00000040, # > 25% of (PSF-weighted) pixels land on convpoor
'PASS1_SRC'               : 0x00000080, # source detected in first pass analysis
'HAS_BRIGHTER_NEIGHBOR'   : 0x00000100, # peak is not the brightest in its footprint
'BRIGHT_NEIGHBOR_1'       : 0x00000200, # flux_n / (r^2 flux_p) > 1
'BRIGHT_NEIGHBOR_10'      : 0x00000400, # flux_n / (r^2 flux_p) > 10
'DIFF_SELF_MATCH'  	  : 0x00000800, # positive detection match is probably this source 
'SATSTAR_PROFILE'         : 0x00001000, # saturated source is modeled with a radial profile
}
for x,y in flags2.items():
    flags2[y] = x


    
def unique(obj):
    """Gives an array indexing the unique elements in the sorted list obj.

    Returns the list of indices of the last elements in each group
    of equal-comparing items in obj
    """
    nobj = len(obj)
    if nobj == 0:
        return numpy.zeros(0, dtype='i8')
    if nobj == 1:
        return numpy.zeros(1, dtype='i8')
    out = numpy.zeros(nobj, dtype=numpy.bool)
    out[0:nobj-1] = (obj[0:nobj-1] != obj[1:nobj])
    out[nobj-1] = True
    return numpy.sort(numpy.flatnonzero(out))
    

class subslices:
    "Iterator for looping over subsets of an array"
    def __init__(self, data, uind=None):
        if uind is None:
            self.uind = unique(data)
        else:
            self.uind = uind.copy()
        self.ind = 0
    def __iter__(self):
        return self
    def __len__(self):
        return len(self.uind)
    def next(self):
        if self.ind == len(self.uind):
            raise StopIteration
        if self.ind == 0:
            first = 0
        else:
            first = self.uind[self.ind-1]+1
        last = self.uind[self.ind]+1
        self.ind += 1
        return first, last
#        return slice(first, self.uind[self.ind]+1)

    

                                   

def call_c_code(all_rows_input, code, support_code):
   
		            
    # Sort by object ID and then by MJD
    if 'obj_id' in all_rows_input.dtype.names:
   
   
        # Prepare the output array 
      
        objs, idx = numpy.unique(all_rows_input['obj_id'], return_index=True) #this gives the number of objects in all_rows_input, should be 1
         
        obj_id, psf_inst_mag, psf_inst_mag_sig, ra, dec, detra, detdec, filterid, mjd_obs, x_psf, y_psf, \
        ap_mag, psf_qf, psf_qf_perfect, ucalmag, flags, flags2, zp = ( all_rows_input[attr].copy().astype(all_rows_input[attr].dtype) for attr in ['obj_id', 'psf_inst_mag', 'psf_inst_mag_sig', 'ra', 'dec', 'detra', 'detdec', 'filterid', 'mjd_obs', 'x_psf', 'y_psf', 'ap_mag', 'psf_qf', 'psf_qf_perfect', 'ucalmag', 'flags', 'flags2', 'zp'] )
  
        
        out = numpy.zeros(len(objs),
            dtype=[
                ('processed_flag', 'i2'),
                ('obj_id_out', 'u8'),         #set i8 to get signed (like in fits files), set u8 to get unsigned (like in LSD)
                ('ra_out', 'f8'),
                ('dec_out', 'f8'),
                ('maxloglik_out', 'f8'),
                ('q_out', 'f8'),
                ('stellar_locus_out', 'i2'),
                ('omega_best_out', 'f8'),
                ('tau_best_out', 'f8'),
                ('errorweightedmeanmag_g_out', 'f8'),
                ('errorweightedmeanmag_r_out', 'f8'),
                ('errorweightedmeanmag_i_out', 'f8'),
                ('errorweightedmeanmag_z_out', 'f8'),
                ('errorweightedmeanmag_y_out', 'f8'),
                ('strucfuncmeanmag_g_out', 'f8'),
                ('strucfuncmeanmag_r_out', 'f8'),
                ('strucfuncmeanmag_i_out', 'f8'),
                ('strucfuncmeanmag_z_out', 'f8'),
                ('strucfuncmeanmag_y_out', 'f8'),
                ('gbandepochs_cleaned_out', 'i2'),
                ('rbandepochs_cleaned_out', 'i2'),
                ('ibandepochs_cleaned_out', 'i2'),
                ('zbandepochs_cleaned_out', 'i2'),
                ('ybandepochs_cleaned_out', 'i2')                
                
            ]
                
            )

        (processed_flag, obj_id_out, ra_out, dec_out, maxloglik_out, q_out
        , stellar_locus_out, omega_best_out, tau_best_out, 
        errorweightedmeanmag_g_out, errorweightedmeanmag_r_out, 
        errorweightedmeanmag_i_out, errorweightedmeanmag_z_out
        , errorweightedmeanmag_y_out, strucfuncmeanmag_g_out, strucfuncmeanmag_r_out, strucfuncmeanmag_i_out, strucfuncmeanmag_z_out
        , strucfuncmeanmag_y_out, 
        gbandepochs_cleaned_out, rbandepochs_cleaned_out, ibandepochs_cleaned_out
        , zbandepochs_cleaned_out, ybandepochs_cleaned_out) = [out[name0].copy() for name0 in out.dtype.names]
       
       
        # Convert filterid to index
        band = numpy.empty(len(all_rows_input), dtype='u2')
        
  
        for i in range(0,len(all_rows_input)):
          if (filterid[i] == 'g.0000'):
            band[i]=0
          if (filterid[i] == 'r.0000'):
            band[i]=1
          if (filterid[i] == 'i.0000'):
            band[i]=2
          if (filterid[i] == 'z.0000'):
            band[i]=3
          if (filterid[i] == 'y.0000'):
            band[i]=4
        
        
#	import pdb
	from scipy import weave
#	pdb.set_trace()
        inline(code,
            ['obj_id', 'psf_inst_mag', 'psf_inst_mag_sig', 'ra', 'dec', 'band', 'mjd_obs','ap_mag', 'psf_qf_perfect', 'ucalmag'
            , 'flags', 'zp'
            , 'processed_flag', 'obj_id_out','ra_out','dec_out','maxloglik_out','q_out'
            , 'stellar_locus_out','omega_best_out','tau_best_out',
            'errorweightedmeanmag_g_out','errorweightedmeanmag_r_out',
            'errorweightedmeanmag_i_out','errorweightedmeanmag_z_out', 'errorweightedmeanmag_y_out',
            'strucfuncmeanmag_g_out', 'strucfuncmeanmag_r_out', 'strucfuncmeanmag_i_out', 'strucfuncmeanmag_z_out',
            'strucfuncmeanmag_y_out',
            'gbandepochs_cleaned_out','rbandepochs_cleaned_out','ibandepochs_cleaned_out'
            , 'zbandepochs_cleaned_out','ybandepochs_cleaned_out'     
            ],
            headers=['<string>', '<sstream>', '<cmath>', '<vector>', '<iostream>', '<algorithm>', '<gsl/gsl_vector.h>', '<gsl/gsl_matrix.h>', '<gsl/gsl_blas.h>', '<gsl/gsl_linalg.h>', '<gsl/gsl_sort_vector.h>', '<gsl/gsl_permute_vector.h>', '<gsl/gsl_sort.h>', '<gsl/gsl_permute.h>', '<gsl/gsl_sort.h>', '<cassert>', '<gsl/gsl_statistics_double.h>'],
            support_code=support_code,
            libraries=['gsl', 'gslcblas'],
            include_dirs=['.','/usr/local/gsl/1.14/include'],      
            verbose=2,
            force=1,  # set force=1 for testing, set force=0 for run
            undef_macros=['NDEBUG'])
            
        	            
	out['processed_flag']=processed_flag
	out['obj_id_out']=obj_id_out
	out['ra_out']=ra_out
	out['dec_out']=dec_out
	out['maxloglik_out']=maxloglik_out
	out['q_out']=q_out
	out['stellar_locus_out']=stellar_locus_out
	out['omega_best_out']=omega_best_out
	out['tau_best_out']=tau_best_out
	out['errorweightedmeanmag_g_out']=errorweightedmeanmag_g_out
	out['errorweightedmeanmag_r_out']=errorweightedmeanmag_r_out
	out['errorweightedmeanmag_i_out']=errorweightedmeanmag_i_out
	out['errorweightedmeanmag_z_out']=errorweightedmeanmag_z_out
	out['errorweightedmeanmag_y_out']=errorweightedmeanmag_y_out
	out['strucfuncmeanmag_g_out']=strucfuncmeanmag_g_out	
	out['strucfuncmeanmag_r_out']=strucfuncmeanmag_r_out
	out['strucfuncmeanmag_i_out']=strucfuncmeanmag_i_out
	out['strucfuncmeanmag_z_out']=strucfuncmeanmag_z_out
	out['strucfuncmeanmag_y_out']=strucfuncmeanmag_y_out
	out['gbandepochs_cleaned_out']=gbandepochs_cleaned_out
	out['rbandepochs_cleaned_out']=rbandepochs_cleaned_out
	out['ibandepochs_cleaned_out']=ibandepochs_cleaned_out
	out['zbandepochs_cleaned_out']=zbandepochs_cleaned_out
	out['ybandepochs_cleaned_out']=ybandepochs_cleaned_out
	
        return out    
            
################            

def structure_function(all_rows_lightcurves, lcfilnam,code, support_code):

	                outfile = open('files_out/maxloglik_params_drw_%s.txt'  % (lcfilnam),'w')   # 'a' means append, 'w' means truncate and write
	 
		        outfile.write('objid ra dec longitude_l latitude_b maxloglik q stellar_locus omega_best tau_best errorweightedmeanmag_g errorweightedmeanmag_r errorweightedmeanmag_i errorweightedmeanmag_z errorweightedmeanmag_y strucfuncmeanmag_g strucfuncmeanmag_r strucfuncmeanmag_i strucfuncmeanmag_z strucfuncmeanmag_y gbandepochs_cleaned rbandepochs_cleaned ibandepochs_cleaned zbandepochs_cleaned ybandepochs_cleaned w1mpro w2mpro w1sigmpro w2sigmpro wisera wisedec sigrawise sigdecwise \n')
			outfile.flush()

			# read in and apply solution
			usol = ubercal_flat_bigflat_chmask.read_flat_solution('ucalqy_bigflat_chmask.fits', chmask=True)
			zp = ubercal_flat_bigflat_chmask.apply_flat_solution(all_rows_lightcurves, usol)
 
		        all_rows_lightcurves = rec_append_fields(all_rows_lightcurves, 'zp', zp.copy())
		        all_rows_lightcurves = rec_append_fields(all_rows_lightcurves, 'ucalmag', all_rows_lightcurves['psf_inst_mag']+zp)
		   
	
	

		        for f, l in subslices(all_rows_lightcurves['obj_id']):
			  
			    out = call_c_code(all_rows_lightcurves[f:l].copy(), code, support_code)
		           

		            if  (out["processed_flag"] == 1) :
		                longitude_l, latitude_b = eq2gal( out["ra_out"], out["dec_out"])
		                outfile.write("%i %.13f %.13f %.13f %.13f %.13f %.13f %i %.13f %.13f %.13f %.13f %.13f %.13f %.13f %.13f %.13f %.13f %.13f %.13f %i %i %i %i %i %.13f %.13f %.13f %.13f %.13f %.13f %.13f %.13f\n" % (out["obj_id_out"], out["ra_out"]
		                , out["dec_out"] , longitude_l, latitude_b , out["maxloglik_out"] , out["q_out"]
		                , out["stellar_locus_out"] , out["omega_best_out"] , out["tau_best_out"] 
		                , out["errorweightedmeanmag_g_out"] , out["errorweightedmeanmag_r_out"] 
		                , out["errorweightedmeanmag_i_out"], out["errorweightedmeanmag_z_out"], out["errorweightedmeanmag_y_out"]
		                , out["strucfuncmeanmag_g_out"], out["strucfuncmeanmag_r_out"] 
		                , out["strucfuncmeanmag_i_out"], out["strucfuncmeanmag_z_out"], out["strucfuncmeanmag_y_out"]
		                , out["gbandepochs_cleaned_out"], out["rbandepochs_cleaned_out"]
		                , out["ibandepochs_cleaned_out"], out["zbandepochs_cleaned_out"], out["ybandepochs_cleaned_out"],
		                all_rows_lightcurves['w1mpro'][f], all_rows_lightcurves['w2mpro'][f], 
		                all_rows_lightcurves['w1sigmpro'][f], all_rows_lightcurves['w2sigmpro'][f], 
		                all_rows_lightcurves['wisera'][f], all_rows_lightcurves['wisedec'][f], 
		                all_rows_lightcurves['sigrawise'][f], all_rows_lightcurves['sigdecwise'][f]
		                ))    
		                outfile.flush()
	  
		        outfile.close()    
		        #add filename to processed_fitsfiles.log
		        logfile = open('processed_fitsfiles.log','a')
		        
		        logfile.write("%s\n" % lcfilnam)
		        logfile.close()
		        
		       

def main():
	import sys
	numpy.seterr(invalid='ignore')
	
	
        # C++ code to connect
        f = open('connect_lsd.cpp')
        code = f.read()
        f.close()

        # support C code
        f = open('estimate_multiband_variability.cpp')
        support_code = f.read()
        f.close()
      
				
	lcfilnam=sys.argv[1]
	
	print 'lcfilnam ', lcfilnam
	
	
	#####################
	
	#get filename from argument
	
	# make a log logfile (it's for all; later, write filename to logfile when file is processed)
	
	#if cell logfile doesn't exist: create it
	logfile = open('processed_fitsfiles.log','a')
	logfile.close()	
	
	
        lightcurves =pyfits.getdata('files_in/'+lcfilnam,0)

        all_rows_lightcurves = numpy.zeros(lightcurves.size,
        dtype=[
            ('psf_inst_mag', 'f8'), ('psf_inst_mag_sig', 'f8'), ('ra', 'f8'),
            ('dec', 'f8'), ('detra', 'f8'), ('detdec', 'f8'), ('filterid', 'S6'),
            ('mjd_obs', 'f8'), ('x_psf', 'f8'), ('y_psf', 'f8'), 
            ('ap_mag', 'f8'), ('psf_qf', 'f8'), ('psf_qf_perfect', 'f8'),
            ('w1mpro', 'f8'), ('w2mpro', 'f8'), ('w1sigmpro', 'f8'), 
            ('w2sigmpro', 'f8'), ('wisera', 'f8'), ('wisedec', 'f8'), ('sigrawise', 'f8'),
            ('sigdecwise', 'f8'), ('flags', 'u8'), ('flags2', 'u8'),('obj_id', 'i8'), ('chip_id', 'i8'),
        ]
        
        )
   
         
         
        for k in lightcurves.dtype.names:
            all_rows_lightcurves[k][:] = lightcurves[k]
        
        
        test = structure_function(all_rows_lightcurves,lcfilnam,code, support_code)
    
    
if __name__ == "__main__":
  main()
  
  


	   
