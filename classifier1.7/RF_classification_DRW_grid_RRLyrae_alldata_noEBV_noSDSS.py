

####################################################################################################################################
# 
# Random Forest Classification for QSO candidates
# 
# to execute the code, use classify_QSO_run_noSDSS.sh
# 
####################################################################################################################################


import os
import numpy as np
import numpy.random
from esutil.numpy_util import match

import learn_class_csc_classweight

from dust import getval
from esutil.coords import eq2gal
from pyfits import getdata
from numpy.lib import recfunctions
import esutil

import cPickle   #need for serialization

from random import randint

import sys



make_trainingset=int(sys.argv[1])
###########################################
# make_trainingset:
# 0	 applying the already trained RFC to classify sources
# 1	 training the classifier on S82, prepare RFC and pickle it
###########################################


print 'make_trainingset ', make_trainingset

if(make_trainingset==0):     #we're not making the training set, we are classifying!
	   indirectory=sys.argv[2]
	   outdirectory=sys.argv[3]
           sdssdirectory =sys.argv[4]
           
	   print 'indirectory ', indirectory
	   print 'outdirectory ', outdirectory
	   print 'sdssdirectory ', sdssdirectory


#PREPARE THE TRAINING SET

# Eddie's dust map
fname = 'ps1-ebv-4.5kpc.fits'
fname = os.path.join(os.environ['DUST_DIR'], fname)

dustmap = getdata(fname)

#the params we need for the classifier
par_keys = ['omega_best', 'tau_best', 'q', 'stellar_locus','gr', 'ri', 'iz', 'zy', 'W12', 'iW1','corrmeanmag_r']


#variability only
#par_keys = ['omega_best', 'tau_best', 'q']

#color only
#par_keys = ['stellar_locus','gr', 'ri', 'iz', 'zy', 'ugSDSS', 'grSDSS', 'riSDSS', 'izSDSS', 'W12', 'iW1', 'corrmeanmag_r', 'EBV']


if(make_trainingset==1):
  
 
#PREPARE THE TRAINING SET

   params_trainingset = np.genfromtxt('/a41233d1/hernitschek/_results_juelich/_results_classifier/_trainingset_S82/results_s82_trainingset_testcase3_5sig.txt', \
   names = 'objid,ra,dec,longitude_l,latitude_b,maxloglik,q,stellar_locus,omega_best,tau_best,errorweightedmeanmag_g,errorweightedmeanmag_r,errorweightedmeanmag_i,errorweightedmeanmag_z,errorweightedmeanmag_y,strucfuncmeanmag_g, strucfuncmeanmag_r, strucfuncmeanmag_i, strucfuncmeanmag_z, strucfuncmeanmag_y, gbandepochs_cleaned,rbandepochs_cleaned,ibandepochs_cleaned,zbandepochs_cleaned,ybandepochs_cleaned,w1mpro,w2mpro,w1sigmpro,w2sigmpro,wisera,wisedec,sigrawise,sigdecwise,SDSStype,errorweightedmeanmag_g_sigma,errorweightedmeanmag_r_sigma,errorweightedmeanmag_i_sigma,errorweightedmeanmag_z_sigma,errorweightedmeanmag_y_sigma', \
   dtype='|S20, f8, f8, f8, f8, f8, f8, u1, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, u2, u2, u2, u2, u2, f8, f8, f8, f8, f8, f8, f8, f8, |S7, f8, f8, f8, f8, f8', skip_header=1)


# restrict to a certain region
# in the same way, one could restrict to a q or color cut...

   good = ((params_trainingset['ra'] > 360-52) | (params_trainingset['ra'] < 60)) & (np.abs(params_trainingset['dec']) < 1.26)
   params_trainingset = params_trainingset[good]

   # load SDSS
   sdss_s82 = np.genfromtxt('/a41233d1/hernitschek/_sdssdata_s82_newsdss/s82_sdss.txt', \
   names = 'objid,ps1_ra,ps1_dec,sdss_ra,sdss_dec,uPsf,gPsf,rPsf,iPsf,zPsf,uPsfErr,gPsfErr,rPsfErr,iPsfErr,zPsfErr', \
   dtype='|S20, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8', skip_header=1)


   good = (     (sdss_s82['uPsf'] > 0) &    (sdss_s82['gPsf'] > 0) &  (sdss_s82['rPsf'] > 0) &  (sdss_s82['iPsf'] > 0) &  (sdss_s82['zPsf'] > 0)       )
   sdss_s82 = sdss_s82[good]


# match SDSS and params using RA and Dec
   h = esutil.htm.HTM(10)
   m1, m2, d12 = h.match(params_trainingset['ra'], params_trainingset['dec'], sdss_s82['sdss_ra'],sdss_s82['sdss_dec'], 1.0/3600., maxmatch=1)

# prepare arrays for SDSS magnitudes and its uncertainty
   s82_sdss_mag_u = np.zeros(params_trainingset.size) - 9999.99
   s82_sdss_mag_g = np.zeros(params_trainingset.size) - 9999.99
   s82_sdss_mag_r = np.zeros(params_trainingset.size) - 9999.99
   s82_sdss_mag_i = np.zeros(params_trainingset.size) - 9999.99
   s82_sdss_mag_z = np.zeros(params_trainingset.size) - 9999.99
   s82_sdss_magerr_u = np.zeros(params_trainingset.size) - 9999.99
   s82_sdss_magerr_g = np.zeros(params_trainingset.size) - 9999.99
   s82_sdss_magerr_r = np.zeros(params_trainingset.size) - 9999.99
   s82_sdss_magerr_i = np.zeros(params_trainingset.size) - 9999.99
   s82_sdss_magerr_z = np.zeros(params_trainingset.size) - 9999.99

# fill in the values where they exist

   s82_sdss_mag_u[m1] = sdss_s82['uPsf'][m2]
   s82_sdss_mag_g[m1] = sdss_s82['gPsf'][m2]
   s82_sdss_mag_r[m1] = sdss_s82['rPsf'][m2]
   s82_sdss_mag_i[m1] = sdss_s82['iPsf'][m2]
   s82_sdss_mag_z[m1] = sdss_s82['zPsf'][m2]


   s82_sdss_magerr_u[m1] =   sdss_s82['uPsfErr'][m2]    
   s82_sdss_magerr_g[m1] =   sdss_s82['gPsfErr'][m2]  
   s82_sdss_magerr_r[m1] =   sdss_s82['rPsfErr'][m2]  
   s82_sdss_magerr_i[m1] =   sdss_s82['iPsfErr'][m2]  
   s82_sdss_magerr_z[m1] =   sdss_s82['zPsfErr'][m2]  
   
  
# galactic coordinates
   gl_trainingset, gb_trainingset = eq2gal(params_trainingset['ra'], params_trainingset['dec'])
   
   EBV_trainingset = hp.get_interp_val(dustmap['ebv'], (90.-gb_trainingset)*np.pi/180., gl_trainingset*np.pi/180.)
   bad = (~np.isfinite(EBV_trainingset)) | (EBV_trainingset <= 0)
   EBV_trainingset[bad] = getval(gl_trainingset[bad], gb_trainingset[bad])

   
# dereddening PS1

   params_trainingset = recfunctions.append_fields(params_trainingset, names=('corrmeanmag_g', 'corrmeanmag_r', 'corrmeanmag_i', 'corrmeanmag_z', 'corrmeanmag_y', 'EBV'),
   data=(params_trainingset['errorweightedmeanmag_g'], params_trainingset['errorweightedmeanmag_r'], params_trainingset['errorweightedmeanmag_i'],params_trainingset['errorweightedmeanmag_z'],params_trainingset['errorweightedmeanmag_y'],EBV_trainingset),
   usemask=False, asrecarray=False)


   for band, c in zip(['g', 'r', 'i', 'z', 'y'], [3.172, 2.271, 1.682, 1.322, 1.087]):
       params_trainingset['corrmeanmag_%s' % band] -= c*EBV_trainingset

    
#dereddening SDSS


   params_trainingset = recfunctions.append_fields(params_trainingset, names=('corr_sdss_mag_u', 'corr_sdss_mag_g', 'corr_sdss_mag_r', 'corr_sdss_mag_i', 'corr_sdss_mag_z', \
   'sdss_magerr_u', 'sdss_magerr_g', 'sdss_magerr_r', 'sdss_magerr_i', 'sdss_magerr_z'),
   data=(s82_sdss_mag_u, s82_sdss_mag_g, s82_sdss_mag_r, s82_sdss_mag_i, s82_sdss_mag_z,\
   s82_sdss_magerr_u, s82_sdss_magerr_g, s82_sdss_magerr_r, s82_sdss_magerr_i, s82_sdss_magerr_z),
   usemask=False, asrecarray=False)


   for band, c in zip(['u', 'g', 'r', 'i', 'z'], [4.239, 3.303, 2.285, 1.698, 1.263]):
       params_trainingset['corr_sdss_mag_%s' % band] -= c*EBV_trainingset

   
   params_trainingset = recfunctions.append_fields(params_trainingset, names=('gr', 'ri', 'iz', 'zy', 'gi', 'W12', 'W12Err', 'ugSDSS', 'grSDSS', 'riSDSS', 'izSDSS'),
   data=(params_trainingset['corrmeanmag_g']-params_trainingset['corrmeanmag_r'], params_trainingset['corrmeanmag_r']-params_trainingset['corrmeanmag_i'], 
   params_trainingset['corrmeanmag_i']-params_trainingset['corrmeanmag_z'], 
   params_trainingset['corrmeanmag_z']-params_trainingset['corrmeanmag_y'], params_trainingset['corrmeanmag_g']-params_trainingset['corrmeanmag_i'],
   params_trainingset['w1mpro']-params_trainingset['w2mpro'], np.sqrt(params_trainingset['w1sigmpro']**2+params_trainingset['w2sigmpro']**2),
   params_trainingset['corr_sdss_mag_u']-params_trainingset['corr_sdss_mag_g'], params_trainingset['corr_sdss_mag_g']-params_trainingset['corr_sdss_mag_r'], 
   params_trainingset['corr_sdss_mag_r']-params_trainingset['corr_sdss_mag_i'], params_trainingset['corr_sdss_mag_i']-params_trainingset['corr_sdss_mag_z']),
   usemask=False, asrecarray=False)
    

# training set: reset missing and unreliable WISE data 

   good = ( (params_trainingset['w1mpro'] != 0.0) &  (params_trainingset['w2mpro'] != 0) &  (params_trainingset['w1sigmpro'] >0) & (params_trainingset['w2sigmpro'] >0) &  (params_trainingset['w1sigmpro'] <0.3) & (params_trainingset['w2sigmpro'] <0.3 ))


   for i, (w12mag, w12err, w1mag,w2mag,w1err,w2err) in enumerate(zip(params_trainingset['W12'], params_trainingset['W12Err'], params_trainingset['w1mpro'], params_trainingset['w2mpro'], params_trainingset['w1sigmpro'], params_trainingset['w2sigmpro'])):
     if( (w1mag == 0.0) |  (w2mag == 0) |  (np.isnan(w12mag)) | (np.isnan(w1mag)) |  (np.isnan(w2mag)) | (w1err > 0.3) |   (w2err > 0.3) | (w1err < 0.0) | (w2err < 0.0) |   (w2err > 0.3)):   # (w1err > 0.1) |   (w2err > 0.1) |
  
         params_trainingset['W12'][i]= -9999.99
         params_trainingset['W12Err'][i] = -9999.99
     
         params_trainingset['w1mpro'][i]=-9999.99
         params_trainingset['w1sigmpro'][i]=-9999.99
         params_trainingset['w2mpro'][i]=-9999.99
         params_trainingset['w2sigmpro'][i]=-9999.99
      
# training set: reset missing and unreliable PS1 data 

     

   good = ( (params_trainingset['corrmeanmag_g'] >0.0)  &   (~np.isnan(params_trainingset['corrmeanmag_g'])) & (np.isfinite(params_trainingset['corrmeanmag_g']))  & (params_trainingset['errorweightedmeanmag_g_sigma'] >0.0) & (~np.isnan(params_trainingset['errorweightedmeanmag_g_sigma'])) & (np.isfinite(params_trainingset['errorweightedmeanmag_g_sigma']))   )
   good_ps1_g = params_trainingset[good]

   good = ( (params_trainingset['corrmeanmag_r'] >0.0)  &   (~np.isnan(params_trainingset['corrmeanmag_r'])) & (np.isfinite(params_trainingset['corrmeanmag_r']))  & (params_trainingset['errorweightedmeanmag_r_sigma'] >0.0) & (~np.isnan(params_trainingset['errorweightedmeanmag_r_sigma'])) & (np.isfinite(params_trainingset['errorweightedmeanmag_r_sigma']))   )
   good_ps1_r = params_trainingset[good]

   good = ( (params_trainingset['corrmeanmag_i'] >0.0)  &   (~np.isnan(params_trainingset['corrmeanmag_i'])) & (np.isfinite(params_trainingset['corrmeanmag_i']))  & (params_trainingset['errorweightedmeanmag_i_sigma'] >0.0) & (~np.isnan(params_trainingset['errorweightedmeanmag_i_sigma'])) & (np.isfinite(params_trainingset['errorweightedmeanmag_i_sigma']))   )
   good_ps1_i = params_trainingset[good]

   good = ( (params_trainingset['corrmeanmag_z'] >0.0)  &   (~np.isnan(params_trainingset['corrmeanmag_z'])) & (np.isfinite(params_trainingset['corrmeanmag_z']))  & (params_trainingset['errorweightedmeanmag_z_sigma'] >0.0) & (~np.isnan(params_trainingset['errorweightedmeanmag_z_sigma'])) & (np.isfinite(params_trainingset['errorweightedmeanmag_z_sigma']))   )
   good_ps1_z = params_trainingset[good]

   good = ( (params_trainingset['corrmeanmag_y'] >0.0)  &   (~np.isnan(params_trainingset['corrmeanmag_y'])) & (np.isfinite(params_trainingset['corrmeanmag_y']))  & (params_trainingset['errorweightedmeanmag_y_sigma'] >0.0) & (~np.isnan(params_trainingset['errorweightedmeanmag_y_sigma'])) & (np.isfinite(params_trainingset['errorweightedmeanmag_y_sigma']))   )
   good_ps1_y = params_trainingset[good]

   inputation_g_ps1 = -9999.99 
   inputation_r_ps1 = -9999.99 
   inputation_i_ps1 = -9999.99 
   inputation_z_ps1 = -9999.99 
   inputation_y_ps1 = -9999.99 

   inputation_magerr_g_ps1 = -9999.99 
   inputation_magerr_r_ps1 = -9999.99
   inputation_magerr_i_ps1 = -9999.99 
   inputation_magerr_z_ps1 = -9999.99
   inputation_magerr_y_ps1 = -9999.99 

#   print("inputation g_ps1: %.8f " % (inputation_g_ps1))
#   print("inputation r_ps1: %.8f " % (inputation_r_ps1))
#   print("inputation i_ps1: %.8f " % (inputation_i_ps1))
#   print("inputation z_ps1: %.8f " % (inputation_z_ps1))
#   print("inputation y_ps1: %.8f " % (inputation_y_ps1))

#   print("inputation magerr_g_ps1: %.8f " % (inputation_magerr_g_ps1))
#   print("inputation magerr_r_ps1: %.8f " % (inputation_magerr_r_ps1))
#   print("inputation magerr_i_ps1: %.8f " % (inputation_magerr_i_ps1))
#   print("inputation magerr_z_ps1: %.8f " % (inputation_magerr_z_ps1))
#   print("inputation magerr_y_ps1: %.8f " % (inputation_magerr_y_ps1))

   for i, (corrmeanmag_g,corrmeanmag_r,corrmeanmag_i,corrmeanmag_z,corrmeanmag_y,errorweightedmeanmag_g_sigma,errorweightedmeanmag_r_sigma,errorweightedmeanmag_i_sigma,errorweightedmeanmag_z_sigma,errorweightedmeanmag_y_sigma) in enumerate(zip(params_trainingset['corrmeanmag_g'], params_trainingset['corrmeanmag_r'],params_trainingset['corrmeanmag_i'],params_trainingset['corrmeanmag_z'],params_trainingset['corrmeanmag_y'], params_trainingset['errorweightedmeanmag_g_sigma'], params_trainingset['errorweightedmeanmag_r_sigma'],params_trainingset['errorweightedmeanmag_i_sigma'],params_trainingset['errorweightedmeanmag_z_sigma'],params_trainingset['errorweightedmeanmag_y_sigma'])):      
     if( (corrmeanmag_g <= 0.0)  |   (np.isnan(corrmeanmag_g)) | (np.isinf(corrmeanmag_g)) | (errorweightedmeanmag_g_sigma <= 0.0) |(np.isnan(errorweightedmeanmag_g_sigma)) | (np.isinf(errorweightedmeanmag_g_sigma))        ):
         params_trainingset['corrmeanmag_g'][i]= inputation_g_ps1
         params_trainingset['errorweightedmeanmag_g'][i]= inputation_g_ps1
         params_trainingset['gr'][i]=inputation_g_ps1-inputation_r_ps1
         params_trainingset['gi'][i]=inputation_g_ps1-inputation_i_ps1
         params_trainingset['errorweightedmeanmag_g_sigma'][i]= inputation_magerr_g_ps1     
      
     if( (corrmeanmag_r <= 0.0)  |   (np.isnan(corrmeanmag_r)) | (np.isinf(corrmeanmag_r)) | (errorweightedmeanmag_r_sigma <= 0.0) |(np.isnan(errorweightedmeanmag_r_sigma)) | (np.isinf(errorweightedmeanmag_r_sigma))        ):
         params_trainingset['corrmeanmag_r'][i]= inputation_r_ps1
         params_trainingset['errorweightedmeanmag_r'][i]= inputation_r_ps1
         params_trainingset['gr'][i]=inputation_g_ps1-inputation_r_ps1
         params_trainingset['ri'][i]=inputation_r_ps1-inputation_i_ps1
         params_trainingset['errorweightedmeanmag_r_sigma'][i]= inputation_magerr_r_ps1           
      
     if( (corrmeanmag_i <= 0.0)  |   (np.isnan(corrmeanmag_i)) | (np.isinf(corrmeanmag_i)) | (errorweightedmeanmag_i_sigma <= 0.0) |(np.isnan(errorweightedmeanmag_i_sigma)) | (np.isinf(errorweightedmeanmag_i_sigma))        ):
         params_trainingset['corrmeanmag_i'][i]= inputation_i_ps1
         params_trainingset['errorweightedmeanmag_i'][i]= inputation_i_ps1
         params_trainingset['ri'][i]=inputation_r_ps1-inputation_i_ps1
         params_trainingset['iz'][i]=inputation_i_ps1-inputation_z_ps1
         params_trainingset['gi'][i]=inputation_g_ps1-inputation_i_ps1
         params_trainingset['errorweightedmeanmag_i_sigma'][i]= inputation_magerr_i_ps1           
      
     if( (corrmeanmag_z <= 0.0)  |   (np.isnan(corrmeanmag_z)) | (np.isinf(corrmeanmag_z)) | (errorweightedmeanmag_z_sigma <= 0.0) |(np.isnan(errorweightedmeanmag_z_sigma)) | (np.isinf(errorweightedmeanmag_z_sigma))        ):
         params_trainingset['corrmeanmag_z'][i]= inputation_z_ps1
         params_trainingset['errorweightedmeanmag_z'][i]= inputation_z_ps1
         params_trainingset['iz'][i]=inputation_i_ps1-inputation_z_ps1
         params_trainingset['zy'][i]=inputation_z_ps1-inputation_y_ps1
         params_trainingset['errorweightedmeanmag_z_sigma'][i]= inputation_magerr_z_ps1      
            
     if( (corrmeanmag_y <= 0.0)  |   (np.isnan(corrmeanmag_y)) | (np.isinf(corrmeanmag_y)) | (errorweightedmeanmag_y_sigma <= 0.0) |(np.isnan(errorweightedmeanmag_y_sigma)) | (np.isinf(errorweightedmeanmag_y_sigma))        ):
         params_trainingset['corrmeanmag_y'][i]= inputation_y_ps1
         params_trainingset['errorweightedmeanmag_y'][i]= inputation_y_ps1
         params_trainingset['zy'][i]=inputation_z_ps1-inputation_y_ps1
         params_trainingset['errorweightedmeanmag_y_sigma'][i]= inputation_magerr_y_ps1      
      
# training set: reset missing and unreliable SDSS data           



   good = ( (params_trainingset['corr_sdss_mag_u'] >0.0)  &   (~np.isnan(params_trainingset['corr_sdss_mag_u'])) & (np.isfinite(params_trainingset['corr_sdss_mag_u']))  & (params_trainingset['sdss_magerr_u'] >0.0) & (params_trainingset['sdss_magerr_u'] <= 0.2) & (~np.isnan(params_trainingset['sdss_magerr_u'])) & (np.isfinite(params_trainingset['sdss_magerr_u']))   )
   good_sdss_u = params_trainingset[good]

   good = ( (params_trainingset['corr_sdss_mag_g'] >0.0)  &   (~np.isnan(params_trainingset['corr_sdss_mag_g'])) & (np.isfinite(params_trainingset['corr_sdss_mag_g']))  & (params_trainingset['sdss_magerr_g'] >0.0) & (params_trainingset['sdss_magerr_g'] <= 0.2) & (~np.isnan(params_trainingset['sdss_magerr_g'])) & (np.isfinite(params_trainingset['sdss_magerr_g']))   )
   good_sdss_g = params_trainingset[good]

   good = ( (params_trainingset['corr_sdss_mag_r'] >0.0)  &   (~np.isnan(params_trainingset['corr_sdss_mag_r'])) & (np.isfinite(params_trainingset['corr_sdss_mag_r']))  & (params_trainingset['sdss_magerr_r'] >0.0) & (params_trainingset['sdss_magerr_r'] <= 0.2) & (~np.isnan(params_trainingset['sdss_magerr_r'])) & (np.isfinite(params_trainingset['sdss_magerr_r']))   )
   good_sdss_r = params_trainingset[good]

   good = ( (params_trainingset['corr_sdss_mag_i'] >0.0)  &   (~np.isnan(params_trainingset['corr_sdss_mag_i'])) & (np.isfinite(params_trainingset['corr_sdss_mag_i']))  & (params_trainingset['sdss_magerr_i'] >0.0) & (params_trainingset['sdss_magerr_i'] <= 0.2) & (~np.isnan(params_trainingset['sdss_magerr_i'])) & (np.isfinite(params_trainingset['sdss_magerr_i']))   )
   good_sdss_i = params_trainingset[good]

   good = ( (params_trainingset['corr_sdss_mag_z'] >0.0)  &   (~np.isnan(params_trainingset['corr_sdss_mag_z'])) & (np.isfinite(params_trainingset['corr_sdss_mag_z']))  & (params_trainingset['sdss_magerr_z'] >0.0) & (params_trainingset['sdss_magerr_z'] <= 0.2) & (~np.isnan(params_trainingset['sdss_magerr_z'])) & (np.isfinite(params_trainingset['sdss_magerr_z']))   )
   good_sdss_z = params_trainingset[good]     

   for i, (corr_sdss_mag_u,corr_sdss_mag_g,corr_sdss_mag_r,corr_sdss_mag_i,corr_sdss_mag_z, sdss_magerr_u,sdss_magerr_g,sdss_magerr_r,sdss_magerr_i,sdss_magerr_z) in enumerate(zip(params_trainingset['corr_sdss_mag_u'], params_trainingset['corr_sdss_mag_g'],params_trainingset['corr_sdss_mag_r'],params_trainingset['corr_sdss_mag_i'],params_trainingset['corr_sdss_mag_z'], params_trainingset['sdss_magerr_u'], params_trainingset['sdss_magerr_g'],params_trainingset['sdss_magerr_r'],params_trainingset['sdss_magerr_i'],params_trainingset['sdss_magerr_z'])):      
     if( (corr_sdss_mag_u <= 0.0)  |   (np.isnan(corr_sdss_mag_u)) | (np.isinf(corr_sdss_mag_u)) | (sdss_magerr_u <= 0.0) | (sdss_magerr_u > 0.2) |(np.isnan(sdss_magerr_u)) | (np.isinf(sdss_magerr_u))        ):
         params_trainingset['corr_sdss_mag_u'][i]= -9999.99
         s82_sdss_mag_u[i]= -9999.99
         params_trainingset['ugSDSS'][i]=-9999.99
         params_trainingset['sdss_magerr_u'][i]= -9999.99  
     if( (corr_sdss_mag_g <= 0.0)  |   (np.isnan(corr_sdss_mag_g)) | (np.isinf(corr_sdss_mag_g))  |(sdss_magerr_g <= 0.0) | (sdss_magerr_g > 0.2) |(np.isnan(sdss_magerr_g)) | (np.isinf(sdss_magerr_g))        ):
         params_trainingset['corr_sdss_mag_g'][i]= -9999.99
         s82_sdss_mag_g[i]= -9999.99
         params_trainingset['ugSDSS'][i]=-9999.99
         params_trainingset['grSDSS'][i]=-9999.99
         params_trainingset['sdss_magerr_g'][i]= -9999.99
     if( (corr_sdss_mag_r <= 0.0)  |   (np.isnan(corr_sdss_mag_r)) | (np.isinf(corr_sdss_mag_r)) |(sdss_magerr_r <= 0.0) | (sdss_magerr_r > 0.2) |(np.isnan(sdss_magerr_r)) | (np.isinf(sdss_magerr_r))         ):
         params_trainingset['corr_sdss_mag_r'][i]= -9999.99
         s82_sdss_mag_r[i]= -9999.99
         params_trainingset['grSDSS'][i]=-9999.99
         params_trainingset['riSDSS'][i]=-9999.99
         params_trainingset['sdss_magerr_r'][i]= -9999.99
     if( (corr_sdss_mag_i <= 0.0)  |   (np.isnan(corr_sdss_mag_i)) | (np.isinf(corr_sdss_mag_i))   |(sdss_magerr_i <= 0.0) | (sdss_magerr_i > 0.2) |(np.isnan(sdss_magerr_i)) | (np.isinf(sdss_magerr_i))       ):
         params_trainingset['corr_sdss_mag_i'][i]= -9999.99
         s82_sdss_mag_i[i]= -9999.99
         params_trainingset['riSDSS'][i]=-9999.99
         params_trainingset['izSDSS'][i]=-9999.99
         params_trainingset['sdss_magerr_i'][i]= -9999.99
     if( (corr_sdss_mag_z <= 0.0)  |   (np.isnan(corr_sdss_mag_z)) | (np.isinf(corr_sdss_mag_z)) |(sdss_magerr_z <= 0.0) | (sdss_magerr_z > 0.2) |(np.isnan(sdss_magerr_z)) | (np.isinf(sdss_magerr_z))         ):
         params_trainingset['corr_sdss_mag_z'][i]= -9999.99
         s82_sdss_mag_z[i]= -9999.99
         params_trainingset['izSDSS'][i]=-9999.99
         params_trainingset['sdss_magerr_z'][i]= -9999.99

      
      
   params_trainingset = recfunctions.append_fields(params_trainingset, names=('iW1'),
   data=(params_trainingset['corrmeanmag_i']-params_trainingset['w1mpro']), usemask=False, asrecarray=False)




   label = np.zeros(params_trainingset.size, dtype='u2')


   label[params_trainingset['SDSStype'] == 'QSO'] = 1   
   label[params_trainingset['SDSStype'] == 'RRLyrae'] = 2






#re-sample within colors and its errors: W1, W2, PS1, SDSS

#10 samples from each object in the training set (additionally to original one)

#errors are not longer parameters

   number_copies=5 


# this should have the the same number of columns as params, and number_copies times the rows of the training set
#newtrainingset = np.zeros(( ((np.floor(params.size/2.))+1)*number_copies, len(par_keys)))
   newtrainingset = np.zeros(( (len(params_trainingset))*(number_copies+1), len(par_keys)))

   print("len(params_trainingset) %d " % (len(params_trainingset)))

   print("len(newtrainingset) %d " % (len(newtrainingset)))

   print("len(par_keys) %d " % (len(par_keys)))


   row=0
   for i in range(0,len(params_trainingset)):
     for j, key in enumerate(par_keys):
          newtrainingset[row,j]=params_trainingset[i][key]
     row=row+1
  
     for k in range(0,number_copies):
    
       #E(B-V)_{sample} drawn from G (E(B-V)_{catalog}, delta E(B-V)=0.1 E(B-V)_{catalog})        
       EBV_sample = numpy.random.normal(params_trainingset[i]['EBV'],0.1*params_trainingset[i]['EBV'])
    
       # 5 % chance that E(B-V)=0, irrespective of catalog entry
       randnumber=randint(1, 20)
    
       if(randnumber ==1):
         EBV_sample=0
      
    #sample new not dereddened grizy for PS1, ugriz for SDSS

       if( (params_trainingset[i]['errorweightedmeanmag_g']<0) | (params_trainingset[i]['errorweightedmeanmag_g_sigma']<0)):
          dereddened_sample_g_ps1 = -9999.99
       else:   
          sample_g_ps1 = numpy.random.normal(params_trainingset[i]['errorweightedmeanmag_g'],params_trainingset[i]['errorweightedmeanmag_g_sigma'])
          #deredden by E(B-V)_{sample} 
          dereddened_sample_g_ps1 = sample_g_ps1 - 3.172*EBV_sample
    
       if( (params_trainingset[i]['errorweightedmeanmag_r']<0) | (params_trainingset[i]['errorweightedmeanmag_r_sigma']<0)):
          dereddened_sample_r_ps1 = -9999.99
       else:   
          sample_r_ps1 = numpy.random.normal(params_trainingset[i]['errorweightedmeanmag_r'],params_trainingset[i]['errorweightedmeanmag_r_sigma'])
          dereddened_sample_r_ps1 = sample_r_ps1 - 2.271*EBV_sample
       if( (params_trainingset[i]['errorweightedmeanmag_i']<0) | (params_trainingset[i]['errorweightedmeanmag_i_sigma']<0)):
          dereddened_sample_i_ps1 = -9999.99
       else:   
          sample_i_ps1 = numpy.random.normal(params_trainingset[i]['errorweightedmeanmag_i'],params_trainingset[i]['errorweightedmeanmag_i_sigma'])    
          dereddened_sample_i_ps1 = sample_i_ps1 - 1.682*EBV_sample
       if( (params_trainingset[i]['errorweightedmeanmag_z']<0) | (params_trainingset[i]['errorweightedmeanmag_z_sigma']<0)):
          dereddened_sample_z_ps1 = -9999.99
       else:   
          sample_z_ps1 = numpy.random.normal(params_trainingset[i]['errorweightedmeanmag_z'],params_trainingset[i]['errorweightedmeanmag_z_sigma'])
          dereddened_sample_z_ps1 = sample_z_ps1 - 1.322*EBV_sample
        
       if( (params_trainingset[i]['errorweightedmeanmag_y']<0) | (params_trainingset[i]['errorweightedmeanmag_y_sigma']<0)):
          dereddened_sample_y_ps1 = -9999.99
       else:   
          sample_y_ps1 = numpy.random.normal(params_trainingset[i]['errorweightedmeanmag_y'],params_trainingset[i]['errorweightedmeanmag_y_sigma'])    
          dereddened_sample_y_ps1 = sample_y_ps1 - 1.087*EBV_sample 
    
    

       if( (s82_sdss_mag_u[i]<0) | (params_trainingset[i]['sdss_magerr_u']<0)):
          dereddened_sample_u_sdss = -9999.99
       else:   
          sample_u_sdss = numpy.random.normal(s82_sdss_mag_u[i],params_trainingset[i]['sdss_magerr_u'])
          dereddened_sample_u_sdss = sample_u_sdss - 4.239*EBV_sample  
       if( (s82_sdss_mag_g[i]<0) | (params_trainingset[i]['sdss_magerr_g']<0)):
          dereddened_sample_g_sdss = -9999.99
       else:   
          sample_g_sdss = numpy.random.normal(s82_sdss_mag_g[i],params_trainingset[i]['sdss_magerr_g'])
          dereddened_sample_g_sdss = sample_g_sdss - 3.303*EBV_sample 
       if( (s82_sdss_mag_r[i]<0) | (params_trainingset[i]['sdss_magerr_r']<0)):
          dereddened_sample_r_sdss = -9999.99
       else:   
          sample_r_sdss = numpy.random.normal(s82_sdss_mag_r[i],params_trainingset[i]['sdss_magerr_r'])  
          dereddened_sample_r_sdss = sample_r_sdss - 2.285*EBV_sample
       
       if( (s82_sdss_mag_i[i]<0) | (params_trainingset[i]['sdss_magerr_i']<0)):
          dereddened_sample_i_sdss = -9999.99
       else:   
          sample_i_sdss = numpy.random.normal(s82_sdss_mag_i[i],params_trainingset[i]['sdss_magerr_i']) 
          dereddened_sample_i_sdss = sample_i_sdss - 1.698*EBV_sample   
       if( (s82_sdss_mag_z[i]<0) | (params_trainingset[i]['sdss_magerr_z']<0)):
          dereddened_sample_z_sdss = -9999.99
       else:   
          sample_z_sdss = numpy.random.normal(s82_sdss_mag_z[i],params_trainingset[i]['sdss_magerr_z'])       
          dereddened_sample_z_sdss = sample_z_sdss - 1.263*EBV_sample  
       
     
    
       #bright it so that m_r after dereddening by E(B-V)_{sample} is the same as after dereddening by E(B-V)_{catalog}
       if((dereddened_sample_r_ps1 <0 )|(corrmeanmag_r<0)):
    
           brightened_sample_g_ps1 = -9999.99
           brightened_sample_r_ps1 = -9999.99
           brightened_sample_i_ps1 = -9999.99
           brightened_sample_z_ps1 = -9999.99
           brightened_sample_y_ps1 = -9999.99
        
       else:
           bright_factor_ps1 = dereddened_sample_r_ps1 / corrmeanmag_r
           #print("bright_factor_ps1 %.8f " % (bright_factor_ps1))

           brightened_sample_g_ps1 = bright_factor_ps1*dereddened_sample_g_ps1
           brightened_sample_r_ps1 = bright_factor_ps1*dereddened_sample_r_ps1
           brightened_sample_i_ps1 = bright_factor_ps1*dereddened_sample_i_ps1
           brightened_sample_z_ps1 = bright_factor_ps1*dereddened_sample_z_ps1
           brightened_sample_y_ps1 = bright_factor_ps1*dereddened_sample_y_ps1
  
 
       if((dereddened_sample_r_sdss <0 )|(corr_sdss_mag_r<0)):
    
           brightened_sample_u_sdss = -9999.99
           brightened_sample_g_sdss = -9999.99
           brightened_sample_r_sdss = -9999.99
           brightened_sample_i_sdss = -9999.99
           brightened_sample_z_sdss = -9999.99
        
       else:
           bright_factor_sdss = dereddened_sample_r_sdss / corr_sdss_mag_r
           #print("bright_factor_sdss %.8f " % (bright_factor_sdss))
    
           brightened_sample_u_sdss = bright_factor_sdss*dereddened_sample_u_sdss
           brightened_sample_g_sdss = bright_factor_sdss*dereddened_sample_g_sdss
           brightened_sample_r_sdss = bright_factor_sdss*dereddened_sample_r_sdss
           brightened_sample_i_sdss = bright_factor_sdss*dereddened_sample_i_sdss
           brightened_sample_z_sdss = bright_factor_sdss*dereddened_sample_z_sdss
    
    
    #assign a slightly wrong E(B-V): I had used E(B-V)_{catalog}, this was assigned by        else:   sample = tablevalue

    
    #and for WISE
    
    #print("w1mpro %.8f  w1sigmpro %.8f" % (params_trainingset[i]['w1mpro'],params_trainingset[i]['w1sigmpro']))
    #print("W12 %.8f  W12Err %.8f" % (params_trainingset[i]['W12'],params_trainingset[i]['W12Err']))
    
       if( (params_trainingset[i]['w1mpro']<0) | (params_trainingset[i]['w1sigmpro']<0)):
          sample_w1 = -9999.99
       else:        
          sample_w1 =numpy.random.normal(params_trainingset[i]['w1mpro'],params_trainingset[i]['w1sigmpro'])
    
       if( (params_trainingset[i]['W12']<0) | (params_trainingset[i]['W12Err']<0)):
          sample_w12 = -9999.99
       else:     
          sample_w12 = numpy.random.normal(params_trainingset[i]['W12'],params_trainingset[i]['W12Err'])
    
    
            
       for j, key in enumerate(par_keys):
	
          tablevalue = params_trainingset[i][key]
              #the ones for which we sample, 'gr', 'ri', 'iz', 'zy', 'ugSDSS', 'grSDSS', 'riSDSS', 'izSDSS', 'W12', 'iW1', 'corrmeanmag_r'
          if( (key=='gr') | (key=='ri') | (key=='iz') | (key=='zy') | (key=='ugSDSS') | (key=='grSDSS') | (key=='riSDSS') | (key=='izSDSS') | (key=='W12') | (key=='iW1') | (key=='corrmeanmag_r')):
             if (key=='gr'):
            
               #make the new colors
               sample = brightened_sample_g_ps1-brightened_sample_r_ps1
               if( (brightened_sample_g_ps1<0) | (brightened_sample_r_ps1<0)):
	         sample = -9999.99
             if (key=='ri'):
               sample = brightened_sample_r_ps1-brightened_sample_i_ps1
               if( (brightened_sample_r_ps1<0) | (brightened_sample_i_ps1<0)):
	         sample = -9999.99
             if (key=='iz'):
               sample = brightened_sample_i_ps1-brightened_sample_z_ps1
               if( (brightened_sample_i_ps1<0) | (brightened_sample_z_ps1<0)):
	         sample = -9999.99
             if (key=='zy'):
               sample = brightened_sample_z_ps1-brightened_sample_y_ps1
               if( (brightened_sample_z_ps1<0) | (brightened_sample_y_ps1<0)):
	         sample = -9999.99
             if (key=='ugSDSS'):
               sample = brightened_sample_u_sdss-brightened_sample_g_sdss
               if( (brightened_sample_u_sdss<0) | (brightened_sample_g_sdss<0)):
	         sample = -9999.99
             if (key=='grSDSS'):
               sample = brightened_sample_g_sdss-brightened_sample_r_sdss
               if( (brightened_sample_g_sdss<0) | (brightened_sample_r_sdss<0)):
	         sample = -9999.99
             if (key=='riSDSS'):
               sample = brightened_sample_r_sdss-brightened_sample_i_sdss
               if( (brightened_sample_r_sdss<0) | (brightened_sample_i_sdss<0)):
	         sample = -9999.99
             if (key=='izSDSS'):
               sample = brightened_sample_i_sdss-brightened_sample_z_sdss
               if( (brightened_sample_i_sdss<0) | (brightened_sample_z_sdss<0)):
	         sample = -9999.99
             if (key=='W12'):
               sample = sample_w12
             if (key=='iW1'):
               sample = brightened_sample_i_ps1-sample_w1
             if( (brightened_sample_i_ps1<0) | (sample_w1<0)):
	         sample = -9999.99
             if (key=='corrmeanmag_r'):
               sample = brightened_sample_r_ps1
            
          #print("key: %s   tablevalue %.8f    tablevalue_error %.8f " % (key, tablevalue, tablevalue_error))
          
          else:
             sample = tablevalue
          newtrainingset[row,j]=sample
       row=row+1
    
   
  

   print("params_trainingset is done")
 
   cPickle.dump(newtrainingset, open('newtrainingset.pickle', 'wb'), protocol=2)
   
   print("newtrainingset.pickle")
#also, copy labels
   newlabel = np.zeros(len(params_trainingset)*(number_copies+1), dtype='u2')

   row=0
   for i in range(0,len(params_trainingset)):
  
     newlabel[row]=label[i]
 
     row=row+1
  
     for k in range(0,number_copies):
  
      newlabel[row]=label[i]
      row=row+1

   print("newlabel is done")


   cPickle.dump(newlabel, open('newlabel.pickle', 'wb'), protocol=2)
   print("newlabel.pickle")
   newtrainingset = cPickle.load(open('newtrainingset.pickle', 'rb'))
   newlabel = cPickle.load(open('newlabel.pickle', 'rb'))
   
   Y = np.where(newlabel == 2, 1, 0)    
   
   from sklearn.ensemble import RandomForestClassifier
   classifer = RandomForestClassifier(n_estimators=500, n_jobs=-1,class_weight='subsample')
   rf_qso = classifer.fit(newtrainingset, Y)

   print("train_rf is done")


   #serialize the trained Random Forest Classifier
   cPickle.dump(rf_qso, open('/a41233d1/hernitschek/_results_juelich/_results_classifier/_pickled_RFCs/RFC_RRLyrae_reset-9999_noEBV_noSDSS_classweight.pickle', 'wb'), protocol=2)
   print("RF serialized")
   
   
#####################

#this is if newtrainingset and newlabel where done and serialized correctly, but it had crashed before rf_qso was trained and serialized
#if(repeat_train==1):   
#
#   print("repeat_train=1")
#   newtrainingset = cPickle.load(open('newtrainingset.pickle', 'rb'))
#
#   newlabel = cPickle.load(open('newlabel.pickle', 'rb'))
#   
#   Y = np.where(newlabel == 2, 1, 0)    
#   rf_qso = learn_class_csc.train_rf(newtrainingset, Y)
#
#   print("train_rf is done")#
#
#
#  #serialize the trained Random Forest Classifier
#   cPickle.dump(rf_qso, open('RFC_RRLyrae_reset-9999.pickle', 'wb'), protocol=2)
#   print("RF serialized")


#load the trained Random Forest Classifier

if(make_trainingset==0):
  
  rf_qso = cPickle.load(open('/a41233d1/hernitschek/_results_juelich/_results_classifier/_pickled_RFCs/RFC_RRLyrae_reset-9999_noEBV_noSDSS_classweight.pickle', 'rb')) #This function automatically determines whether the data stream was written in binary mode or not.
  print("RF loaded")

   

#FOR ALL THE TARGET SET FILES

#####################



#FOR ALL THE TARGET SET FILES



#sdssdata/sdssdata_l65.00to70.00_b-30.00to-25.00_cell13552457182959632383.txt
#drwdata/
#maxloglik_params_drw_lightcurves_l65.00to70.00_b-30.00to-25.00_cell13552457182959632383_det10_1.fits.txt
#maxloglik_params_drw_lightcurves_l65.00to70.00_b-30.00to-25.00_cell13552457182959632383_det10.fits.txt



#for all files in drw
#  get filename
# prepare filename for loading SDSS 
# in the end: write to a file that has corresponding name 

  folderList = os.listdir(indirectory)

  print("%i " % (len(folderList)))
  for filenum in range(0,len(folderList)):

   print("%s " % (folderList[filenum]))
   filename_drw = folderList[filenum]
   namstring=filename_drw[20:]

    #search position of "d" in namstring
   pos=namstring.find("d")
  
   filename_sdss = "sdssdata"+namstring[12:pos-1]+".txt"
   filename_result = "result_"+filename_drw
   print("%s " % (filename_drw))
   print("%s " % (filename_sdss))
   print("%s " % (filename_result))  
  
  
    #load DRW output
   try:
       params_classify = np.genfromtxt(indirectory+filename_drw, names = 'objid, ra, dec, longitude_l, latitude_b, maxloglik, q, stellar_locus, omega_best, tau_best, errorweightedmeanmag_g, errorweightedmeanmag_r, errorweightedmeanmag_i, errorweightedmeanmag_z, errorweightedmeanmag_y, strucfuncmeanmag_g, strucfuncmeanmag_r, strucfuncmeanmag_i, strucfuncmeanmag_z, strucfuncmeanmag_y, gbandepochs_cleaned, rbandepochs_cleaned, ibandepochs_cleaned, zbandepochs_cleaned, ybandepochs_cleaned, w1mpro, w2mpro, w1sigmpro, w2sigmpro, wisera, wisedec, sigrawise, sigdecwise', dtype='|S20, f8, f8, f8, f8, f8, f8, u1, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, u2, u2, u2, u2, u2, f8, f8, f8, f8, f8, f8, f8, f8', skip_header=1)
   except ValueError:
       print 'input file damaged'
   else:
      
    already_proc=0
    #if (os.path.isfile(outdirectory+filename_result)):
#	    already_proc=1
 
 
    #if not in list: was not fully processed --> redo!
    
    
    
    if (os.path.isfile(outdirectory+filename_result)):
	    already_proc=1
     
      
    
    if( (params_classify.size>1) & (already_proc==0)):
   
   
            # prepare arrays for SDSS magnitudes and its uncertainty
      targetset_sdss_mag_u = np.zeros(params_classify.size) - 9999.99
      targetset_sdss_mag_g = np.zeros(params_classify.size) - 9999.99
      targetset_sdss_mag_r = np.zeros(params_classify.size) - 9999.99
      targetset_sdss_mag_i = np.zeros(params_classify.size) - 9999.99
      targetset_sdss_mag_z = np.zeros(params_classify.size) - 9999.99
      targetset_sdss_magerr_u = np.zeros(params_classify.size) - 9999.99
      targetset_sdss_magerr_g = np.zeros(params_classify.size) - 9999.99
      targetset_sdss_magerr_r = np.zeros(params_classify.size) - 9999.99
      targetset_sdss_magerr_i = np.zeros(params_classify.size) - 9999.99
      targetset_sdss_magerr_z = np.zeros(params_classify.size) - 9999.99

#    print("params_classify.size: %i" % (params_classify.size))
    
      SDSSavailable=0

    #if sdss file exists and is larger than 0: 
      if (os.path.isfile(sdssdirectory+filename_sdss)):
        if (os.path.getsize(sdssdirectory+filename_sdss)>0)        :    #only 1 line, means only header available, no data
           with open(sdssdirectory+filename_sdss) as f:
            line_count = len(f.readlines())
            if(line_count>1):
               SDSSavailable=1
  
      print("SDSS file available: %i" % (SDSSavailable))
      if(SDSSavailable==1):

      
    #load SDSS 
        sdss_patch = np.genfromtxt(sdssdirectory+filename_sdss, \
        names = 'objid,ps1_ra,ps1_dec,sdss_ra,sdss_dec,uPsf,gPsf,rPsf,iPsf,zPsf,uPsfErr,gPsfErr,rPsfErr,iPsfErr,zPsfErr', \
        dtype='|S20, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8', skip_header=1)
        
       
	good = (     (sdss_patch['uPsf'] > 0) &    (sdss_patch['gPsf'] > 0) &  (sdss_patch['rPsf'] > 0) &  (sdss_patch['iPsf'] > 0) &  (sdss_patch['zPsf'] > 0)       )
	sdss_patch = sdss_patch[good]


        print 'reading SDSS successfully'
  
  
  # for patch: match SDSS and params using RA and Dec
        h = esutil.htm.HTM(10)
        m1, m2, d12 = h.match(params_classify['ra'], params_classify['dec'], sdss_patch['sdss_ra'],sdss_patch['sdss_dec'], 1.0/3600., maxmatch=1)

        print("merge done")
	#print("len m1: %i" % (len(m1)))
        #print("len m2: %i" % (len(m2)))
	if( (len(m1)==0) | (len(m2)==0) ): #merge is empty
	  print("m1 or m2 empty")
	  SDSSavailable=0
  
      print("SDSS data available: %i" % (SDSSavailable))
      if(SDSSavailable==1):

       
        # fill in the values where they exist


          targetset_sdss_mag_u[m1] = sdss_patch['uPsf'][m2]
          targetset_sdss_mag_g[m1] = sdss_patch['gPsf'][m2]
          targetset_sdss_mag_r[m1] = sdss_patch['rPsf'][m2]
          targetset_sdss_mag_i[m1] = sdss_patch['iPsf'][m2]
          targetset_sdss_mag_z[m1] = sdss_patch['zPsf'][m2]

          targetset_sdss_magerr_u[m1] =   sdss_patch['uPsfErr'][m2]    
          targetset_sdss_magerr_g[m1] =   sdss_patch['gPsfErr'][m2]  
          targetset_sdss_magerr_r[m1] =   sdss_patch['rPsfErr'][m2]  
          targetset_sdss_magerr_i[m1] =   sdss_patch['iPsfErr'][m2]  
          targetset_sdss_magerr_z[m1] =   sdss_patch['zPsfErr'][m2]  



  # galactic coordinates
      gl_classify, gb_classify = eq2gal(params_classify['ra'], params_classify['dec'])

  # Eddie's dust map

      EBV_classify = hp.get_interp_val(dustmap['ebv'], (90.-gb_classify)*np.pi/180., gl_classify*np.pi/180.)
      bad = (~np.isfinite(EBV_classify)) | (EBV_classify <= 0.0)
      EBV_classify[bad] = getval(gl_classify[bad], gb_classify[bad])

# dereddening PS1


      params_classify = recfunctions.append_fields(params_classify, names=('corrmeanmag_g', 'corrmeanmag_r', 'corrmeanmag_i', 'corrmeanmag_z', 'corrmeanmag_y', 'EBV'),
      data=(params_classify['errorweightedmeanmag_g'], params_classify['errorweightedmeanmag_r'], params_classify['errorweightedmeanmag_i'],params_classify['errorweightedmeanmag_z'],params_classify['errorweightedmeanmag_y'],EBV_classify),
      usemask=False, asrecarray=False)


      for band, c in zip(['g', 'r', 'i', 'z', 'y'], [3.172, 2.271, 1.682, 1.322, 1.087]):
          params_classify['corrmeanmag_%s' % band] -= c*EBV_classify
    
    
    
      params_classify = recfunctions.append_fields(params_classify, names=('corr_sdss_mag_u', 'corr_sdss_mag_g', 'corr_sdss_mag_r', 'corr_sdss_mag_i', 'corr_sdss_mag_z', \
      'sdss_magerr_u', 'sdss_magerr_g', 'sdss_magerr_r', 'sdss_magerr_i', 'sdss_magerr_z'),
      data=(targetset_sdss_mag_u, targetset_sdss_mag_g, targetset_sdss_mag_r, targetset_sdss_mag_i, targetset_sdss_mag_z,\
      targetset_sdss_magerr_u, targetset_sdss_magerr_g, targetset_sdss_magerr_r, targetset_sdss_magerr_i, targetset_sdss_magerr_z),
      usemask=False, asrecarray=False)

#dereddening SDSS

#classify    
    
      if(SDSSavailable==1): 



        for band, c in zip(['u', 'g', 'r', 'i', 'z'], [4.239, 3.303, 2.285, 1.698, 1.263]):
          params_classify['corr_sdss_mag_%s' % band] -= c*EBV_classify

   
        params_classify = recfunctions.append_fields(params_classify, names=('gr', 'ri', 'iz', 'zy', 'gi', 'W12', 'W12Err', 'ugSDSS', 'grSDSS', 'riSDSS', 'izSDSS'),
        data=(params_classify['corrmeanmag_g']-params_classify['corrmeanmag_r'], params_classify['corrmeanmag_r']-params_classify['corrmeanmag_i'], 
        params_classify['corrmeanmag_i']-params_classify['corrmeanmag_z'], 
        params_classify['corrmeanmag_z']-params_classify['corrmeanmag_y'], params_classify['corrmeanmag_g']-params_classify['corrmeanmag_i'],
        params_classify['w1mpro']-params_classify['w2mpro'], np.sqrt(params_classify['w1sigmpro']**2+params_classify['w2sigmpro']**2),
        params_classify['corr_sdss_mag_u']-params_classify['corr_sdss_mag_g'], params_classify['corr_sdss_mag_g']-params_classify['corr_sdss_mag_r'], 
        params_classify['corr_sdss_mag_r']-params_classify['corr_sdss_mag_i'], params_classify['corr_sdss_mag_i']-params_classify['corr_sdss_mag_z']),
        usemask=False, asrecarray=False)
 
      else:   
    #in this case, targetset_sdss_mag_u contains -9999.99
        
        targetset_sdss_mag_u = np.zeros(params_classify.size) - 9999.99
        targetset_sdss_mag_g = np.zeros(params_classify.size) - 9999.99
        targetset_sdss_mag_r = np.zeros(params_classify.size) - 9999.99
        targetset_sdss_mag_i = np.zeros(params_classify.size) - 9999.99
        targetset_sdss_mag_z = np.zeros(params_classify.size) - 9999.99
        targetset_sdss_magerr_u = np.zeros(params_classify.size) - 9999.99
        targetset_sdss_magerr_g = np.zeros(params_classify.size) - 9999.99
        targetset_sdss_magerr_r = np.zeros(params_classify.size) - 9999.99
        targetset_sdss_magerr_i = np.zeros(params_classify.size) - 9999.99
        targetset_sdss_magerr_z = np.zeros(params_classify.size) - 9999.99

        
        params_classify = recfunctions.append_fields(params_classify, names=('gr', 'ri', 'iz', 'zy', 'gi', 'W12', 'W12Err', 'ugSDSS', 'grSDSS', 'riSDSS', 'izSDSS'),
        data=(params_classify['corrmeanmag_g']-params_classify['corrmeanmag_r'], params_classify['corrmeanmag_r']-params_classify['corrmeanmag_i'], 
        params_classify['corrmeanmag_i']-params_classify['corrmeanmag_z'], 
        params_classify['corrmeanmag_z']-params_classify['corrmeanmag_y'], params_classify['corrmeanmag_g']-params_classify['corrmeanmag_i'],
        params_classify['w1mpro']-params_classify['w2mpro'], np.sqrt(params_classify['w1sigmpro']**2+params_classify['w2sigmpro']**2),
        targetset_sdss_mag_u, targetset_sdss_mag_u, 
        targetset_sdss_mag_u, targetset_sdss_mag_u),
        usemask=False, asrecarray=False)
 
 
# target set: reset missing and unreliable WISE data

 
      for i, (w12mag, w12err, w1mag,w2mag,w1err,w2err) in enumerate(zip(params_classify['W12'], params_classify['W12Err'], params_classify['w1mpro'], params_classify['w2mpro'], params_classify['w1sigmpro'], params_classify['w2sigmpro'])):
        if( (w1mag == 0.0) |  (w2mag == 0) |  (np.isnan(w12mag)) | (np.isnan(w1mag)) |  (np.isnan(w2mag)) | (w1err <0) |   (w2err <0) | (w1err > 0.3) |   (w2err > 0.3)):   # (w1err > 0.1) |   (w2err > 0.1) |
            
            params_classify['W12'][i]= -9999.99
            params_classify['W12Err'][i] = -9999.99
      
            params_classify['w1mpro'][i]=-9999.99
            params_classify['w1sigmpro'][i]=-9999.99
            params_classify['w2mpro'][i]=-9999.99
            params_classify['w2sigmpro'][i]=-9999.99
     
# target set: reset missing and unreliable PS1 data
     
      good = ( (params_classify['corrmeanmag_g'] >0.0)  &   (~np.isnan(params_classify['corrmeanmag_g'])) & (np.isfinite(params_classify['corrmeanmag_g']))  )
      good_ps1_g = params_classify[good]

      good = ( (params_classify['corrmeanmag_r'] >0.0)  &   (~np.isnan(params_classify['corrmeanmag_r'])) & (np.isfinite(params_classify['corrmeanmag_r']))   )
      good_ps1_r = params_classify[good]

      good = ( (params_classify['corrmeanmag_i'] >0.0)  &   (~np.isnan(params_classify['corrmeanmag_i'])) & (np.isfinite(params_classify['corrmeanmag_i']))   )
      good_ps1_i = params_classify[good]

      good = ( (params_classify['corrmeanmag_z'] >0.0)  &   (~np.isnan(params_classify['corrmeanmag_z'])) & (np.isfinite(params_classify['corrmeanmag_z']))   )
      good_ps1_z = params_classify[good]

      good = ( (params_classify['corrmeanmag_y'] >0.0)  &   (~np.isnan(params_classify['corrmeanmag_y'])) & (np.isfinite(params_classify['corrmeanmag_y']))  )
      good_ps1_y = params_classify[good]

      inputation_g_ps1 = -9999.99
      inputation_r_ps1 = -9999.99
      inputation_i_ps1 = -9999.99 
      inputation_z_ps1 = -9999.99 
      inputation_y_ps1 = -9999.99 

      for i, (corrmeanmag_g,corrmeanmag_r,corrmeanmag_i,corrmeanmag_z,corrmeanmag_y) in enumerate(zip(params_classify['corrmeanmag_g'], params_classify['corrmeanmag_r'],params_classify['corrmeanmag_i'],params_classify['corrmeanmag_z'],params_classify['corrmeanmag_y'])):      
        if( (corrmeanmag_g <= 0.0)  |   (np.isnan(corrmeanmag_g)) | (np.isinf(corrmeanmag_g))   ):
          params_classify['corrmeanmag_g'][i]= inputation_g_ps1
          params_classify['errorweightedmeanmag_g'][i]= inputation_g_ps1
          params_classify['gr'][i]=inputation_g_ps1-inputation_r_ps1
          params_classify['gi'][i]=inputation_g_ps1-inputation_i_ps1
      
        if( (corrmeanmag_r <= 0.0)  |   (np.isnan(corrmeanmag_r)) | (np.isinf(corrmeanmag_r)) ):
          params_classify['corrmeanmag_r'][i]= inputation_r_ps1
          params_classify['errorweightedmeanmag_r'][i]= inputation_r_ps1
          params_classify['gr'][i]=inputation_g_ps1-inputation_r_ps1
          params_classify['ri'][i]=inputation_r_ps1-inputation_i_ps1
      
        if( (corrmeanmag_i <= 0.0)  |   (np.isnan(corrmeanmag_i)) | (np.isinf(corrmeanmag_i))   ):
          params_classify['corrmeanmag_i'][i]= inputation_i_ps1
          params_classify['errorweightedmeanmag_i'][i]= inputation_i_ps1
          params_classify['ri'][i]=inputation_r_ps1-inputation_i_ps1
          params_classify['iz'][i]=inputation_i_ps1-inputation_z_ps1
          params_classify['gi'][i]=inputation_g_ps1-inputation_i_ps1

        if( (corrmeanmag_z <= 0.0)  |   (np.isnan(corrmeanmag_z)) | (np.isinf(corrmeanmag_z))     ):
          params_classify['corrmeanmag_z'][i]= inputation_z_ps1
          params_classify['errorweightedmeanmag_z'][i]= inputation_z_ps1
          params_classify['iz'][i]=inputation_i_ps1-inputation_z_ps1
          params_classify['zy'][i]=inputation_z_ps1-inputation_y_ps1
            
        if( (corrmeanmag_y <= 0.0)  |   (np.isnan(corrmeanmag_y)) | (np.isinf(corrmeanmag_y))      ):
	  params_classify['errorweightedmeanmag_y'][i]= inputation_y_ps1
          params_classify['corrmeanmag_y'][i]= inputation_y_ps1
          params_classify['zy'][i]=inputation_z_ps1-inputation_y_ps1
      
      
      
# target set: reset missing and unreliable SDSS data     

      if(SDSSavailable==1): 
  
        good = ( (params_classify['corr_sdss_mag_u'] >0.0)  &   (~np.isnan(params_classify['corr_sdss_mag_u'])) & (np.isfinite(params_classify['corr_sdss_mag_u']))  & (params_classify['sdss_magerr_u'] >0.0) & (params_classify['sdss_magerr_u'] <= 0.2) & (~np.isnan(params_classify['sdss_magerr_u'])) & (np.isfinite(params_classify['sdss_magerr_u']))   )
        good_sdss_u = params_classify[good]

        good = ( (params_classify['corr_sdss_mag_g'] >0.0)  &   (~np.isnan(params_classify['corr_sdss_mag_g'])) & (np.isfinite(params_classify['corr_sdss_mag_g']))  & (params_classify['sdss_magerr_g'] >0.0) & (params_classify['sdss_magerr_g'] <= 0.2) & (~np.isnan(params_classify['sdss_magerr_g'])) & (np.isfinite(params_classify['sdss_magerr_g']))   )
        good_sdss_g = params_classify[good]

        good = ( (params_classify['corr_sdss_mag_r'] >0.0)  &   (~np.isnan(params_classify['corr_sdss_mag_r'])) & (np.isfinite(params_classify['corr_sdss_mag_r']))  & (params_classify['sdss_magerr_r'] >0.0) & (params_classify['sdss_magerr_r'] <= 0.2) & (~np.isnan(params_classify['sdss_magerr_r'])) & (np.isfinite(params_classify['sdss_magerr_r']))   )
        good_sdss_r = params_classify[good]

        good = ( (params_classify['corr_sdss_mag_i'] >0.0)  &   (~np.isnan(params_classify['corr_sdss_mag_i'])) & (np.isfinite(params_classify['corr_sdss_mag_i']))  & (params_classify['sdss_magerr_i'] >0.0) & (params_classify['sdss_magerr_i'] <= 0.2) & (~np.isnan(params_classify['sdss_magerr_i'])) & (np.isfinite(params_classify['sdss_magerr_i']))   )
        good_sdss_i = params_classify[good]

        good = ( (params_classify['corr_sdss_mag_z'] >0.0)  &   (~np.isnan(params_classify['corr_sdss_mag_z'])) & (np.isfinite(params_classify['corr_sdss_mag_z']))  & (params_classify['sdss_magerr_z'] >0.0) & (params_classify['sdss_magerr_z'] <= 0.2) & (~np.isnan(params_classify['sdss_magerr_z'])) & (np.isfinite(params_classify['sdss_magerr_z']))   )
        good_sdss_z = params_classify[good]     

           
      
        inputation_u = -9999.99 
        inputation_g = -9999.99 
        inputation_r = -9999.99
        inputation_i = -9999.99 
        inputation_z = -9999.99 

        inputation_magerr_u = -9999.99 
        inputation_magerr_g = -9999.99 
        inputation_magerr_r = -9999.99
        inputation_magerr_i = -9999.99 
        inputation_magerr_z = -9999.99 


        for i, (corr_sdss_mag_u,corr_sdss_mag_g,corr_sdss_mag_r,corr_sdss_mag_i,corr_sdss_mag_z, sdss_magerr_u,sdss_magerr_g,sdss_magerr_r,sdss_magerr_i,sdss_magerr_z) in enumerate(zip(params_classify['corr_sdss_mag_u'], params_classify['corr_sdss_mag_g'],params_classify['corr_sdss_mag_r'],params_classify['corr_sdss_mag_i'],params_classify['corr_sdss_mag_z'], params_classify['sdss_magerr_u'], params_classify['sdss_magerr_g'],params_classify['sdss_magerr_r'],params_classify['sdss_magerr_i'],params_classify['sdss_magerr_z'])):      
          if( (corr_sdss_mag_u <= 0.0)  |   (np.isnan(corr_sdss_mag_u)) | (np.isinf(corr_sdss_mag_u)) | (sdss_magerr_u <= 0.0) | (sdss_magerr_u > 0.2) |(np.isnan(sdss_magerr_u)) | (np.isinf(sdss_magerr_u))        ):
            params_classify['corr_sdss_mag_u'][i]=  -9999.99
            targetset_sdss_mag_u[i]=  -9999.99
            params_classify['ugSDSS'][i]= -9999.99
            params_classify['sdss_magerr_u'][i]=  -9999.99      
          if( (corr_sdss_mag_g <= 0.0)  |   (np.isnan(corr_sdss_mag_g)) | (np.isinf(corr_sdss_mag_g))  |(sdss_magerr_g <= 0.0) | (sdss_magerr_g > 0.2) |(np.isnan(sdss_magerr_g)) | (np.isinf(sdss_magerr_g))        ):
            params_classify['corr_sdss_mag_g'][i]=  -9999.99
            targetset_sdss_mag_g[i]=  -9999.99
            params_classify['ugSDSS'][i]= -9999.99
            params_classify['grSDSS'][i]= -9999.99
            params_classify['sdss_magerr_g'][i]=  -9999.99
          if( (corr_sdss_mag_r <= 0.0)  |   (np.isnan(corr_sdss_mag_r)) | (np.isinf(corr_sdss_mag_r)) |(sdss_magerr_r <= 0.0) | (sdss_magerr_r > 0.2) |(np.isnan(sdss_magerr_r)) | (np.isinf(sdss_magerr_r))         ):
            params_classify['corr_sdss_mag_r'][i]=  -9999.99
            targetset_sdss_mag_r[i]=  -9999.99
            params_classify['grSDSS'][i]= -9999.99
            params_classify['riSDSS'][i]= -9999.99
            params_classify['sdss_magerr_r'][i]=  -9999.99
          if( (corr_sdss_mag_i <= 0.0)  |   (np.isnan(corr_sdss_mag_i)) | (np.isinf(corr_sdss_mag_i))   |(sdss_magerr_i <= 0.0) | (sdss_magerr_i > 0.2) |(np.isnan(sdss_magerr_i)) | (np.isinf(sdss_magerr_i))       ):
            params_classify['corr_sdss_mag_i'][i]=  -9999.99
            targetset_sdss_mag_i[i]=  -9999.99
            params_classify['riSDSS'][i]= -9999.99
            params_classify['izSDSS'][i]= -9999.99
            params_classify['sdss_magerr_i'][i]=  -9999.99
          if( (corr_sdss_mag_z <= 0.0)  |   (np.isnan(corr_sdss_mag_z)) | (np.isinf(corr_sdss_mag_z)) |(sdss_magerr_z <= 0.0) | (sdss_magerr_z > 0.2) |(np.isnan(sdss_magerr_z)) | (np.isinf(sdss_magerr_z))         ):
            params_classify['corr_sdss_mag_z'][i]=  -9999.99
            targetset_sdss_mag_z[i]= -9999.99
            params_classify['izSDSS'][i]= -9999.99
            params_classify['sdss_magerr_z'][i]=  -9999.99   


      params_classify = recfunctions.append_fields(params_classify, names=('iW1'),
    
      data=(params_classify['corrmeanmag_i']-params_classify['w1mpro']), usemask=False, asrecarray=False)
      
      params_classify['longitude_l']=gl_classify
      params_classify['latitude_b']=gb_classify
      
      
#classify      
    
    
      classifyparams = np.zeros((params_classify.size, len(par_keys)))
 
      for i, key in enumerate(par_keys):
        classifyparams[:, i] = params_classify[key]
      
#trainingset: params_trainingset
# to classify: params_classify

      probas_qso = rf_qso.predict_proba(classifyparams)
      probas_qso = probas_qso[:, 1]

      #probas[:,1] selects probabilities that each sample is in class 1

# output classification probability p_QSO for each object

     
      f = open(outdirectory+filename_result,'w')
      f.write('objid ra dec longitude_l latitude_b maxloglik q stellar_locus omega_best tau_best errorweightedmeanmag_g errorweightedmeanmag_r errorweightedmeanmag_i errorweightedmeanmag_z errorweightedmeanmag_y strucfuncmeanmag_g strucfuncmeanmag_r strucfuncmeanmag_i strucfuncmeanmag_z strucfuncmeanmag_y gbandepochs_cleaned rbandepochs_cleaned ibandepochs_cleaned zbandepochs_cleaned ybandepochs_cleaned w1mpro w2mpro w1sigmpro w2sigmpro wisera wisedec sigrawise sigdecwise W12 W12Err sdss_mag_u sdss_mag_g sdss_mag_r sdss_mag_i sdss_mag_z  sdss_magerr_u sdss_magerr_g sdss_magerr_r sdss_magerr_i sdss_magerr_z p_RRLyrae EBV\n')  
      
    
      for objid, ra, dec, longitude_l, latitude_b, maxloglik, q, stellar_locus, omega_best, tau_best, errorweightedmeanmag_g, errorweightedmeanmag_r, errorweightedmeanmag_i, \
      errorweightedmeanmag_z, errorweightedmeanmag_y, strucfuncmeanmag_g, strucfuncmeanmag_r, strucfuncmeanmag_i, strucfuncmeanmag_z, \
      strucfuncmeanmag_y, gbandepochs_cleaned, rbandepochs_cleaned, ibandepochs_cleaned, zbandepochs_cleaned, ybandepochs_cleaned, \
      w1mpro, w2mpro, w1sigmpro, w2sigmpro, wisera, wisedec, sigrawise, sigdecwise, W12, W12Err, \
      corr_sdss_mag_u, corr_sdss_mag_g, corr_sdss_mag_r, corr_sdss_mag_i, corr_sdss_mag_z, sdss_magerr_u, sdss_magerr_g, sdss_magerr_r, \
      sdss_magerr_i, sdss_magerr_z, prob, ebv \
      in zip(params_classify["objid"], params_classify["ra"],params_classify["dec"],params_classify["longitude_l"],params_classify["latitude_b"],params_classify["maxloglik"],params_classify["q"], \
      params_classify["stellar_locus"],params_classify["omega_best"],params_classify["tau_best"],params_classify["errorweightedmeanmag_g"], \
      params_classify["errorweightedmeanmag_r"],params_classify["errorweightedmeanmag_i"],params_classify["errorweightedmeanmag_z"], \
      params_classify["errorweightedmeanmag_y"], \
      params_classify["strucfuncmeanmag_g"], \
      params_classify["strucfuncmeanmag_r"],params_classify["strucfuncmeanmag_i"],params_classify["strucfuncmeanmag_z"],params_classify["strucfuncmeanmag_y"], \
      params_classify["gbandepochs_cleaned"],params_classify["rbandepochs_cleaned"],params_classify["ibandepochs_cleaned"], \
      params_classify["zbandepochs_cleaned"],params_classify["ybandepochs_cleaned"],\
      params_classify["w1mpro"],params_classify["w2mpro"],params_classify["w1sigmpro"],params_classify["w2sigmpro"],\
      params_classify["wisera"],params_classify["wisedec"],params_classify["sigrawise"],params_classify["sigdecwise"],\
      params_classify["W12"],params_classify["W12Err"], params_classify["corr_sdss_mag_u"], params_classify["corr_sdss_mag_g"],\
      params_classify["corr_sdss_mag_r"],params_classify["corr_sdss_mag_i"],params_classify["corr_sdss_mag_z"],\
      params_classify["sdss_magerr_u"], \
      params_classify["sdss_magerr_g"],params_classify["sdss_magerr_r"],params_classify["sdss_magerr_i"],params_classify["sdss_magerr_z"] ,\
      probas_qso,params_classify["EBV"]):
        f.write("%s %.8f %.8f %.8f %.8f %.8f %.8f %d %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %d %d %d %d %d %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f \n" % \
        (objid, ra, dec, longitude_l, latitude_b, maxloglik, q, stellar_locus, omega_best, tau_best, errorweightedmeanmag_g, errorweightedmeanmag_r, errorweightedmeanmag_i, errorweightedmeanmag_z, errorweightedmeanmag_y ,\
        strucfuncmeanmag_g, strucfuncmeanmag_r, strucfuncmeanmag_i, strucfuncmeanmag_z, strucfuncmeanmag_y,\
        gbandepochs_cleaned, rbandepochs_cleaned, ibandepochs_cleaned, zbandepochs_cleaned, ybandepochs_cleaned, \
        w1mpro, w2mpro, w1sigmpro, w2sigmpro, wisera, wisedec, sigrawise, sigdecwise, W12, W12Err, \
        corr_sdss_mag_u, corr_sdss_mag_g, corr_sdss_mag_r, corr_sdss_mag_i, corr_sdss_mag_z, \
        sdss_magerr_u, sdss_magerr_g, sdss_magerr_r, sdss_magerr_i, sdss_magerr_z, \
        prob,ebv ))
      f.close()  
      



