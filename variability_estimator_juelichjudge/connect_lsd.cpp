/* connect_lsd.cpp as part of variability_estimator_lsd0.9
this code is used to prepare for estimate_multiband_variabiliy.cpp
*/

#line 2 "connect_lsd.cpp"

assert(Sobj_id[0] == sizeof(*obj_id));    // Make sure we've got a contiguous array


// stream through the input arrays

int number_of_objects = Nobj_id[0];  

//std::cout<< "number of lines: "<<number_of_objects<<std::endl;



//these are the temporary vectors, from input - they have to be preprocessed

  gsl_vector *psf_inst_mag_vec_tmp=0;
  gsl_vector *psf_inst_mag_sig_vec_tmp=0;
  gsl_vector *ra_vec_tmp=0;
  gsl_vector *dec_vec_tmp=0;
  gsl_vector *filterid_vec_tmp=0;
  gsl_vector *smf_fn_vec_tmp=0;
  gsl_vector *exptime_vec_tmp=0;
  gsl_vector *mjd_obs_vec_tmp=0;
  gsl_vector *ap_mag_vec_tmp=0;
  gsl_vector *psf_qf_perfect_vec_tmp=0;
  gsl_vector *ucalmag_vec_tmp=0;
  gsl_vector *zp_vec_tmp=0;	    
  std::vector<int> flags_vec_tmp;	  
    

//these are the vectors that are processed

  gsl_vector *gband;
  gsl_vector *rband;
  gsl_vector *iband;
  gsl_vector *zband;
  gsl_vector *yband;
  gsl_vector *gbanderr;
  gsl_vector *rbanderr;
  gsl_vector *ibanderr; 
  gsl_vector *zbanderr; 
  gsl_vector *ybanderr;
  gsl_vector *gbandtime_vec; 
  gsl_vector *rbandtime_vec; 
  gsl_vector *ibandtime_vec;  
  gsl_vector *zbandtime_vec; 
  gsl_vector *ybandtime_vec;  
  std::vector<int> gbandflag_vec;  
  std::vector<int> rbandflag_vec;
  std::vector<int> ibandflag_vec; 
  std::vector<int> zbandflag_vec; 
  std::vector<int> ybandflag_vec;
  gsl_vector *gband_psf_qf_perfect_vec;
  gsl_vector *rband_psf_qf_perfect_vec;
  gsl_vector *iband_psf_qf_perfect_vec;
  gsl_vector *zband_psf_qf_perfect_vec;
  gsl_vector *yband_psf_qf_perfect_vec;
  gsl_vector *gband_ap_mag;
  gsl_vector *rband_ap_mag;
  gsl_vector *iband_ap_mag;
  gsl_vector *zband_ap_mag;
  gsl_vector *yband_ap_mag;
  gsl_vector *gband_psf_inst_mag;
  gsl_vector *rband_psf_inst_mag; 
  gsl_vector *iband_psf_inst_mag; 
  gsl_vector *zband_psf_inst_mag; 
  gsl_vector *yband_psf_inst_mag; 
  double object_ra;
  double object_dec;

  std::string obj_id_string;
  int dim; //this is the number of epochs for the object

  
  //this is the output

  
  int processed=0;
  double ra_output=0;
  double dec_output=0;
  double maxloglik_output=0;
  double q_output=0;
  int stellar_locus_output=0;
  double omega_best_output=0;
  double tau_best_output=0;
  double errorweightedmeanmag_g_output=0;
  double errorweightedmeanmag_r_output=0;
  double errorweightedmeanmag_i_output=0;
  double errorweightedmeanmag_z_output=0;
  double errorweightedmeanmag_y_output=0;
  double strucfuncmeanmag_g_output=0;
  double strucfuncmeanmag_r_output=0;
  double strucfuncmeanmag_i_output=0;
  double strucfuncmeanmag_z_output=0;
  double strucfuncmeanmag_y_output=0;
	
  int gbandepochs_cleaned_output=0;
  int rbandepochs_cleaned_output=0;
  int ibandepochs_cleaned_output=0;
  int zbandepochs_cleaned_output=0;
  int ybandepochs_cleaned_output=0;
  
  int idx1=0;

  int gband_dim=0;
  int rband_dim=0;	
  int iband_dim=0;	
  int zband_dim=0;	
  int yband_dim=0;	
  int idx=0;

  int g_idx=0;
  int r_idx=0;
  int i_idx=0;
  int z_idx=0;
  int y_idx=0;
	
  int idx2=0;
  int i=0;
  int j=0;

	
	
  for(i = 0, j = 0; i != number_of_objects; i = j) //size: all lines, means all rows from input file, can contain many objects
  {
	j = i;
	
	dim=0;
	
	
	gband=0;
	rband=0;
	iband=0;
	zband=0;
	yband=0;
	gbanderr=0;
	rbanderr=0;
	ibanderr=0;
	zbanderr=0;
	ybanderr=0;
	gbandtime_vec=0;
	rbandtime_vec=0;
	ibandtime_vec=0;
	zbandtime_vec=0;
	ybandtime_vec=0;
	  
	gband_psf_qf_perfect_vec=0;
	rband_psf_qf_perfect_vec=0;
	iband_psf_qf_perfect_vec=0;
	zband_psf_qf_perfect_vec=0;
	yband_psf_qf_perfect_vec=0;
	gband_ap_mag=0;
	rband_ap_mag=0;
	iband_ap_mag=0;
	zband_ap_mag=0;
	yband_ap_mag=0;
	gband_psf_inst_mag=0;
	rband_psf_inst_mag=0;
	iband_psf_inst_mag=0;
	zband_psf_inst_mag=0;
	yband_psf_inst_mag=0;
	ra=0;
	dec=0;
		
		
	psf_inst_mag_vec_tmp = gsl_vector_calloc(300);
	psf_inst_mag_sig_vec_tmp = gsl_vector_calloc(300);
	ra_vec_tmp = gsl_vector_calloc(300);
	dec_vec_tmp = gsl_vector_calloc(300);
	filterid_vec_tmp = gsl_vector_calloc(300);
	smf_fn_vec_tmp = gsl_vector_calloc(300);
	mjd_obs_vec_tmp = gsl_vector_calloc(300);
	ap_mag_vec_tmp = gsl_vector_calloc(300);
	psf_qf_perfect_vec_tmp = gsl_vector_calloc(300);
	ucalmag_vec_tmp = gsl_vector_calloc(300);
	zp_vec_tmp = gsl_vector_calloc(300);
	flags_vec_tmp.clear();
	
		  	  
	std::stringstream ss; 
	ss << obj_id[i];
	obj_id_string=ss.str();

//	std::cout<<"OBJ_ID "<<obj_id[i]<<std::endl;
	
//	std::cout<<"OBJID "<<obj_id_string<<std::endl;
	
	
	idx1=0;
        while(j != number_of_objects && obj_id[i] == obj_id[j])
	{
	    
	  //get the data
	  
	  gsl_vector_set(psf_inst_mag_vec_tmp,idx1, PSF_INST_MAG1(j));
	  gsl_vector_set(psf_inst_mag_sig_vec_tmp,idx1, PSF_INST_MAG_SIG1(j));
	  gsl_vector_set(ra_vec_tmp,idx1, RA1(j));
	  gsl_vector_set(dec_vec_tmp,idx1, DEC1(j));
	  gsl_vector_set(filterid_vec_tmp,idx1, BAND1(j));
	  gsl_vector_set(mjd_obs_vec_tmp,idx1, MJD_OBS1(j));
	  gsl_vector_set(ap_mag_vec_tmp,idx1, AP_MAG1(j));
	  gsl_vector_set(psf_qf_perfect_vec_tmp,idx1, PSF_QF_PERFECT1(j));
	  gsl_vector_set(ucalmag_vec_tmp,idx1, UCALMAG1(j));
	  gsl_vector_set(zp_vec_tmp,idx1, ZP1(j));	     
	  flags_vec_tmp.push_back(FLAGS1(j));
	   
	
	  idx1++;
	  j++;
	}
	
	//preprocess them, output to vector for light curve
	gband_dim=0;
	rband_dim=0;	
	iband_dim=0;	
	zband_dim=0;	
	yband_dim=0;	
	  

        for(idx=0;idx<mjd_obs_vec_tmp->size;idx++)
	{   
	    
	//    std::cout<<gsl_vector_get(mjd_obs_vec_tmp,idx)<<" "<<gsl_vector_get(ucalmag_vec_tmp,idx)<<" "<<gsl_vector_get(psf_inst_mag_sig_vec_tmp,idx)<<std::endl;
	    
	      if((gsl_vector_get(mjd_obs_vec_tmp,idx)!=0) && (gsl_vector_get(ucalmag_vec_tmp,idx)>1) 
		&& (gsl_vector_get(psf_inst_mag_sig_vec_tmp,idx)>0) 
	        && (gsl_vector_get(zp_vec_tmp,idx)>10) 
	        && (gsl_vector_get(zp_vec_tmp,idx)<40) 
		&& (!isnan(gsl_vector_get(zp_vec_tmp,idx)))
	      )         
		//removing entries with error code, -99.99 mag, and bad zp
	      {
		
		//calculate dim
		
			if(gsl_vector_get(filterid_vec_tmp,idx)==0)
			{
			  gband_dim++;
			}
			if(gsl_vector_get(filterid_vec_tmp,idx)==1)
			{
			  rband_dim++;
			}
			if(gsl_vector_get(filterid_vec_tmp,idx)==2)
			{
			  iband_dim++;
			}
			if(gsl_vector_get(filterid_vec_tmp,idx)==3)
			{
			  zband_dim++;
			}
			if(gsl_vector_get(filterid_vec_tmp,idx)==4)
			{
			  yband_dim++;
			}
		}    
	}

    if(gband_dim+rband_dim+iband_dim+zband_dim+yband_dim>0)
    {    
	if(gband_dim>0)
	{
	  gband = gsl_vector_calloc(gband_dim);  
	  gbanderr = gsl_vector_calloc(gband_dim); 
	  gbandtime_vec = gsl_vector_calloc(gband_dim);
	  gband_psf_qf_perfect_vec = gsl_vector_calloc(gband_dim);
	  gband_ap_mag = gsl_vector_calloc(gband_dim);
	  gband_psf_inst_mag = gsl_vector_calloc(gband_dim);
	  gbandflag_vec.resize(gband_dim);
	}
    
	if(rband_dim>0)
	{
	  rband = gsl_vector_calloc(rband_dim);  
	  rbanderr = gsl_vector_calloc(rband_dim); 
	  rbandtime_vec = gsl_vector_calloc(rband_dim);
	  rband_psf_qf_perfect_vec = gsl_vector_calloc(rband_dim);
	  rband_ap_mag = gsl_vector_calloc(rband_dim);
	  rband_psf_inst_mag = gsl_vector_calloc(rband_dim);
	  rbandflag_vec.resize(rband_dim);
	}
	
	if(iband_dim>0)
	{
	  iband = gsl_vector_calloc(iband_dim);  
	  ibanderr = gsl_vector_calloc(iband_dim); 
	  ibandtime_vec = gsl_vector_calloc(iband_dim);
	  iband_psf_qf_perfect_vec = gsl_vector_calloc(iband_dim);
	  iband_ap_mag = gsl_vector_calloc(iband_dim);
	  iband_psf_inst_mag = gsl_vector_calloc(iband_dim);
	  ibandflag_vec.resize(iband_dim);
	}
	
	if(zband_dim>0)
	{
	  zband = gsl_vector_calloc(zband_dim);  
	  zbanderr = gsl_vector_calloc(zband_dim); 
	  zbandtime_vec = gsl_vector_calloc(zband_dim);
	  zband_psf_qf_perfect_vec = gsl_vector_calloc(zband_dim);
	  zband_ap_mag = gsl_vector_calloc(zband_dim);
	  zband_psf_inst_mag = gsl_vector_calloc(zband_dim);
	  zbandflag_vec.resize(zband_dim);
	}
	
	if(yband_dim>0)
	{
	  yband = gsl_vector_calloc(yband_dim);  
	  ybanderr = gsl_vector_calloc(yband_dim); 
	  ybandtime_vec = gsl_vector_calloc(yband_dim);
	  yband_psf_qf_perfect_vec = gsl_vector_calloc(yband_dim);
	  yband_ap_mag = gsl_vector_calloc(yband_dim);
	  yband_psf_inst_mag = gsl_vector_calloc(yband_dim);
	  ybandflag_vec.resize(yband_dim);
	}
	
	g_idx=0;
	r_idx=0;
	i_idx=0;
	z_idx=0;
	y_idx=0;
	
	
	  
	for(idx2=0;idx2<mjd_obs_vec_tmp->size;idx2++)
	{
    //removing entries with error code, -99.99 mag
	//removing entries with   zp < 10 || zp > 40
	  
	  if(gsl_vector_get(filterid_vec_tmp,idx2)==0)
	  {
	  
	      if((gsl_vector_get(mjd_obs_vec_tmp,idx2)!=0) && (gsl_vector_get(ucalmag_vec_tmp,idx2)>1) 
		&& (gsl_vector_get(psf_inst_mag_sig_vec_tmp,idx2)>0) 
	        && (gsl_vector_get(zp_vec_tmp,idx2)>10) 
	        && (gsl_vector_get(zp_vec_tmp,idx2)<40) 
		&& (!isnan(gsl_vector_get(zp_vec_tmp,idx2)))
	      )         
	      {
		gsl_vector_set(gband,g_idx,gsl_vector_get(ucalmag_vec_tmp,idx2));
		gsl_vector_set(gbanderr,g_idx,2.5* log10(1+gsl_vector_get(psf_inst_mag_sig_vec_tmp,idx2)));
		gsl_vector_set(gbandtime_vec,g_idx,gsl_vector_get(mjd_obs_vec_tmp,idx2));
		gsl_vector_set(gband_psf_qf_perfect_vec,g_idx,gsl_vector_get(psf_qf_perfect_vec_tmp,idx2));
		gsl_vector_set(gband_ap_mag,g_idx,gsl_vector_get(ap_mag_vec_tmp,idx2));
		gsl_vector_set(gband_psf_inst_mag,g_idx,gsl_vector_get(psf_inst_mag_vec_tmp,idx2));
		gbandflag_vec[g_idx]= flags_vec_tmp[idx2];
		g_idx++;
	      }
	  
	  } 
      
	  if(gsl_vector_get(filterid_vec_tmp,idx2)==1)
	  {
	 
	      if((gsl_vector_get(mjd_obs_vec_tmp,idx2)!=0) && (gsl_vector_get(ucalmag_vec_tmp,idx2)>1) 
		&& (gsl_vector_get(psf_inst_mag_sig_vec_tmp,idx2)>0) 
	        && (gsl_vector_get(zp_vec_tmp,idx2)>10) 
	        && (gsl_vector_get(zp_vec_tmp,idx2)<40) 
		&& (!isnan(gsl_vector_get(zp_vec_tmp,idx2)))
	      )         
	      {
		gsl_vector_set(rband,r_idx,gsl_vector_get(ucalmag_vec_tmp,idx2));
		gsl_vector_set(rbanderr,r_idx,2.5* log10(1+gsl_vector_get(psf_inst_mag_sig_vec_tmp,idx2)));
		gsl_vector_set(rbandtime_vec,r_idx,gsl_vector_get(mjd_obs_vec_tmp,idx2));
		gsl_vector_set(rband_psf_qf_perfect_vec,r_idx,gsl_vector_get(psf_qf_perfect_vec_tmp,idx2));
		gsl_vector_set(rband_ap_mag,r_idx,gsl_vector_get(ap_mag_vec_tmp,idx2));
		gsl_vector_set(rband_psf_inst_mag,r_idx,gsl_vector_get(psf_inst_mag_vec_tmp,idx2));
		rbandflag_vec[r_idx]= flags_vec_tmp[idx2];
		r_idx++;
	      }
	  
	  } 
	 
  
	  if(gsl_vector_get(filterid_vec_tmp,idx2)==2)
	  {
	      if((gsl_vector_get(mjd_obs_vec_tmp,idx2)!=0) && (gsl_vector_get(ucalmag_vec_tmp,idx2)>1) 
		&& (gsl_vector_get(psf_inst_mag_sig_vec_tmp,idx2)>0) 
	        && (gsl_vector_get(zp_vec_tmp,idx2)>10) 
	        && (gsl_vector_get(zp_vec_tmp,idx2)<40) 	
		&& (!isnan(gsl_vector_get(zp_vec_tmp,idx2)))
	      )         
	      {
		gsl_vector_set(iband,i_idx,gsl_vector_get(ucalmag_vec_tmp,idx2));
		gsl_vector_set(ibanderr,i_idx,2.5* log10(1+gsl_vector_get(psf_inst_mag_sig_vec_tmp,idx2)));
		gsl_vector_set(ibandtime_vec,i_idx,gsl_vector_get(mjd_obs_vec_tmp,idx2));
		gsl_vector_set(iband_psf_qf_perfect_vec,i_idx,gsl_vector_get(psf_qf_perfect_vec_tmp,idx2));
		gsl_vector_set(iband_ap_mag,i_idx,gsl_vector_get(ap_mag_vec_tmp,idx2));
		gsl_vector_set(iband_psf_inst_mag,i_idx,gsl_vector_get(psf_inst_mag_vec_tmp,idx2));
		ibandflag_vec[i_idx]= flags_vec_tmp[idx2];
		i_idx++;
	      }
	  
	  } 

	  if(gsl_vector_get(filterid_vec_tmp,idx2)==3)
	  {
	  
	      if((gsl_vector_get(mjd_obs_vec_tmp,idx2)!=0) && (gsl_vector_get(ucalmag_vec_tmp,idx2)>1) 
		&& (gsl_vector_get(psf_inst_mag_sig_vec_tmp,idx2)>0) 
	        && (gsl_vector_get(zp_vec_tmp,idx2)>10) 
	        && (gsl_vector_get(zp_vec_tmp,idx2)<40) 
		&& (!isnan(gsl_vector_get(zp_vec_tmp,idx2)))
	      )         
	      {
		gsl_vector_set(zband,z_idx,gsl_vector_get(ucalmag_vec_tmp,idx2));
		gsl_vector_set(zbanderr,z_idx,2.5* log10(1+gsl_vector_get(psf_inst_mag_sig_vec_tmp,idx2)));
		gsl_vector_set(zbandtime_vec,z_idx,gsl_vector_get(mjd_obs_vec_tmp,idx2));
		gsl_vector_set(zband_psf_qf_perfect_vec,z_idx,gsl_vector_get(psf_qf_perfect_vec_tmp,idx2));
		gsl_vector_set(zband_ap_mag,z_idx,gsl_vector_get(ap_mag_vec_tmp,idx2));
		gsl_vector_set(zband_psf_inst_mag,z_idx,gsl_vector_get(psf_inst_mag_vec_tmp,idx2));
		zbandflag_vec[z_idx]= flags_vec_tmp[idx2];
		z_idx++;
	      }
	  
	  } 
	
	  if(gsl_vector_get(filterid_vec_tmp,idx2)==4)
	  {
	  
	      if((gsl_vector_get(mjd_obs_vec_tmp,idx2)!=0) && (gsl_vector_get(ucalmag_vec_tmp,idx2)>1) 
		&& (gsl_vector_get(psf_inst_mag_sig_vec_tmp,idx2)>0) 
		&& (gsl_vector_get(zp_vec_tmp,idx2)>10) 
	        && (gsl_vector_get(zp_vec_tmp,idx2)<40) 	
		&& (!isnan(gsl_vector_get(zp_vec_tmp,idx2)))
	      )         
	      {
		gsl_vector_set(yband,y_idx,gsl_vector_get(ucalmag_vec_tmp,idx2));
		gsl_vector_set(ybanderr,y_idx,2.5* log10(1+gsl_vector_get(psf_inst_mag_sig_vec_tmp,idx2)));
		gsl_vector_set(ybandtime_vec,y_idx,gsl_vector_get(mjd_obs_vec_tmp,idx2));
		gsl_vector_set(yband_psf_qf_perfect_vec,y_idx,gsl_vector_get(psf_qf_perfect_vec_tmp,idx2));
		gsl_vector_set(yband_ap_mag,y_idx,gsl_vector_get(ap_mag_vec_tmp,idx2));
		gsl_vector_set(yband_psf_inst_mag,y_idx,gsl_vector_get(psf_inst_mag_vec_tmp,idx2));
		ybandflag_vec[y_idx]= flags_vec_tmp[idx2];
		y_idx++;
	      }
	  
	  } 
	 
	}
    
    
	object_ra = gsl_vector_get(ra_vec_tmp,0);
	object_dec = gsl_vector_get(dec_vec_tmp,0);
		  
	
		
	if(psf_inst_mag_vec_tmp!=0)
	{
	  gsl_vector_free(psf_inst_mag_vec_tmp);
	}
	if(psf_inst_mag_sig_vec_tmp!=0)
	{
	  gsl_vector_free(psf_inst_mag_sig_vec_tmp);
	}
	if(ra_vec_tmp!=0)
	{
	  gsl_vector_free(ra_vec_tmp);
	}
	if(dec_vec_tmp!=0)
	{
	  gsl_vector_free(dec_vec_tmp);
	}
	if(filterid_vec_tmp!=0)
	{
	  gsl_vector_free(filterid_vec_tmp);
	}
	if(mjd_obs_vec_tmp!=0)
	{
	  gsl_vector_free(mjd_obs_vec_tmp);
	}
	if(ap_mag_vec_tmp!=0)
	{
	  gsl_vector_free(ap_mag_vec_tmp);
	}
	if(psf_qf_perfect_vec_tmp!=0)
	{
	   gsl_vector_free(psf_qf_perfect_vec_tmp);
	}
	if(ucalmag_vec_tmp!=0)
	{
	  gsl_vector_free(ucalmag_vec_tmp);
	}
	if(zp_vec_tmp!=0)
	{
	  gsl_vector_free(zp_vec_tmp);
	}	
	if(flags_vec_tmp.size()>0)
	{
	  flags_vec_tmp.clear();
	}

	
	//calculate DRW structure function parameters and other values used for classification		

	calculate_classification_params(obj_id_string, gband, rband, iband, zband, yband, 
				     gbanderr, rbanderr, ibanderr, zbanderr, ybanderr,
				     gbandtime_vec, rbandtime_vec, ibandtime_vec, zbandtime_vec, ybandtime_vec, gbandflag_vec,
				     rbandflag_vec, ibandflag_vec, zbandflag_vec, ybandflag_vec, 
				     gband_psf_qf_perfect_vec, rband_psf_qf_perfect_vec, iband_psf_qf_perfect_vec, zband_psf_qf_perfect_vec,
				     yband_psf_qf_perfect_vec,
				     gband_ap_mag, rband_ap_mag, iband_ap_mag, zband_ap_mag, 
				     yband_ap_mag,
				     gband_psf_inst_mag, rband_psf_inst_mag, iband_psf_inst_mag, 
				     zband_psf_inst_mag, yband_psf_inst_mag, 
				     object_ra, object_dec,
				     processed,ra_output,dec_output,maxloglik_output,
				     q_output,stellar_locus_output,omega_best_output,tau_best_output,
				     errorweightedmeanmag_g_output, errorweightedmeanmag_r_output,
				     errorweightedmeanmag_i_output,errorweightedmeanmag_z_output,errorweightedmeanmag_y_output,
				     strucfuncmeanmag_g_output, strucfuncmeanmag_r_output,strucfuncmeanmag_i_output,
				     strucfuncmeanmag_z_output, strucfuncmeanmag_y_output,
				     gbandepochs_cleaned_output,
				     rbandepochs_cleaned_output,ibandepochs_cleaned_output,zbandepochs_cleaned_output,
				     ybandepochs_cleaned_output	);
	
	// fill the output
	
	PROCESSED_FLAG1(0) = processed;
	
	
	OBJ_ID_OUT1(0) = obj_id[0];
	RA_OUT1(0)=ra_output;
	DEC_OUT1(0)=dec_output;
	MAXLOGLIK_OUT1(0)=maxloglik_output;
	Q_OUT1(0)=q_output;
	STELLAR_LOCUS_OUT1(0)=stellar_locus_output;
	OMEGA_BEST_OUT1(0)=omega_best_output;
	TAU_BEST_OUT1(0)=tau_best_output;
	ERRORWEIGHTEDMEANMAG_G_OUT1(0)=errorweightedmeanmag_g_output;
	ERRORWEIGHTEDMEANMAG_R_OUT1(0)=errorweightedmeanmag_r_output;
	ERRORWEIGHTEDMEANMAG_I_OUT1(0)=errorweightedmeanmag_i_output;
	ERRORWEIGHTEDMEANMAG_Z_OUT1(0)=errorweightedmeanmag_z_output;
	ERRORWEIGHTEDMEANMAG_Y_OUT1(0)=errorweightedmeanmag_y_output;
	STRUCFUNCMEANMAG_G_OUT1(0)=strucfuncmeanmag_g_output;
	STRUCFUNCMEANMAG_R_OUT1(0)=strucfuncmeanmag_r_output;
	STRUCFUNCMEANMAG_I_OUT1(0)=strucfuncmeanmag_i_output;
	STRUCFUNCMEANMAG_Z_OUT1(0)=strucfuncmeanmag_z_output;
	STRUCFUNCMEANMAG_Y_OUT1(0)=strucfuncmeanmag_y_output;
	GBANDEPOCHS_CLEANED_OUT1(0)=gbandepochs_cleaned_output;
	RBANDEPOCHS_CLEANED_OUT1(0)=rbandepochs_cleaned_output;
	IBANDEPOCHS_CLEANED_OUT1(0)=ibandepochs_cleaned_output;
	ZBANDEPOCHS_CLEANED_OUT1(0)=zbandepochs_cleaned_output;
	YBANDEPOCHS_CLEANED_OUT1(0)=ybandepochs_cleaned_output;
	
  }
	
	//clean up now
	
	
	if(gband!=0)
	{
	  gsl_vector_free(gband);
	}
	if(gbanderr!=0)
	{
	  gsl_vector_free(gbanderr);
	}
	if(gbandtime_vec!=0)
	{
	  gsl_vector_free(gbandtime_vec);
	}
	if(gband_psf_qf_perfect_vec!=0)
	{
	  gsl_vector_free(gband_psf_qf_perfect_vec);
	}
	if(gband_ap_mag!=0)
	{
	  gsl_vector_free(gband_ap_mag);
	}
	if(gband_psf_inst_mag!=0)
	{
	  gsl_vector_free(gband_psf_inst_mag);
	}

	if(rband!=0)
	{
	  gsl_vector_free(rband);
	}
	
	  
	if(rbanderr!=0)
	{
	  gsl_vector_free(rbanderr);
	}
	if(rbandtime_vec!=0)
	{
	  gsl_vector_free(rbandtime_vec);
	}
	if(rband_psf_qf_perfect_vec!=0)
	{
	  gsl_vector_free(rband_psf_qf_perfect_vec);
	}
	if(rband_ap_mag!=0)
	{
	  gsl_vector_free(rband_ap_mag);
	}
	if(rband_psf_inst_mag!=0)
	{
	  gsl_vector_free(rband_psf_inst_mag);
	}
	if(iband!=0)
	{
	  gsl_vector_free(iband);
	}
	if(ibanderr!=0)
	{
	  gsl_vector_free(ibanderr);
	}
	if(ibandtime_vec!=0)
	{
	  gsl_vector_free(ibandtime_vec);
	}
	  
	if(iband_psf_qf_perfect_vec!=0)
	{
	  gsl_vector_free(iband_psf_qf_perfect_vec);
	}
	if(iband_ap_mag!=0)
	{
	  gsl_vector_free(iband_ap_mag);
	}
	if(iband_psf_inst_mag!=0)
	{
	  gsl_vector_free(iband_psf_inst_mag);
	}
	if(zband!=0)
	{
	  gsl_vector_free(zband);
	}
	if(zbanderr!=0)
	{
	  gsl_vector_free(zbanderr);
	}
	if(zbandtime_vec!=0)
	{
	  gsl_vector_free(zbandtime_vec);
	}
	  
	  
	if(zband_ap_mag!=0)
	{
	  gsl_vector_free(zband_ap_mag);
	}
	if(zband_psf_inst_mag!=0)
	{
	  gsl_vector_free(zband_psf_inst_mag);
	}

	if(yband!=0)
	{
	  gsl_vector_free(yband);
	}  
	 
	  
	if(ybanderr!=0)
	{
	  gsl_vector_free(ybanderr);
	}
	if(ybandtime_vec!=0)
	{
	  gsl_vector_free(ybandtime_vec);
	}
	if(yband_psf_qf_perfect_vec!=0)
	{
	  gsl_vector_free(yband_psf_qf_perfect_vec);
	}
	if(yband_ap_mag!=0)
	{
	  gsl_vector_free(yband_ap_mag);
	}
	if(yband_psf_inst_mag!=0)
	{
	  gsl_vector_free(yband_psf_inst_mag);
	}
	if(gbandflag_vec.size()>0)
	{
	  gbandflag_vec.clear();
	}
	if(rbandflag_vec.size()>0)
	{
	  rbandflag_vec.clear();
	}
	if(ibandflag_vec.size()>0)
	{
	  ibandflag_vec.clear();
	}
	if(zbandflag_vec.size()>0)
	{
	  zbandflag_vec.clear();
	}
	if(ybandflag_vec.size()>0)
	{
	  ybandflag_vec.clear();
	}
	
	
	

  gband=0;
  rband=0;
  iband=0;
  zband=0;
  yband=0;
  gbanderr=0;
  rbanderr=0;
  ibanderr=0; 
  zbanderr=0; 
  ybanderr=0;
  gbandtime_vec=0; 
  rbandtime_vec=0; 
  ibandtime_vec=0;  
  zbandtime_vec=0; 
  ybandtime_vec=0;  
  gband_psf_qf_perfect_vec=0;
  rband_psf_qf_perfect_vec=0;
  iband_psf_qf_perfect_vec=0;
  zband_psf_qf_perfect_vec=0;
  yband_psf_qf_perfect_vec=0;
  gband_ap_mag=0;
  rband_ap_mag=0;
  iband_ap_mag=0;
  zband_ap_mag=0;
  yband_ap_mag=0;
  gband_psf_inst_mag=0;
  rband_psf_inst_mag=0; 
  iband_psf_inst_mag=0; 
  zband_psf_inst_mag=0; 
  yband_psf_inst_mag=0; 
  
  
	  
	  
  }

  
  
 