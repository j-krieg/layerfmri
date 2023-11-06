# fMRI Preprocessing Pipeline #

The current preprocessing pipeline comprises the steps of:

- undistortion
- bias field correction
- motion correction

<br>


In this particular version, the estimation of the bias field has been improved:

    def biasfieldCorrection (input, output):
	    biasfield = naming (input, "biasFieldEstimation", temp_path)
	    biasfield_smooth = naming(input, "smooth", temp_path)
	    
	    firstimg = naming (input, "firstimg", temp_path)# first image from functional
	    firstbiascorr= naming(input, "firstbiascorr", temp_path)			# first image bias corrected
	    firstbiascorr_restore = firstbiascorr.replace(".nii", "_restore.nii")   	# naming convention comes from fast
	    
	    
	    tasks = [   "fslroi <input> <firstimg> 0 1",				# select first image from time series (0 refers to the index, 1 refers to the number of images to be extracted)
	    		"fast -B -o <firstbiascorr> <firstimg>",			# output is (...)_restore.nii.gz, needs to be renamed in next steps
	   		"rm <firstbiascorr>",
	    		"mv <firstbiascorr_restore> <firstbiascorr>",
	    		"fslmaths <firstimg> -div <firstbiascorr> <biasfield>", 	# Want to calculate bias field to apply it too all other images in time series. Bias field is multiplicative
	    		"fslmaths <biasfield> -s 1 <biasfield_smooth>", 		# needed to not depend on noise in first image of the time series (-s 3 refers to 3 voxels)
	    		"fslmaths <input> -div <biasfield_smooth> <output>" 		# Remove bias field from entire time series.
		    ]
      
	    for task in tasks:
	        task = task.replace("<input>", input)
	        task = task.replace("<output>", output)
	        task = task.replace("<biasfield>", biasfield)
	        task = task.replace("<biasfield_smooth>", biasfield_smooth)
	        task = task.replace("<firstimg>", firstimg)
	        task = task.replace("<firstbiascorr>", firstbiascorr)
	        task = task.replace("<firstbiascorr_restore>", firstbiascorr_restore)
	
	        run(task)
	    



As a result, the preprocessing pipeline performs the following steps:


	jobs = []
	
	fieldmap_mag        = fieldmap_path + "sub-s574_acq-0p8mm_magnitude2.nii.gz"
	fieldmap_phase      = fieldmap_path + "sub-s574_acq-0p8mm_phasediff.nii.gz"
	functional          = functional_path +"sub-s574_task-rest_acq-0p8mm_bold.nii.gz"
	jobs.append([anatomicalT1, anatomicalT2, fieldmap_mag, fieldmap_phase, functional])
	
	fieldmap_mag        = fieldmap_path + "sub-s574_acq-1p3mm_magnitude2.nii.gz"
	fieldmap_phase      = fieldmap_path + "sub-s574_acq-1p3mm_phasediff.nii.gz"
	functional          = functional_path +"sub-s574_task-rest_acq-1p3mm_bold.nii.gz"
	jobs.append([anatomicalT1, anatomicalT2, fieldmap_mag, fieldmap_phase, functional])




	for job in jobs:
	    [anatomicalT1, anatomicalT2, fieldmap_mag, fieldmap_phase, functional] = job
	
	    # 1. undistortion
	    functional_undistorted = naming(functional, "undistorted", bias_path)
	    biasfieldUndistortion (functional, fieldmap_mag, fieldmap_phase, functional_undistorted)
	
	
	    # 2. bias field correction
	    functional_biascorr = naming(functional_undistorted, "biascorr", bias_path)
	    biasfieldCorrection(functional_undistorted, functional_biascorr)
	
	
	    # 3. motion correction
	    functional_motioncorr = naming (functional_biascorr, "motioncorr", output_path)
	    preprocessing_motionCorrection(functional_biascorr, functional_motioncorr, matrix_motion)
    
