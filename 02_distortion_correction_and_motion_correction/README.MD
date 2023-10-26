# Distortion Correction and Motion Correction#


Simple demo to:

-	1. perform distortion correction with **fsl_prepare_fieldmap** and **fugue**
-	2. perform motion correction with **3dvolreg**
-	calculate a concatenated matrix with **cat_matvec** to perform one-step undistortion with **3dAllineate**



<br>
<br>

The dataset comes from s574-session01 and comprises the following images:





#### Anatomical Images: ####
- 	sub-s574_t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_UNI_DEN.nii




####Fieldmapping Images:####

- 	sub-s574_acq-0p8mm_magnitude2.nii.gz
- 	sub-s574_acq-0p8mm_phasediff.nii.gz
- 	sub-s574_acq-1p3mm_magnitude2.nii.gz
- 	sub-s574_acq-1p3mm_phasediff.nii.gz


####Functional Images:####

- 	sub-s574_task-rest_acq-0p8mm_bold.nii.gz
- 	sub-s574_task-rest_acq-1p3mm_bold.nii.gz


<br>

Distortion correction is performed by estimating a distortion map with ***fsl_prepare_fieldmap*** and applying it to the image by using ***fugue***. The detailed processing steps comprise:

- preparing fieldmap_mag
	- bias correcting: **fsl_anat --strongbias --nocrop --noreg --nosubcortseg --noseg -i input -o output**
	- skullstripping: **bet input output -R**
	- eroding: **fslmaths input -ero output**
	- registration to the input image: **epi_reg --epi=input --t1=structural --t1brain=structuralBet --out=output**
	
<br>

- preparing fieldmap_phase:
	- convert to radians: **slmaths input -div 2048 -sub 1 -mul 3.14159 output -odt float**
	- unwrapping the phase image: **prelude -c input -o output**
	- registration to the input image: **epi_reg --epi=input --t1=structural --t1brain=structuralBet --out=output**



<br>

- calculating a fieldmap: **fsl\_prepare\_fieldmap SIEMENS phaseImg magImg output deltaTE**
- unwarping with the fieldmap: **fugue -i inputEPI --dwell=dwell\_time --loadfmap=fieldmap -u output**
- motion correction: **3dvolreg -twopass -1Dmatrix\_save matrix\_output -prefix epi_motionCorrected epi\_input**
- matrix concatenation with **cat\_matvec -ONELINE anat_ns+tlrc::WARP_DATA -I matrix\_distortion -I matrix\_motion > concatenated\_matrix**
- one-step unwarping: **3dAllineate -base anatomical -1Dmatrix\_apply concentrated\_matrix -mast\_dxyz 3 -float -prefix epi\_output epi\_input**


<br>
The dwell time is determined by extracting it from a .json file associated with a .nii file.

