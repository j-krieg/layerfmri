# Simple Distortion Correction Demo with FSL #



The dataset comes from s574-session01 and comprises the following images:



Structural Images

- 	5-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_INV1.nii
- 	5-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_INV1.nii
- 	6-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_T1_Images.nii
- 	8-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_UNI_DEN.nii
- 	9-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_INV2.nii

<br>


Fieldmapping Images

- 	17-gre_field_mapping_HCP_0pp8mm3_e2.nii
- 	18-gre_field_mapping_HCP_0pp8mm3_e2_ph.nii


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
- unwarping witht the fieledmap: **fugue -i inputEPI --dwell=dwell\_time --loadfmap=fieldmap -u output**



<br>
The dwell time can be determined by
>     dwell_time = 1 / (bandwidth * phase_encoding_steps)
    
wherein bandwidth and the number of phase encoding steps can be extracted from the DICOM header:

     (0018,0095)  PixelBandwidth: Reciprocal of the total sampling period, in hertz per pixel.
     (0018,0089)  NumberOfPhaseEncodingSteps
     (0018, 9058) MR Acquisition Frequency Encoding Steps              320
     (0018, 9231) MR Acquisition Phase Encoding Steps in-plane         320
     (0018, 9232) MR Acquisition Phase Encoding Steps out-of-plane     240
