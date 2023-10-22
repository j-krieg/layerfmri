# Layer-fMRI Processing Pipeline #


This is a collection of code for the layer-fMRI processing pipeline.

## Tasks ##
Current Status:

- Minimal Example of FSL-based distortion correction
<br>
<br>

To-Do:

- registration to the standard MNI152 template
- noise estimation with NORDIC
- motion correction
- BOLD response


## Sequence Overview ##

- CEST (chemical exchange saturation transfer)
	- CEST\_MIMOSA (Multiple interleaved mode saturation, See [Whole-brain quantitative CEST MRI at 7T ](https://pubmed.ncbi.nlm.nih.gov/33634505/)	
	- CEST\_WASABI (simultaneous mapping of the water shift and B1)
	- CEST\_T1_SATREC (control saturation recovery) See [Calibration of arterial spin labeling data](https://onlinelibrary.wiley.com/doi/pdfdirect/10.1002/mrm.28000)
<br><br>
- T2s Faruk Gulban, See [A scalable method to improve gray matter segmentation at ultra high field MRI](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0198335) and [Introduction (p. 24)](http://www.81bones.net/mri/mri_introSlides.pdf)



- T2-star-mapping: used to map inhomogenities


- MP2RAGE (Magnetization Prepared - RApid Gradient Echo)
	[MP-RAGE vs. MP2RAGE](https://mriquestions.com/mp-rage-v-mr2rage.html) :"The MP2RAGE sequence, by comparison, uses two Turbo-FLASH GRE readouts between each inversion pulse."
	

- T2 SPACE: [mrimaster](https://mrimaster.com/characterise-image-3d-tse/#:~:text=This%20makes%20the%20SPACE%20sequence,echo%20sequence%20used%20in%20MRI.): "SPACE (Sampling Perfection with Application optimized Contrasts using different flip-angle Evolutions) is a 3D TSE (turbo spin echo) sequence used in magnetic resonance imaging (MRI). It’s known for providing high-resolution isotropic 3D images, which can be reformatted in any plane without loss of image quality. This makes the SPACE sequence particularly useful for imaging structures with complex anatomy or for cases where multiplanar reconstructions are needed."


- field mapping methods ( [source](https://andysbrainbook.readthedocs.io/en/latest/OpenScience/OS/BIDS_Overview.html) )
	- GRE field mapping (gradient echo): field mapping with magnitude and phase difference
	- fmap-SE-AP & fmap-SE-PA (spin echoes with opposite phase encoding directions)

<be>

## Available Datasets ##

    +---s574_session01
    |   +---anatomical
    |   |   10-t2_space_sag_p9_test_grappa.nii
    |   |   11-t2_space_sag_p9_test_ir.nii
    |   |   5-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_INV1.nii
    |   |   6-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_T1_Images.nii
    |   |   8-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_UNI_DEN.nii
    |   |   9-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_INV2.nii
    |   |
    |   +---fieldmap
    |   |   12-gre_field_mapping_HCP_1pp3mm3_e1.nii
    |   |   12-gre_field_mapping_HCP_1pp3mm3_e2.nii
    |   |   13-gre_field_mapping_HCP_1pp3mm3_e2_ph.nii
    |   |   17-gre_field_mapping_HCP_0pp8mm3_e1.nii
    |   |   17-gre_field_mapping_HCP_0pp8mm3_e2.nii
    |   |   18-gre_field_mapping_HCP_0pp8mm3_e2_ph.nii
    |   |
    |   +---functional
    |   |   14-GRE_ep2d_bold_tra_p2_dirAP_phasemap_PEpolar.nii
    |   |   15-GRE_ep2d_bold_tra_p2_dirPA_phasemap_PEpolar.nii
    |   |   16-GRE_ep2d_bold_tra_p2_sms4_tr1700_te30_fa75_peAP_bw1984_res1ppmm3.nii
    |   |   19-GRE_ep2d_bold_tra_p2_sms4_highres_dirAP_phasemap_PEpolar.nii
    |   |   20-GRE_ep2d_bold_tra_p2_sms4_highres_dirPA_phasemap_PEpolar.nii
    |   |   21-GRE_ep2d_bold_tra_p2_sms4_tr1830_te23_fa72_peAP_res0pp8mm3.nii
    |   |
    |   \---others
    |   1-AAhead_scout.nii
    |   2-AAhead_scout_MPR_sag.nii
    |   3-AAhead_scout_MPR_cor.nii
    |   4-AAhead_scout_MPR_tra.nii
    |
    |
    |
    |
    |  
    \---s574_session02
    +---anatomical
    |   10-t2_space_sag_p9_0pp8mm3.nii
    |   11-T2s_0.8mm_tra_3TE_AP_FarukGulban2022_Images_e1.nii
    |   11-T2s_0.8mm_tra_3TE_AP_FarukGulban2022_Images_e2.nii
    |   11-T2s_0.8mm_tra_3TE_AP_FarukGulban2022_Images_e3.nii
    |   11-T2s_0.8mm_tra_3TE_AP_FarukGulban2022_Images_e4.nii
    |   11-T2s_0.8mm_tra_3TE_AP_FarukGulban2022_Images_e5.nii
    |   11-T2s_0.8mm_tra_3TE_AP_FarukGulban2022_Images_e6.nii
    |   5-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_INV1.nii
    |   6-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_T1_Images.nii
    |   8-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_UNI-DEN.nii
    |   9-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_INV2.nii
    |
    +---fieldmap
    |   12-T2s_0.8mm_tra_3TE_AP_FarukGulban2022_e1_ph.nii
    |   12-T2s_0.8mm_tra_3TE_AP_FarukGulban2022_e2_ph.nii
    |   12-T2s_0.8mm_tra_3TE_AP_FarukGulban2022_e3_ph.nii
    |   12-T2s_0.8mm_tra_3TE_AP_FarukGulban2022_e4_ph.nii
    |   12-T2s_0.8mm_tra_3TE_AP_FarukGulban2022_e5_ph.nii
    |   12-T2s_0.8mm_tra_3TE_AP_FarukGulban2022_e6_ph.nii
    |   13-T2s_0.8mm_tra_3TE_AP_FarukGulban2022_T2Star_Maps_.nii
    |   15-WB_B1map_CP_B1.nii
    |   16-WB_B1map_EP_B1.nii
    |
    +---functional
    |   17-CEST_MIMOSA_1p00uT.nii
    |   18-CEST_WASABI.nii
    |   19-CEST_MIMOSA_1p90uT.nii
    |   20-CEST_MIMOSA_3p10uT.nii
    |   21-CEST_MIMOSA_3p75uT.nii
    |   22-CEST_MIMOSA_6p25uT.nii
    |   23-CEST_T1_SATREC.nii
    |
    \---others
    1-AAhead_scout.nii
    2-AAhead_scout_MPR_sag.nii
    3-AAhead_scout_MPR_cor.nii
    4-AAhead_scout_MPR_tra.nii
