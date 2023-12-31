# BET Algorithm Comparison for 7T T1-Weighted Images #

Ultra High-Field Magnetic Resonance Imaging enables layer-dependent functional analysis. Yet, segmentation algorithms "are not maximally effective on 7T volumes, due to increased data complexity, voxel intensity inhomogeneity (intra‐/inter‐sites) and site‐specific artefacts." ([Svanera et al. (2021)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8559470/)) Therefore, the particular problem is that "fMRI studies crucially depend on accurate and precise delineations of the gray matter (GM) ribbon both at the inner (white matter; WM)." ([Gulban et al. (2018)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5991408/))

<br><br>Brain extraction is the first step before applying atlas-based extraction of specific brain areas.
Five different configurations for brain extraction were tested:

- FSL's bet with T1 only
- FSL's bet with T1 and T2
- ROBEX with T1 only
- ROBEX with T1 and T2
- HD-BET with T1 only

<br><br>
ROBEX is provided on https://www.nitrc.org/projects/robex with the claim:

> Many methods have been proposed in the literature, but they often: 1. work well on certain datasets but fail on others; 2. require case-specific parameter tuning. ROBEX aims for robust skull-stripping across datasets with no parameter settings.

<br><br>
HD-BET is provided on https://github.com/MIC-DKFZ/HD-BET with the claim:
> HD-BET outperformed five publicly available brain extraction algorithms (FSL BET, AFNI 3DSkullStrip, Brainsuite BSE, ROBEX and BEaST)



<br><br>
They can be applied to an image by the following prompts:


    # 1. classic bet with only T1
    bet <T1> <output> -m
    
    # 2. FSL bet with T1 and T2
    bet <T1> <output> -A2 <T2_regT1>
    
    # 3. ROBEX only T1 
    ./others/ROBEX/runROBEX.sh -i <T1> -o out.nii.gz
    mv others/ROBEX/out.nii.gz <output>
    
    # 4. ROBEX T1 and T2
    ./others/ROBEX/runROBEX.sh -i <T1> -o out.nii.gz -t <T2_regT1>
    mv others/ROBEX/out.nii.gz <output>
    
    # 5. HD-bet (only T1)
    hd-bet -i <T1> -o <output> -device cpu -mode fast -tta 0


Tests with native and bias-corrected T1-weighted iamges were performed.

Results: See video on FAUbox.
