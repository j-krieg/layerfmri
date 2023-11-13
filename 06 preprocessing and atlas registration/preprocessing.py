import subprocess
import os
import shutil
import nibabel as nib
import numpy as np
import json
import time


import ants

# git clone https://github.com/ANTsX/ANTsPy.git
# cd ANTsPy
# python setup.py install


output_path       = "output/"
temp_path         = "temp/"
reg_path          = "registration/" # folder to save registration results (processing takes some time)
bias_path         = "biascorr/"     # folder to save bias corrected images (first step in the procssing pipeline, routine task)

if not os.path.exists(output_path):
    os.makedirs(output_path)

if not os.path.exists(temp_path):
    os.makedirs(temp_path)
else:
    task = "rm -rf " + temp_path
    #subprocess.Popen(task.split(), stdout=subprocess.PIPE)

    

anatomical_path     = "data/sub-s574/anat/"
fieldmap_path       = "data/sub-s574/fmap/"
functional_path     = "data/sub-s574/func/"



anatomicalT1         = anatomical_path  + "8-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_UNI_DEN.nii"
anatomicalT2         = anatomical_path  + "sub-s574_acq-grappatest_T2w.nii.gz"
jobs = []

if False:
    fieldmap_mag        = fieldmap_path + "sub-s574_acq-0p8mm_magnitude2.nii.gz"
    fieldmap_phase      = fieldmap_path + "sub-s574_acq-0p8mm_phasediff.nii.gz"
    functional          = functional_path +"sub-s574_task-rest_acq-0p8mm_bold.nii.gz"
    jobs.append([anatomicalT1, anatomicalT2, fieldmap_mag, fieldmap_phase, functional])

if True:
    fieldmap_mag        = fieldmap_path + "sub-s574_acq-1p3mm_magnitude2.nii.gz"
    fieldmap_phase      = fieldmap_path + "sub-s574_acq-1p3mm_phasediff.nii.gz"
    functional          = functional_path +"sub-s574_task-rest_acq-1p3mm_bold.nii.gz"
    jobs.append([anatomicalT1, anatomicalT2, fieldmap_mag, fieldmap_phase, functional])





def run(cmd):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

def wait():
    print("Press any key to continue...")
    input()
    print("Script continued!")


def bet_demo (T1, T2, outputNii): # demo of different brain extraction methods. Performs bet on T1. Sme methods require auxiliary T2 image. Demonstrate bet problem for 7T
    betfolder = os.path.dirname(outputNii)



    if ((subprocess.run(['dpkg', '-s', "hd-bet"], capture_output=True)).returncode == 0):
        print(
            "hd-bet not installed. \nrun to install:\n\ngit clone https://github.com/MIC-DKFZ/HD-BET\ncd HD-BET\npip install -e .")


    jobs = []

    if False:
        # perform registration of T2 to T1:
        T2_regT1 = betfolder + "T2_regT1.nii"
        job = "fnirt --in=<T2> --ref=<T1> --iout=<T2_regT1>"
        jobs.append(job)
    else:
        T2_regT1 = "T2_regT1.nii.gz"



    # 1. classic bet with only T1
    output = outputNii.replace(".nii", "_01-FSLbetT1.nii")
    job = "bet <T1> <output> -m"
    job = job.replace("<output>", output)
    jobs.append(job)
    # 14 seconds

    # 2. FSL bet with T1 and T2
    output = outputNii.replace(".nii", "_02-FSLbetT1T2.nii")
    job = "bet <T1> <output> -A2 <T2_regT1>"
    job = job.replace("<output>", output)
    jobs.append(job)

    # 3. ROBEX only T1                                                        https://www.nitrc.org/projects/robex
    # claim:    Many methods have been proposed in the literature, but they often:
    #               - work well on certain datasets but fail on others.
    #               - require case-specific parameter tuning
    #           ROBEX aims for robust skull-stripping across datasets with no parameter settings.
    # download (~600 MB): https://www.nitrc.org/frs/download.php/5994/ROBEXv12.linux64.tar.gz//?i_agree=1&download_now=1
    output = outputNii.replace(".nii", "_03-robexT1.nii")
    job = "./others/ROBEX/runROBEX.sh -i <T1> -o out.nii.gz"
    jobs.append (job)
    job = "mv others/ROBEX/out.nii.gz <output>"
    job = job.replace("<output>", output)
    jobs.append(job)
    # 75 seconds

    # 4. ROBEX T1 and T2
    output = outputNii.replace(".nii", "_04-robexT1T2.nii")
    job = "./others/ROBEX/runROBEX.sh -i <T1> -o out.nii.gz -t <T2_regT1>"
    jobs.append (job)
    job = "mv others/ROBEX/out.nii.gz <output>"
    job = job.replace("<output>", output)
    jobs.append(job)
    # 88 seconds

    # 5. HD-bet (only T1)                                                       https://github.com/MIC-DKFZ/HD-BET
    # claim: "HD-BET outperformed five publicly available brain extraction algorithms (FSL BET, AFNI 3DSkullStrip, Brainsuite BSE, ROBEX and BEaST)"
    output = outputNii.replace(".nii", "_05-hd-bet.nii")
    job = "hd-bet -i <T1> -o <output> -device cpu -mode fast -tta 0"
    job = job.replace("<output>", output)
    jobs.append(job)
    # 398 seconds

    # 5. HD-bet not fast                                                       https://github.com/MIC-DKFZ/HD-BET
    # claim: "HD-BET outperformed five publicly available brain extraction algorithms (FSL BET, AFNI 3DSkullStrip, Brainsuite BSE, ROBEX and BEaST)"
    output = outputNii.replace(".nii", "_06-hd-bet-long.nii")
    job = "hd-bet -i <T1> -o <output> -device cpu"
    job = job.replace("<output>", output)
    jobs.append(job)



    for job in jobs:
        job = job.replace("<T1>",         T1)
        job = job.replace("<T2>",         T2)
        job = job.replace("<T2_regT1>",   T2_regT1)

        print("\n\n" + job + "\n")
        start_time = time.time()
        run(job)
        print(f"Elapsed Time: {round(time.time() - start_time)} seconds")





def get_dwell_time(nii):
    # Extracts the dwell time from a .json file associated with an .nii file


    # Alternative: extract from DICOM header
    # https://nipy.org/nibabel/dicom/dicom_niftiheader.html
    # https://dicom.innolitics.com/ciods/mr-image/mr-image/00180095
    # (0018,0095) PixelBandwidth: Reciprocal of the total sampling period, in hertz per pixel.
    # (0018,0089) NumberOfPhaseEncodingSteps
    # (0018, 9058) MR Acquisition Frequency Encoding Steps              320
    # (0018, 9231) MR Acquisition Phase Encoding Steps in-plane         320
    # (0018, 9232) MR Acquisition Phase Encoding Steps out-of-plane     240

    # bandwidth                   = 240                       # Bandwidth per pixel phase encoding steps
    # phase_encoding_steps        = 320                       # Number of phase encoding steps
    # dwell_time = 1 / (bandwidth * phase_encoding_steps)  # seconds per pixel


    jsonfile = (nii.replace(".nii", ".json")).replace(".gz", "")

    try:
        with open(jsonfile) as file:
            data = json.load(file)
            dwell_time = float(data['DwellTime'])
    except (FileNotFoundError, KeyError):
        print("Error: could not extract 'DwellTime' from " + jsonfile)
        dwell_time = 0.0000026
        print("Guessed dwell time. ################################")
        #quit()

    print ("dwell time: " + str(dwell_time))
    return dwell_time


# brain extraction
def run_bet(input_img, output_img):
    #subprocess.check_output(['bet', input_img, output_img])

    task = "bet <input_img> <output_img>"
    task = task.replace("<input_img>",      input_img)
    task = task.replace("<output_img>",     output_img)


    process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    




def run_reg(input_epi, input_structural, output_registered):
    # registration of an EPI to a structural image by using FSL's epi_reg

    # need to perform bet on structural image first:
    input_structural_bet = temp_path + os.path.basename(input_structural.replace(".nii", "_bet.nii"))
    run_bet(input_structural, input_structural_bet)


    task = "epi_reg --epi=<input_epi> --t1=<structural> --t1brain=<structural_bet> --out=<output>"
    task = task.replace ("<input_epi>",         input_epi)
    task = task.replace ("<structural>",        input_structural)
    task = task.replace ("<structural_bet>",    input_structural_bet)
    task = task.replace ("<output>",            output_registered)


    process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    #print ("Done: Registration with FSL epi_reg: " + output_registered)







def preprocessing_estimatefieldmap (fieldmap_mag, fieldmap_phase):
    tasks = ["bet <fieldmap_mag> <fieldmap_mag_bet>",
            "fslmaths <fieldmap_mag_bet> -ero <fieldmap_mag_bet_ero>",
            "fsl_prepare_fieldmap SIEMENS <fieldmap_phase> <fieldmap_mag_bet_ero> <fieldmap_out> 1.02"]

    fieldmap_out = naming(fieldmap_mag, "biasfieldEstimation", output_path)

    fieldmap_mag_bet        = temp_path + os.path.basename(fieldmap_mag.replace(".nii", "_bet.nii"))
    fieldmap_mag_bet_ero    = temp_path + os.path.basename(fieldmap_mag_bet.replace(".nii", "_ero.nii"))


    for task in tasks:
        task = task.replace("<fieldmap_mag>",           fieldmap_mag)
        task = task.replace("<fieldmap_mag_bet>",       fieldmap_mag_bet)
        task = task.replace("<fieldmap_mag_bet_ero>",   fieldmap_mag_bet_ero)
        task = task.replace("<fieldmap_phase>",         fieldmap_phase)
        task = task.replace("<fieldmap_out>",           fieldmap_out)
        process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()

    print("fieldmap computed: " + fieldmap_out)
    return fieldmap_out





def preprocessing_applyfieldmap(input_img, fieldmap_distortion):
    # uses fugue(FMRIB's Utility for Geometrically Unwarping EPIs) to perform an unwarping of an EPI image based on fieldmap data.

    output_corrected = naming(input_img, "undistorted", output_path)

    dwell_time = get_dwell_time(input_img) # taken from json file (alternatively, DICOM header)

    print("- fieldmap application for undistortion...")
    task = "fugue -i <input_epi> --dwell=<dwell_time> --unwarpdir=y- --loadfmap=<fieldmap> -u <result>" # "-savematrix <matrix_distortionCorr>" (--unwarpdir: "Use x, y, z, x-, y- or z- only.")
    #task += " --smooth3=1" # Gaussian smoothing with sigma in mm (optional)
    task = task.replace("<input_epi>",                  input_img)
    task = task.replace("<dwell_time>",                 str(dwell_time))
    task = task.replace("<fieldmap>",                   fieldmap_distortion)
    task = task.replace("<result>",                     output_corrected)
    task += " --verbose"


    process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    print("Done: Application of fieldmap to epi for distortion correction.")
    return output_corrected


def preprocessing_motioncorrection (input):
    output = naming (input, "motioncorr", output_path)
    
    task = "3dvolreg -twopass -prefix <epi_motionCorrected> <epi_input>"

    task = task.replace("<epi_input>", input)
    task = task.replace("<epi_motionCorrected>", output)
    run (task)
    return output



def preprocessing_biasfield4D (input):
    # input is a functional 4D file
    # estimates bias field based on first entry and removes same bias field from all other elements in the time series
    # assumtion: For small distances, the bias field does not change much, in particular, if the functional images are already motion corrected.
    # Please Note: The bias field is not estimated for every element of the time series individually because this is about fMRI data.

    output = naming (input, "biasFieldRemoval", bias_path)
    biasfield = naming (input, "biasFieldEstimation", temp_path)
    biasfield_smooth = naming(input, "smooth", temp_path)

    firstimg = naming (input, "firstimg", temp_path)                        # first image from functional
    firstbiascorr= naming(input, "firstbiascorr", temp_path)                # first image bias corrected
    firstbiascorr_restore = firstbiascorr.replace(".nii", "_restore.nii")   # naming convention comes form fast


    tasks = [   "fslroi <input> <firstimg> 0 1",                            # select first image from time series (0 refers to the index, 1 refers to the number of images to be extracted)
                "fast -B -o <firstbiascorr> <firstimg>",                    # output is (...)_restore.nii.gz, needs to be renamed in next steps
                "rm <firstbiascorr>",
                "mv <firstbiascorr_restore> <firstbiascorr>",
                "fslmaths <firstimg> -div <firstbiascorr> <biasfield>",     # Want to calculate bias field to apply it too all other images in time series. Bias field is multiplicative
                "fslmaths <biasfield> -s 1 <biasfield_smooth>",             # needed to not depend on noise in first image of the time series (-s 3 refers to 3 voxels)
                "fslmaths <input> -div <biasfield_smooth> <output>"]        # Remove bias field from entire time series.



    for task in tasks:
        task = task.replace("<input>", input)
        task = task.replace("<output>", output)
        task = task.replace("<biasfield>", biasfield)
        task = task.replace("<biasfield_smooth>", biasfield_smooth)
        task = task.replace("<firstimg>", firstimg)
        task = task.replace("<firstbiascorr>", firstbiascorr)
        task = task.replace("<firstbiascorr_restore>", firstbiascorr_restore)

        #print (task)
        run(task)

    return output





def naming (source_file, suffix, target_directory):
    # input:    source_file: directory + filename to original file
    #           suffix: new suffix e.g., "biascorr" to name the new file "_biascorr.nii"
    #           target_directory: replace directory in path


    if not target_directory.endswith("/"):
        target_directory += "/"

    result = target_directory + os.path.basename(source_file.replace(".nii", ("_" + suffix + ".nii")))
    if not ".gz" in result:
        result += ".gz"

    return result



def preprocessing_undistortionAnts(input, fieldmap_mag, fieldmap_phase):

    output = naming(input, "undistortedResampledAnts", bias_path)


    # Load the magnitude and phase maps
    magnitude = ants.image_read(fieldmap_mag)
    phase = ants.image_read(fieldmap_phase)

    # Perform phase unwrapping
    unwrapped_phase = ants.unwrap(phase)

    # Generate a field map from the unwrapped phase map
    field_map = ants.fieldmap_from_unwrapped_phase(magnitude, unwrapped_phase)

    # Load the BOLD time series
    bold = ants.image_read(input)

    # Apply distortion correction and resample to 0.3 x 0.3 x 0.3 mm
    distortion_corrected = ants.resample_image(bold, field_map, (0.3, 0.3, 0.3), use_voxel_space=True)

    # Save the distortion-corrected and resampled image
    distortion_corrected.to_filename(output)


    return output





def preprocessing_undistortion(input, fieldmap_mag, fieldmap_phase):
    #   input: functional image
    #   returns filename of output
    #   output: undistorted functional image

    output = naming(input, "undistorted", bias_path)

    # 1. fieldmap estimation
    fieldmap_distortion = preprocessing_estimatefieldmap(fieldmap_mag, fieldmap_phase)  # uses fsl_prepare_fieldmap (and bet for skullstripping, fslmaths -ero for eroding, fsl_anat --strongbias for bias correction, epi-reg for registration, fslmaths -div/sub/mul for radian conversion)
   
    # 2. apply the fieldmap to undistort the input image
    output = preprocessing_applyfieldmap(input, fieldmap_distortion) # uses fugue

    return output



def coregistrationsOriginal (functional, anatomical, MNI, atlas):
    # typical use: [T1w_functionalspace, atlas_functionalspace] = coregistrations (functional_current, anatomicalT1, constant_MNItemplate, constant_atlas)
    #example:   input:
    #               - MNI template
    #               - Glasser atlas
    #               - T1w
    #               - intermediary functional (undistorted, motion corrected, bias field removed)
    #           output:
    #               - T1w in functional space
    #               - Glasser atlas in functional space

    T1w_functionalspace         = output_path   + "anatomical_functionalspace.nii.gz"             # naming(anatomical, "functionalspace", output_path)
    MNI_anatomicalspace         = temp_path     + "MNI_anatomicalspace.nii.gz"                    # temp_path + "MNI_anatomicalspace.nii.gz" # gets too long otherwise
    MNI_functionalspace         = output_path   + "MNI_functionalspace.nii.gz"                    # naming(MNI, "functionalspace", output_path)
    atlas_functionalspace       = output_path   + "atlas_functionalspace.nii.gz"                  # naming(atlas, "functionalspace", output_path)
    anatomical_bet              = temp_path     + "anatomical_bet.nii.gz"                         # gets too long otherwise


    matrix_T1_to_bold = "temp/matrix1.mat"
    matrix_MNI_to_T1 = "temp/matrix2.mat"
    matrix_MNI_to_bold = "temp/matrix3.mat"




    jobs = []


    # Step 1: Register T1.nii.gz to bold.nii.gz with rigid transformation and get matrix

    job = "bash antsRegistrationSyN.sh -d 3 -f <functional> -m <anatomical> -t 'r' -o temp/ANTsOutput"
    jobs.append(job)
    # parameters:
    #               f: fixed image
    #               t: 'r' (rigid), 'a': rigid + affine (2 stages), 's': rigid + affine + deformable syn (3 stages)
    #               d: dimension (2 or 3, not 4)
    #               m: moving image (that needs transform)
    #               o: outputfile
    # output:
    # [outputname]Affine.mat
    # [outputname]Warped.nii.gz

    job = "rm <T1w_functionalspace>"
    jobs.append(job)
    job = "rm <matrix_T1_to_bold>"
    jobs.append(job)
    job = "rm temp/ANTsOutputInverseWarped.nii.gz"
    jobs.append(job)
    job = "cp temp/ANTsOutputWarped.nii.gz <T1w_functionalspace>"
    jobs.append(job)
    job = "cp temp/ANTsOutput0GenericAffine.mat <matrix_T1_to_bold>"
    jobs.append(job)



    # Step 2: Register MNI.nii.gz to T1.nii.gz with nonlinear transformation and get matrix
    # 2a: anatomical image needs bet, otherwise, MNI template is aligned with skull
    
    #job = "bash hd-bet -i <anatomical> -o <anatomical_bet> -device cpu -mode fast -tta 0" # "The options -mode fast and -tta 0 will disable test time data augmentation (speedup of 8x) and use only one model instead of an ensemble of five models for the prediction."
    #jobs.append(job)

    job = "./others/ROBEX/runROBEX.sh -i <anatomical> -o out.nii.gz"
    jobs.append(job)

    job = "rm <anatomical_bet>"
    jobs.append(job)
    job = "cp others/ROBEX/out.nii.gz <anatomical_bet>"
    jobs.append(job)

    # 2b: registrations
    job = "bash antsRegistrationSyN.sh -d 3 -f <anatomical_bet> -m <MNI> -t 's' -o temp/ANTsOutput"
    jobs.append(job)
    # parameters:
    #               f: fixed image
    #               t: 'r' (rigid), 'a': rigid + affine (2 stages), 's': rigid + affine + deformable syn (3 stages)
    #               d: dimension (2 or 3, not 4)
    #               m: moving image (that needs transform)
    #               o: outputfile
    # output:
    # [outputname]Affine.mat
    # [outputname]Warped.nii.gz
    job = "rm <MNI_anatomicalspace>"
    jobs.append(job)
    job = "rm <matrix_MNI_to_T1>"
    jobs.append(job)
    job = "rm temp/ANTsOutputInverseWarped.nii.gz"
    jobs.append(job)
    job = "mv temp/ANTsOutputWarped.nii.gz <MNI_anatomicalspace>"
    jobs.append(job)
    job = "mv temp/ANTsOutput0GenericAffine.mat <matrix_MNI_to_T1>"
    jobs.append(job)


    # Step 3: Concatenate matrix1 and matrix2 (MNI -> anatomical -> functional) and apply to atlas
    job = "antsApplyTransforms -d 3 -i <atlas> -o temp/ANTsOutputFINAL -r <functional> -t <matrix_MNI_to_T1> -t <matrix_T1_to_bold>"
    jobs.append(job)#
    # parameters:
    #               -t transformation matrix (can use multiple, each preceded by a -t)
    #               -i input
    #               -o output
    #               -r reference image: "For warping input images, the reference image defines the spacing, origin, size, and direction of the output warped image."
    #




    for job in jobs:
        job = job.replace("<functional>",               functional)
        job = job.replace("<anatomical>",               anatomical)
        job = job.replace("<anatomical_bet>",           anatomical_bet)
        job = job.replace("<MNI>",                      MNI)
        job = job.replace("<atlas>",                    atlas)
        job = job.replace("<matrix_MNI_to_T1>",         matrix_MNI_to_T1)
        job = job.replace("<matrix_T1_to_bold>",        matrix_T1_to_bold)
        job = job.replace("<matrix_MNI_to_bold>",       matrix_MNI_to_bold)
        job = job.replace("<MNI_anatomicalspace>",      MNI_anatomicalspace)
        job = job.replace("<MNI_functionalspace>",      MNI_functionalspace)
        job = job.replace("<T1w_functionalspace>",      T1w_functionalspace)
        job = job.replace("<atlas_functionalspace>",    atlas_functionalspace)

        print("\n" + job + "\n")
        wait()
        run(job)


    print("===========================")
    quit()



    result = [T1w_functionalspace, atlas_functionalspace]
    return result



def resampling (input, resolution): # resolution: e.g., 0.3 for 0.3 mm
    output = "temp/upsampled.nii.gz" # naming(input, "upsample", temp_path)
    task = "3dresample -input <input> -prefix <output> -dxyz <res> <res> <res>"

    task = task.replace("<input>", input)
    task = task.replace("<output>", output)
    task = task.replace("<res>", str(resolution))

    print(task)
    run(task)

    wait()
    return output



def coregistrations(functional, anatomical, MNI, atlas):
    # typical use: [T1w_functionalspace, atlas_functionalspace] = coregistrations (functional_current, anatomicalT1, constant_MNItemplate, constant_atlas)
    # example:   input:
    #               - MNI template
    #               - Glasser atlas
    #               - T1w
    #               - intermediary functional (undistorted, motion corrected, bias field removed)
    #           output:
    #               - T1w in functional space
    #               - Glasser atlas in functional space

    T1w_functionalspace = output_path + "anatomical_functionalspace.nii.gz"  # naming(anatomical, "functionalspace", output_path)
    anatomical_functionalspace_bet = output_path + "anatomical_functionalspace_bet.nii.gz"
    MNI_anatomicalspace = temp_path + "MNI_anatomicalspace.nii.gz"  # temp_path + "MNI_anatomicalspace.nii.gz" # gets too long otherwise
    MNI_functionalspace = output_path + "MNI_functionalspace.nii.gz"  # naming(MNI, "functionalspace", output_path)
    atlas_functionalspace = output_path + "atlas_functionalspace.nii.gz"  # naming(atlas, "functionalspace", output_path)
    anatomical_bet = temp_path + "anatomical_bet.nii.gz"  # gets too long otherwise
    functional_bet = temp_path + "functional_bet.nii.gz"
    functional_first = temp_path + "functional_first.nii.gz"
    functional_first_bet = temp_path + "functional_first_bet.nii.gz"
    anatomical_bet_functional = output_path + "anatomical_bet_functionalspace.nii.gz"


    matrix_T1_to_bold = "temp/matrix1.mat"
    matrix_MNI_to_T1 = "temp/matrix2.mat"
    matrix_MNI_to_functional = "temp/matrix3.mat"

    jobs = []

    # Step 1: Register T1.nii.gz to bold.nii.gz with rigid transformation

    job = "fslroi <functional> <functional_first> 0 1"
    jobs.append(job)

    job = "bet <functional_first> <functional_first_bet> -f 0.3 -R" # 0.3 is needed to remove less tissue. If too much brain tissue is remved (default: f=0.5), registration produces bad results
    jobs.append(job)

    job = "hd-bet -i <anatomical> -o <anatomical_bet> -device cpu -mode fast -tta 0"
    jobs.append(job)

    job = "flirt -in <anatomical_bet> -ref <functional_first_bet> -noresample -out <anatomical_bet_functional>"
    jobs.append(job)







    # Step 3: Register MNI.nii.gz to anatomical_functionalspace_bet with nonlinear transformation and get matrix
    job = "bash antsRegistrationSyN.sh -d 3 -f <anatomical_bet_functional> -m <MNI> -t 's' -o temp/ANTsOutput"
    jobs.append(job)
    # parameters:
    #               f: fixed image
    #               t: 'r' (rigid), 'a': rigid + affine (2 stages), 's': rigid + affine + deformable syn (3 stages)
    #               d: dimension (2 or 3, not 4)
    #               m: moving image (that needs transform)
    #               o: outputfile
    # output:
    # [outputname]Affine.mat
    # [outputname]Warped.nii.gz
    job = "rm <MNI_functionalspace>"
    jobs.append(job)
    job = "rm <matrix_MNI_to_functional>"
    jobs.append(job)
    job = "mv temp/ANTsOutputWarped.nii.gz <MNI_functionalspace>"
    jobs.append(job)
    job = "mv temp/ANTsOutput0GenericAffine.mat <matrix_MNI_to_functional>"
    jobs.append(job)




    # Step 3: Apply matrix to atlas
    #job = "antsApplyTransforms -d 3 -i <atlas> -o temp/ANTsOutputFINAL -r <anatomical_functionalspace_bet> -t <matrix_MNI_to_functional>"
    job = "antsApplyTransforms -d 3 -i <atlas> -r <anatomical_bet_functional> -o <atlas_functionalspace> -n NearestNeighbor -t <matrix_MNI_to_functional>"
    jobs.append(job)
    # parameters:
    #               -t transformation matrix (can use multiple, each preceded by a -t)
    #               -i input
    #               -o output
    #               -r reference image: "For warping input images, the reference image defines the spacing, origin, size, and direction of the output warped image."
    #

    for job in jobs:
        job = job.replace("<functional>", functional)
        job = job.replace("<anatomical_bet_functional>", anatomical_bet_functional)
        job = job.replace("<functional_bet>", functional_bet)
        job = job.replace("<functional_first>", functional_first)
        job = job.replace("<functional_first_bet>", functional_first_bet)
        job = job.replace("<anatomical>", anatomical)
        job = job.replace("<anatomical_bet>", anatomical_bet)
        job = job.replace("<anatomical_functionalspace_bet>", anatomical_functionalspace_bet)
        job = job.replace("<MNI>", MNI)
        job = job.replace("<atlas>", atlas)
        job = job.replace("<matrix_MNI_to_T1>", matrix_MNI_to_T1)
        job = job.replace("<matrix_T1_to_bold>", matrix_T1_to_bold)
        job = job.replace("<matrix_MNI_to_functional>", matrix_MNI_to_functional)
        job = job.replace("<MNI_anatomicalspace>", MNI_anatomicalspace)
        job = job.replace("<MNI_functionalspace>", MNI_functionalspace)
        job = job.replace("<T1w_functionalspace>", T1w_functionalspace)
        job = job.replace("<atlas_functionalspace>", atlas_functionalspace)
        job = job.replace("<anatomical_functionalspace_bet>", anatomical_functionalspace_bet)



        print("\n" + job + "\n")
        #wait()
        run(job)

    print("===========================")
    quit()

    result = [T1w_functionalspace, atlas_functionalspace]
    return result


def functional_extractROI (functional_preprocessed, ROI, atlas_functionalspace):
    # typical use: functional_ROI = functional_extractROI (functional_preprocessed, ROI, atlas_functionalspace)
    # input:
    #           - preprocessed functional 4D image timeseries (undistortied, bias field removed, motion corrected)
    #           - ROI: keyword
    #           - atlas in functional space

    glasserparcellation = {                         # source: https://bitbucket.org/dpat/tools/src/master/REF/ATLASES/Glasser_2016_Table.xlsx
        'V1': [1, 'Primary Visual Cortex'],
        'A1': [24, 'Primary Auditory Cortex'],
        '1': [51, 'Area 1'],                        # postcentral gyrus
        '2': [52, 'Area 2'],                        # postcentral gyrus
        '3a': [53, 'Area 3a']                       # postcentral gyrus
    }
    


def deobliqueUpsample (input, targetres):
    output = naming(input, "deobliqueUpsample", temp_path)

    job = "3dWarp -deoblique -NN -newgrid <targetres> -prefix <output> <input>"
    job = job.replace("<targetres>",    str(targetres)) # e.g., 0.3
    job = job.replace("<input>",        input)
    job = job.replace("<output>",       output)

    run(job)
    return output


for job in jobs:
    [anatomicalT1, anatomicalT2, fieldmap_mag, fieldmap_phase, functional] = job

    functional_current = functional

    # 1. undistortion (i.e., correct warping from bias field)
    functional_current =  preprocessing_undistortion (functional_current, fieldmap_mag, fieldmap_phase)

    # 2. deoblique and upsampling
    # targetres = 0.3
    #functional_current = deobliqueUpsample(functional_current, targetres)

    # 3. bias field correction (Bias field correction is performed before motion correction because the bias field is independent of motion.)
    functional_current = preprocessing_biasfield4D(functional_current)


    # 4. motion correction to distortion corrrected functional images (with additional output of transformation matrix):
    functional_current = preprocessing_motioncorrection(functional_current)

    #functional_current = "before/bold_undistorted_motioncorr_biasFieldEstimation.nii.gz"


    #functional_current = resampling (functional_current, 0.3)

    # 4. coregistrations
    # input:
    #   - MNI template
    #   - Glasser atlas
    #   - T1w
    #   - intermediary functional (undistorted, motion corrected, bias field removed)
    # output:
    #   - T1w in functional space
    #   - Glasser atlas in functional space

    constant_MNItemplate = "atlas/MNI152_2009_template.nii.gz"
    constant_atlas = "atlas/HCPMMP1_on_MNI152_ICBM2009a_nlin_hd.nii.gz" # Glasser
    [T1w_functionalspace, atlas_functionalspace] = coregistrations (functional_current, anatomicalT1, constant_MNItemplate, constant_atlas)
    quit()


    # 5. brain extraction: uses structural image (already in functional space) for brain extraction of functional images. Result is preprocessed fMRI (undistortion, motion correction, bias field removed, brain extracted)
#    functional_preprocessed = bet_functional (functional_current, T1w_functionalspace)



    # 6. atlas-based ROI extraction
    ROI = "motor"
#    functional_ROI = functional_extractROI (functional_preprocessed, ROI, atlas_functionalspace)



