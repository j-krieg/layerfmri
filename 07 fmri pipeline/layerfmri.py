import os
import shutil
import numpy as np
import json
import time
import psutil # for memory usage


global_simulated = False # simulate runs to re-use file names without running the whole pipeline


# git clone https://github.com/ANTsX/ANTsPy.git
# cd ANTsPy
# python setup.py install


output_path         = "output/"
temp_path           = "temp/"
bias_path           = "biascorr/"     # folder to save bias corrected images (first step in the procssing pipeline, routine task)
atlas_path          = "../atlas/"
mask_path           = output_path + "masks/"
layniiresult_path   = output_path + "layer_results/"

constant_MNItemplate = atlas_path + "MNI152_2009_template.nii.gz"
constant_atlas = atlas_path + "HCPMMP1_on_MNI152_ICBM2009a_nlin_hd.nii.gz"  # Glasser



if not os.path.exists(output_path):
    os.makedirs(output_path)

if not os.path.exists(bias_path):
    os.makedirs(bias_path)

if not os.path.exists(temp_path):
    os.makedirs(temp_path)
else:
    task = "rm -rf -f " + temp_path
    #run(task)

if not os.path.exists(mask_path):
    os.makedirs(mask_path)

if not os.path.exists(layniiresult_path):
    os.makedirs(layniiresult_path)

    
jobs = []

import individualFiles
jobs.append(individualFiles.indivjob())






def sec2MMSS (duration): # e.g., 72 seconds => 01:12
    minutes = duration // 60
    seconds = duration % 60

    mm = str(minutes)
    if minutes < 10:
        mm = "0" + mm

    ss = str(seconds)
    if seconds < 10:
        ss = "0" + ss

    time_format = mm + ":" + ss  # f"{minutes:02d}:{seconds:02d}"
    return time_format


def run(cmd):
    if not global_simulated: # simulate runs to re-use file names without running the whole pipeline


        start_time = time.time()

        largelist = ["recon-all", "flirt", "3dvolreg", "antsRegistrationSyN.sh", "fugue"]  # list of commands that take a long while => show timer
        large = False



        for title in largelist:
            if title in cmd:
                if False: # not working yet (just a vanity function to display a timer and monitor RAM usage during time-demanding operations)
                    RAMusage = runTimer(cmd, title)
                    cmd += "    # max memory usage: " + str(RAMusage) + " GB"
                    large = True
                else:
                    print("running time-demanding operation: " + cmd)


        if not large: # otherwise already performed
            #result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            result = os.popen(cmd).read()

        end_time = time.time()
        duration = int(end_time - start_time)
        time_format = sec2MMSS (duration)
        print(str(time_format) + "  " + cmd)

        with open("run.txt", "a+") as file:
            file.write(time_format + " " + cmd + "\n")


    else:
        print("simulated: " + cmd)


def runTimer(cmd, title): # same as run, but for larger jobs to show a timer.
    # returns maximal RAM usage in GB

    start_time = time.time()
    #process = subprocess.Popen(cmd, shell=True)
    process = os.popen(cmd)


    maxRAM = 0.0

    print("\n")
    while process.poll() is None:  # Check if the process is still running
        end_time = time.time()
        duration = int(end_time - start_time)
        time_format = sec2MMSS (duration)

        memory_usage = psutil.virtual_memory().used / (1024 ** 3)  # Convert bytes to GB
        RAMtext = f" - {memory_usage:.2f} GB"
        if memory_usage > maxRAM:
            maxRAM = memory_usage

        print(time_format + " running: " + title + RAMtext, end="\r")
        time.sleep(1)

    return round(maxRAM, 2)





def wait():
    print("Press any key to continue...")
    input()
    print("Script continued!")





def get_dwell_time(nii):
    # Extracts the effective echo spacing time from a .json file associated with an .nii file
    # comment 2023-12-28: needs to be EffectiveEchoSpacing and not "DwellTime" (otherwise no difference visible, i.e. value too small). direction "y" (and not "y-") is the correct one
    # EffectiveEchoSpacing": 0.000319989
    # "DwellTime": 2.6e-06, 0.0000026

    jsonfile = (nii.replace(".nii", ".json")).replace(".gz", "")

    try:
        with open(jsonfile) as file:
            data = json.load(file)
            #dwell_time = float(data['DwellTime'])             # "dwell time, Dwell Time Equivalent to echo spacing or time between echoes in k-space. " https://neuroimaging-core-docs.readthedocs.io/en/latest/pages/glossary.html
            dwell_time = float(data['EffectiveEchoSpacing']) # "EffectiveEchoSpacing": 0.000319989,
    except (FileNotFoundError, KeyError):
        print("Error: could not extract 'EffectiveEchoSpacing' from " + jsonfile)
        dwell_time = 0.000319989 # used previously
        print("Guessed dwell time. ################################")
        #quit()



    run ("echo dwell time from JSON file: " + str(dwell_time))
    #print ("dwell time: " + str(dwell_time))
    return dwell_time


# brain extraction
def run_bet(input_img, output_img):
    #subprocess.check_output(['bet', input_img, output_img])

    task = "bet <input_img> <output_img>"
    task = task.replace("<input_img>",      input_img)
    task = task.replace("<output_img>",     output_img)


    run(task)
    




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


    run(task)

    #print ("Done: Registration with FSL epi_reg: " + output_registered)







def preprocessing_estimatefieldmap (fieldmap_mag, fieldmap_phase):
    tasks = [       "mri_synthstrip -i <fieldmap_mag> --no-csf -o <fieldmap_mag_bet>",
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
        run(task)

    #print("fieldmap computed: " + fieldmap_out)
    return fieldmap_out





def preprocessing_applyfieldmap(input_img, fieldmap_distortion):
    # uses fugue(FMRIB's Utility for Geometrically Unwarping EPIs) to perform an unwarping of an EPI image based on fieldmap data.

    output_corrected = naming(input_img, "undistorted", output_path)

    dwell_time = get_dwell_time(input_img) # taken from json file (alternatively, DICOM header)

    #print("- fieldmap application for undistortion...")
    task = "fugue -i <input_epi> --dwell=<dwell_time> --unwarpdir=y --loadfmap=<fieldmap> -u <result>" # "-savematrix <matrix_distortionCorr>" (--unwarpdir: "Use x, y, z, x-, y- or z- only.")
    #task += " --smooth3=1" # Gaussian smoothing with sigma in mm (optional)
    task = task.replace("<input_epi>",                  input_img)
    task = task.replace("<dwell_time>",                 str(dwell_time))
    task = task.replace("<fieldmap>",                   fieldmap_distortion)
    task = task.replace("<result>",                     output_corrected)
    task += " --verbose"


    run(task)

    #print("Done: Application of fieldmap to epi for distortion correction.")
    return output_corrected


def preprocessing_motioncorrection (input):
    checkFiles([input])
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

    checkFiles([input])


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

    checkFiles([input, fieldmap_mag, fieldmap_phase])

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




def resampling (input, resolution): # resolution: e.g., 0.3 for 0.3 mm
    output = "<temp/>upsampled.nii.gz" # naming(input, "upsample", temp_path)
    task = "3dresample -input <input> -prefix <output> -dxyz <res> <res> <res>"

    task = task.replace("<input>", input)
    task = task.replace("<output>", output)
    task = task.replace("<temp/>", output_path)
    task = task.replace("<res>", str(resolution))

    run(task)

    wait()
    return output



def coregistrations(functional, functional_reference, functional_reference_up, anatomical, MNI, atlas):
    # typical use: [anatomical_functionalspace(up), atlas_functionalspace(up), T1w_functionalspace_bet, matrix_atlas_to_functional, matrix_anatomical_to_functional] = coregistrations(functional, functional_reference, functional_reference_up, anatomical, MNI, atlas)
    #               - functional_reference
    #               - functional_reference_up (previously upsampled)
    #               - anatomical T1
    #               - MNI template
    #               - Glasser atlas
    #
    #           output:
    #               - amatomical in functional space (upsampled)
    #               - Glasser atlas in functional space (upsampled)
    #               - anatomical_bet_functional (upsamled)
    #               - matrix_atlas_to_functional
    #               - matrix_anatomical_to_functional

    checkFiles([functional, anatomical, MNI, atlas, functional_reference, functional_reference_up])

    T1w_functionalspace =                           output_path + "anatomical_functionalspace.nii.gz"           # naming(anatomical, "functionalspace", output_path)
    anatomical_functionalspace_bet =                output_path + "anatomical_functionalspace_bet.nii.gz"
    MNI_anatomicalspace =                           temp_path + "MNI_anatomicalspace.nii.gz"                    # temp_path + "MNI_anatomicalspace.nii.gz" # gets too long otherwise
    MNI_functionalspace =                           naming(MNI, "functionalspace", output_path)                 # output_path + "MNI_functionalspace.nii.gz"
    atlas_functionalspace =                         naming(atlas, "functionalspace", output_path)               # output_path + "atlas_functionalspace.nii.gz"
    anatomical_bet =                                temp_path + "anatomical_bet.nii.gz"                         # gets too long otherwise
    functional_bet =                                temp_path + "functional_bet.nii.gz"
    functional_first =                              temp_path + "functional_first.nii.gz"
    functional_first_bet =                          temp_path + "functional_first_bet.nii.gz"
    anatomical_bet_functional =                     output_path + "anatomical_bet_functionalspace.nii.gz"
    matrix_anatomical_to_functional =               output_path + "matrix_anatomical_to_functional.mat"
    matrix_MNI_to_functional =                      output_path + "matrix_MNI_to_functional.mat"
    functional_reference_up_bet =                   naming(functional_reference_up, "bet", output_path)
    anatomical_bet_functional_low =                 naming(anatomical_bet_functional, "low", temp_path)
    functional_reference_bet =                      naming(functional_reference, "bet", output_path)

    jobs = []

    # Step 1: Register T1.nii.gz to bold.nii.gz with rigid transformation

    jobs.append("rm <functional_first>") # fslroi cannot overwrite
    jobs.append("fslroi <functional> <functional_first> 0 1")

    if False:
        jobs.append("bet <functional_first> <functional_first_bet> -f 0.3 -R") # 0.3 is needed to remove less tissue. If too much brain tissue is remved (default: f=0.5), registration produces bad results
    else:
        jobs.append("mri_synthstrip -i <functional_reference_up> --no-csf -o <functional_reference_up_bet>")
        jobs.append("mri_synthstrip -i <functional_reference> --no-csf -o <functional_reference_bet>")

    jobs.append("mri_synthstrip -i <anatomical> --no-csf -o <anatomical_bet>")

    # align anatomical image to functional image and keep the original resolution of the anatomical image. Takes 5-25 minutes.
    jobs.append("flirt -in <anatomical_bet> -ref <functional_reference_bet> -noresample -omat <matrix_anatomical_to_functional> -out <anatomical_bet_functional_low>") # will be used in step 3. Faster if the non-upsampled image is used.

    if False:
        jobs.append("flirt -in <anatomical_bet> -ref <functional_reference_up_bet> -omat <matrix_anatomical_to_functional> -out <anatomical_bet_functional>") # includes upsampling the anatomical image to the resolution of functional_reference_up
    else: # re-use matrix from lower spatial resolution to save time.
        #jobs.append("antsApplyTransforms -d 3 -i <anatomical_bet> -r <functional_reference_up> -o <anatomical_bet_functional> -n NearestNeighbor -t <matrix_anatomical_to_functional>") # not working
        jobs.append("flirt -in <anatomical_bet> -ref <functional_reference_up> -applyxfm -init <matrix_anatomical_to_functional> -out <anatomical_bet_functional>") # only takes a minute




    # Step 3: Register MNI.nii.gz to anatomical_functionalspace_bet with nonlinear transformation and get matrix
    jobs.append("bash antsRegistrationSyN.sh -d 3 -f <anatomical_bet_functional_low> -m <MNI> -t 's' -o <temp/>ANTsOutput") # mostly interested in the transformation matrix, therefore, not using the uplampled image.
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
    job = "mv <temp/>ANTsOutputWarped.nii.gz <MNI_functionalspace>"
    jobs.append(job)
    job = "mv <temp/>ANTsOutput0GenericAffine.mat <matrix_MNI_to_functional>"
    jobs.append(job)




    # Step 4: Apply matrix to atlas
    #job = "antsApplyTransforms -d 3 -i <atlas> -o <temp/>ANTsOutputFINAL -r <anatomical_functionalspace_bet> -t <matrix_MNI_to_functional>"
    jobs.append("antsApplyTransforms -d 3 -i <atlas> -r <anatomical_bet_functional> -o <atlas_functionalspace> -n NearestNeighbor -t <matrix_MNI_to_functional>")
    # parameters:
    #               -t transformation matrix (can use multiple, each preceded by a -t)
    #               -i input
    #               -o output
    #               -r reference image: "For warping input images, the reference image defines the spacing, origin, size, and direction of the output warped image."
    #

    for job in jobs:
        job = job.replace("<functional>", functional)
        job = job.replace("<anatomical_bet_functional>", anatomical_bet_functional)
        job = job.replace("<anatomical_bet_functional_low>", anatomical_bet_functional_low)
        job = job.replace("<functional_bet>", functional_bet)
        job = job.replace("<functional_first>", functional_first)
        job = job.replace("<functional_first_bet>", functional_first_bet)
        job = job.replace("<anatomical>", anatomical)
        job = job.replace("<anatomical_bet>", anatomical_bet)
        job = job.replace("<anatomical_functionalspace_bet>", anatomical_functionalspace_bet)
        job = job.replace("<MNI>", MNI)
        job = job.replace("<atlas>", atlas)
        job = job.replace("<temp/>", temp_path)
        job = job.replace("<matrix_anatomical_to_functional>", matrix_anatomical_to_functional)
        job = job.replace("<matrix_MNI_to_functional>", matrix_MNI_to_functional)
        job = job.replace("<MNI_anatomicalspace>", MNI_anatomicalspace)
        job = job.replace("<MNI_functionalspace>", MNI_functionalspace)
        job = job.replace("<T1w_functionalspace>", T1w_functionalspace)
        job = job.replace("<atlas_functionalspace>", atlas_functionalspace)
        job = job.replace("<anatomical_functionalspace_bet>", anatomical_functionalspace_bet)
        job = job.replace("<functional_reference>", functional_reference)
        job = job.replace("<functional_reference_up>", functional_reference_up)
        job = job.replace("<functional_reference_up_bet>", functional_reference_up_bet)
        job = job.replace("<functional_reference_bet>", functional_reference_bet)

        run(job)

    result = [T1w_functionalspace, atlas_functionalspace, anatomical_bet_functional, matrix_MNI_to_functional, matrix_anatomical_to_functional]
    return result


def functional_extractROI (functional, ROI, atlas_functionalspace):
    # typical use: functional_ROI = functional_extractROI (functional_preprocessed, ROI, atlas_functionalspace)
    # input:
    #           - preprocessed functional 4D image timeseries (undistortied, bias field removed, motion corrected)
    #           - ROI: keyword
    #           - atlas in functional space
    # output:
    #           - functional_ROI

    glasserparcellation = {                         # source: https://bitbucket.org/dpat/tools/src/master/REF/ATLASES/Glasser_2016_Table.xlsx
        'V1': [1, 'visualcortex', 'Primary Visual Cortex'],
        'A1': [24, 'auditorycortex', 'Primary Auditory Cortex'],
        '1': [51, 'postcentral', 'Area 1'],                        # postcentral gyrus
        '2': [52, 'postcentral', 'Area 2'],                        # postcentral gyrus
        '3a': [53, 'postcentral', 'Area 3a']                       # postcentral gyrus
    }

    if ROI in glasserparcellation:
        keyValue = glasserparcellation[ROI][0]
        keyName = glasserparcellation[ROI][1]
    else:
        print("ROI not found: " + str(ROI) + " (e.g. use V1 for primary visual cortex or A1 for primary auditory cortex)")
        quit()





    jobs = []
    jobs.append("fslmaths <atlas_functionalspace> -thr <keyValue> -uthr <keyValue> <mask> -odt double")
    jobs.append("fslmaths <functional> -mas  <mask> <functional_ROI>")
    jobs.append("")

    mask = naming(atlas_functionalspace, "mask" + ROI, temp_path)
    functional_ROI = naming(functional, "ROI-" + ROI, output_path)

    for job in jobs:
        job = job.replace("<functional>", functional)
        job = job.replace("<keyValue>", str(keyValue))
        job = job.replace("<atlas_functionalspace>", atlas_functionalspace)
        job = job.replace("<mask>", mask)
        job = job.replace("<functional_ROI>", functional_ROI)

        run(job)

    return functional_ROI
    


def deobliqueUpsample (input, targetres):
    output = naming(input, "deobliqueUpsample", temp_path)

    job = "3dWarp -deoblique -NN -newgrid <targetres> -prefix <output> <input>"
    job = job.replace("<targetres>",    str(targetres)) # e.g., 0.3
    job = job.replace("<input>",        input)
    job = job.replace("<output>",       output)

    run(job)
    return output




def bet_functional (functional, anatomical_functionalspace_bet):
    # input: functional (can be 4D), anatomicalInFunctionalSpace
    # output: bet of functional (uses anatomical as reference)
    # typical use: functional_preprocessed = bet_functional (functional_current, T1w_functionalspace)

    checkFiles([functional, anatomical_functionalspace_bet])

    jobs = []
    jobs.append("fslmaths <anatomical_functionalspace_bet> -abs -bin <T1w_functionalspace_bet_mask>")
    jobs.append("fslmaths <functional> -mas  <anatomical_functionalspace_bet_mask> <functional_bet>")


    anatomical_functionalspace_bet          = output_path + "anatomical_functionalspace_bet" # naming(T1w_functionalspace, "bet", temp_path)
    functional_bet                          = naming(functional, "bet", output_path)
    anatomical_functionalspace_bet_mask     = naming(anatomical_functionalspace_bet, "mask", output_path)

    for job in jobs:
        job = job.replace("<functional>", functional)
        job = job.replace("<functional_bet>", functional_bet)
        job = job.replace("<anatomical_functionalspace_bet>", anatomical_functionalspace_bet)
        job = job.replace("<anatomical_functionalspace_bet_mask>", anatomical_functionalspace_bet_mask)

        run(job)

    return functional_bet



def createROImasksInFunctionalSpace (ROI, matrix_atlas_to_functional):
    print ("funktioniert noch nirt. brauch tnoch OIR mask von recon-all als input und centers für spheres.")
    # typical use: [ROIgyrus, ROIsphere, ROIcommon] = createROImasksInFunctionalSpace(ROI, matrix_atlas_to_functional)
    # output:
    #    ROIgyrus: mask for gyri at both lateral sides according to Glasser atlas (2 problems: a) typically larger than ROI for layerfmri b) not perfectly aligned with anatomical/functional data)
    #    ROIsphere: mask for sphere at one lateral side accoding to own settings based on the Glasser atlas (intended for better ROI selection for layerfmri (e.g., "motor" + "L" results in position on left postcentral gyrus where activity of thumb and index finger can be observed)
    #    ROIcommon: intersection of the above


    ROIlocations = {"ctx-lh-postcentral": [37, 11, 12],
                    "ctx-rh-postcentral": []}
    ROIxyz = ROIlocations[ROI] # x, y, z coordinates (according to own estimations and Glasser atlas - transformed into functional space)

    jobs = []


    # create sphere als POI in functional space
    # e.g., gyrus postcentralis: left: 37, 11, 12; right:  155, 119, 125
    jobs.append("fslmaths <atlas/>MNI152_2009_template.nii.gz -mul 0 -add 1 -roi <ROIx> 1 <ROIy> 1 <ROIz> 1 0 1 <temp/>POI.nii.gz -odt float")
    jobs.append("fslmaths <temp/>POI.nii.gz -kernel sphere 3 -fmean <temp/>ROIsphere.nii.gz -odt float")
    jobs.append("fslmaths <temp/>ROIsphere.nii.gz -kernel sphere 3.0 -fmean <temp/>ROIsphere.nii.gz -odt float")
    jobs.append("fslmaths <temp/>ROIsphere.nii.gz -abs -bin <temp/>ROIsphere.nii.gz")
    jobs.append("bash antsApplyTransforms -d 3 -i <temp/>ROIsphere.nii.gz -r <output/>anatomical_bet_functionalspace.nii.gz -o <output/>ROIsphere_functionalspace.nii.gz -n NearestNeighbor -t <matrix_atlas_to_functional>")

    # create intersection between ROI mask and sphere in functional space
    jobs.append("fslmaths <output/>mask_ROI_functionalspace.nii.gz -mul <output/>ROIsphere_functionalspace.nii.gz <output/>ROImask_common.nii.gz")

    for job in jobs:
        job = job.replace("<ROIx>", str(ROIxyz[0]))
        job = job.replace("<ROIy>", str(ROIxyz[1]))
        job = job.replace("<ROIz>", str(ROIxyz[2]))

        job = job.replace("<intensityLower>", str(intensityLower))
        job = job.replace("<intensityUpper>", str(intensityUpper))

        job = job.replace("<matrix_atlas_to_functional>", matrix_atlas_to_functional)

        job = job.replace("<atlas/>", atlas_path)
        job = job.replace("<temp/>", temp_path)
        job = job.replace("<output/>", output_path)

        run(job)


    return [output_path + "mask_ROI_functionalspace.nii.gz", output_path + "ROIsphere_functionalspace.nii.gz", output_path + "ROImask_common.nii.gz"]



def WM_GM_CSF_CTX_masks (subjectID, anatomical_functionalspace_bet, matrix_anatomical_to_functional):
    # typical use: [maskWM, maskGM, maskCSF, freesurfer_results, [listOfCortexMasks]] = WM_GM_CSF_masks (subjectID, anatomical_functionalspace_bet, matrix_anatomical_to_functional)
        # freesurfer_results is folder from recon-all, e.g., "/home/user/freesurfer/subject-574/"

    checkFiles ([anatomical_functionalspace_bet, matrix_anatomical_to_functional])

    import datetime
    current_datetime = datetime.datetime.now()
    myid = "subj_" + subjectID # + "_" + str(current_datetime.strftime("%Y-%m-%d_%H.%M"))



    maskWM = mask_path + "mask_wm.nii.gz"
    maskGM = mask_path + "mask_gm.nii.gz"
    maskWMGM = mask_path + "mask_wmgm.nii.gz"
    maskOuter = mask_path + "mask_csf.nii.gz"
    freesurfer_results = "~/freesurfer/" + myid + "/"





    jobs = []

    # get intersection-free masks in functional space for: WM, GM, CSF (outer)
    jobs.append("mri_convert <output/>anatomical_bet_functionalspace.nii.gz <temp/>anatomical_bet_functionalspace.mgh")
    jobs.append("mkdir ~/freesurfer") # location for results of recon-all (avoid problems with access rights by locating this within the home directory)
    jobs.append("rm -rf -f ~/freesurfer/<myid>") # neccessary otherwise recon-all will stop if directory exists
    jobs.append("recon-all -i <temp/>anatomical_bet_functionalspace.mgh -sd ~/freesurfer -subjid <myid> -all")
    jobs.append("mri_convert ~/freesurfer/myid/mri/aparc+aseg.mgz <temp/>aparc+aseg_functionalspace.nii.gz")

    # process masks
    jobs.append("mri_binarize --i <temp/>aparc+aseg_functionalspace.nii.gz --gm --o <maskGM>")
    jobs.append("mri_binarize --i <temp/>aparc+aseg_functionalspace.nii.gz --all-wm --o <maskWM>")

    # create outer mask: dilate GM and substract GM+WM: boarderGM-CSF
    jobs.append("fslmaths <maskWM> -max <maskGM> <temp/>mask_wmgm.nii.gz")
    jobs.append("fslmaths <temp/>mask_wmgm.nii.gz -sub 1 -mul -1 <temp/>invmask_wmgm.nii.gz")
    jobs.append("fslmaths <maskGM> -dilF -kernel sphere 1.0 -fmean <temp/>/gm_extended.nii.gz -odt float")
    jobs.append("fslmaths <temp/>gm_extended.nii.gz -abs -bin <temp/>gm_extended.nii.gz")
    jobs.append("fslmaths <temp/>gm_extended.nii.gz -mas <temp/>invmask_wmgm.nii.gz <maskOuter>")


    # cortex masks
    locations = {}
    locations["ctx-rh-lingual"]             = 2013  # primary visual
    locations["ctx-rh-cuneus"]              = 2005  # more visual
    locations["ctx-rh-pericalcarine"]       = 2021  # more visual
    locations["ctx-rh-lateraloccipital"]    = 2011  # more visual
    locations["ctx-rh-postcentral"]         = 2022  # primary motor
    locations["ctx-rh-superiortemporal"]    = 2030  # primary auditory

    for locR in list(locations.keys()): # create entries for other hemisphere
        locL = locR.replace("ctx-rh", "ctx-lh")
        locations[locL] = locations[locR] - 1000

    cortexMasks = [] # list of ctx masks, e.g., for later use
    for loc in list(locations.keys()):
        job = "fslmaths <temp/>aparc+aseg_functionalspace.nii.gz -thr <intensity> -uthr <intensity> <CTXmaskFile> -odt double"
        CTXmaskFile = mask_path + "mask_" + loc + ".nii.gz"
        job = job.replace("<CTXmaskFile>", CTXmaskFile)
        job = job.replace ("<loc>", loc)
        job = job.replace ("<intensity>", str(locations[loc]))
        jobs.append(job)
        cortexMasks.append(CTXmaskFile)




    for job in jobs:
        job = job.replace("<anatomical_functionalspace_bet>", anatomical_functionalspace_bet)
        job = job.replace("<matrix_anatomical_to_functional>", matrix_anatomical_to_functional)
        job = job.replace("<myid>", myid)
        job = job.replace("<maskGM>", maskGM)
        job = job.replace("<maskWM>", maskWM)
        job = job.replace("<maskOuter>", maskOuter)
        job = job.replace("<temp/>", temp_path)
        job = job.replace("<output/>", output_path)

        run(job)


    return [freesurfer_results, maskWM, maskGM, maskOuter, cortexMasks]


def upsampling (voxelsize_mm, anatomical, functional):
    # typical use: [anatomicalT1_up, functional_up, functional_reference, functional_reference_up] = upsampling(0.3, anatomicalT1, functional_current)

    checkFiles([anatomical, functional])

    isotropic = "{:.2f}".format(voxelsize_mm)
    jobs = []



    # create reference image: one image averaged from frist ten images in the time series and upscaled.
    functional_reference        = output_path + "functional_reference.nii.gz"
    functional_reference_up     = output_path + "functional_reference" + isotropic + ".nii.gz"

    jobs.append("rm <temp/>functional_first_10_volumes.nii.gz")  # fslroi cannot overwrite
    jobs.append("fslroi <functional> <temp/>functional_first_10_volumes.nii.gz 0 10")
    jobs.append("3dresample -dxyz 0.3 0.3 0.3 -prefix <temp/>functional_first_10_volumes_up.nii.gz -inset <temp/>functional_first_10_volumes.nii.gz")
    jobs.append("fslmaths <temp/>functional_first_10_volumes_up.nii.gz -Tmean <functional_reference_up>")
    jobs.append("fslmaths <temp/>functional_first_10_volumes.nii.gz -Tmean <functional_reference>")






    anatomical_up = output_path + "anatomical_up" + isotropic +".nii.gz" # naming(anatomical, "up"+isotropic, output_path)
    jobs.append("3dresample -dxyz <mm> <mm> <mm> -prefix <anatomical_up> -inset <anatomical>")

    if False: # !!! do not upsample function file. memory problems?? response is only "killed"
        functional_up = output_path + "functional_up" + isotropic +".nii.gz" # naming(functional, "up" + isotropic, output_path)
        jobs.append("3dresample -dxyz <mm> <mm> <mm> -prefix <functional_up> -inset <functional>")
    else:
        functional_up = functional # leave unchanged

        

    for job in jobs:
        job = job.replace("<anatomical>", anatomical)
        job = job.replace("<anatomical_up>", anatomical_up)
        job = job.replace("<functional>", functional)
        job = job.replace("<functional_reference>", functional_reference)
        job = job.replace("<functional_reference_up>", functional_reference_up)
        job = job.replace("<functional_up>", functional_up)
        job = job.replace("<mm>", isotropic)
        job = job.replace("<temp/>", temp_path)

        run (job)

    if not os.path.isfile(functional_up) or not os.path.isfile(anatomical_up) or not os.path.isfile(functional_reference):
        if False:
            print("Error: upsampling failed. At least one file is missing, please check. Memory may be too low.")
            quit()
        else:
            print("Warning: upsampling failed. At least one file is missing, please check. Memory may be too low.\n => Continue without upsampling")
            run('echo "Warning: upsampling failed. At least one file is missing, please check. Memory may be too low. Continue without upsampling"')
            return [anatomical, functional, functional_reference, functional_reference]


    return [anatomical_up, functional_up, functional_reference, functional_reference_up]






def checkFiles (filelist): # checks if all files in list exist
    error = False
    for file in filelist:
        if not os.path.isfile(file):
            error = True
            print ("### error: file not found: " + file)

    if error:
        quit()


def checkCommands(): # checks if all commands for the processing pipeline or a part thereof are available
    commands = ["fslmaths",
                "fsl_prepare_fieldmap",
                "fugue",
                "fslroi",
                "fsl-cluster",
                "fast",
                "3dvolreg",
                "3dresample",
                "flirt",
                "antsRegistrationSyN.sh",
                "antsApplyTransforms",
                "recon-all",
                "mri_binarize",
                "mri_convert",
                "mri_synthstrip",
                "LN2_LAYERS",
                "LN_INFO"]



    error = False
    for command in commands:

#        testing = "command -v <cmd> >> <file>"               # if command exist, folder name is returned, otherwise nothing.
        testing = "command -v <cmd>"
        testing = testing.replace("<cmd>", command)
        #result = subprocess.run(testing, shell=True, capture_output=True, text=True)
        #probe = result.stdout.strip()
        probe = os.popen(testing).read()


        #print ("got: " + probe)

        if len(probe) < 2:
            error = True
            print("### error: command not found: " + command)


    if error:
        print ("An error occured. Please make sure that the following are installed: FSL, ANTs, AFNI, freesurfer, laynii")
        quit() # an error occured: command not found, cannot run pipeline or part thereof





def createLayers (dict_label_filenames, numberLayers):
    # typical use: dict_layer_filenames = createLayers (dict_label_filenames, numberLayers)
    # dict_layer_filenames: {"ROI1": [path_equivol_metric, path_equidist_metric, path_equivol_layers, path_equidist_layers], "ROI2": "path2", [...], ...}

    checking = []
    for filename in dict_label_filenames.values():
        checking.append (filename)
    checkFiles(checking)

    dict_layer_filenames = {}
    jobs = []



    for ROIname, ROIlabelfile in dict_label_filenames.items():
        buffer = []

        buffer.append("echo ROI: " + ROIname)
        buffer.append("rm -rf -f <temp/>layniitemp")
        buffer.append("mkdir <temp/>layniitemp")
        buffer.append("cp <ROIlabelfile> <temp/>layniitemp/layniiLabels.nii.gz")                                  # perfom the LN2_LAYERS operation in a temporary folder because it creates many unneccessary files
        buffer.append("LN2_LAYERS -rim <temp/>layniitemp/layniiLabels.nii.gz -nr_layers <nr_layers> -equal_counts -equivol")    # equivol creates equivol layers in addition to equidist. resulting file is, for example, "[layniiLabels]_metric_equidist.nii.gz"

        output_files = []
        for type1 in ["layers", "metric", "layerbins"]:
            for type2 in ["equidist", "equivol"]:  # layerbins needed for profile

                output_file = "<layerresults/>result_" + str(numberLayers) + "-layers_" + ROIname + "_<type1>_<type2>.nii.gz"          # output: e.g., output/layer_results/layniiLabels_layerbins_equivol.nii.gz
                output_file = output_file.replace("<type1>", type1)
                output_file = output_file.replace("<type2>", type2)
                output_file = output_file.replace("<layerresults/>", layniiresult_path)
                output_files.append(output_file)


                job = "rm <output_file>"
                job = job.replace("<output_file>", output_file)
                buffer.append(job)


                job = "cp <temp/>layniitemp/layniiLabels_<type1>_<type2>.nii.gz <output_file>"
                job = job.replace("<output_file>", output_file)
                job = job.replace("<type1>", type1)
                job = job.replace("<type2>", type2)
                buffer.append(job)



        dict_layer_filenames[ROIname] = output_files


        buffer.append("rm -rf -f <temp/>layniitemp")


        for job in buffer:
            job = job.replace("<temp/>", temp_path)
            job = job.replace("<output/>", output_path)
            job = job.replace("<ROIlabelfile>", ROIlabelfile)
            job = job.replace("<nr_layers>", str(numberLayers))  # e.g., 9

            jobs.append(job)


    for job in jobs:
        run (job)

    return dict_layer_filenames


def createFunctionalReference (functional, upsampled_isotropic_voxelsize):
    # create functional reference image (intention: Functional image includes a time series of 300 samples. create average of the first 10 to obtain one image and remove noise. Input must be motion corrected.)
    # typical use:
    #       upsampled_isotropic_voxelsize = 0.3
    #       [functional_reference, functional_reference_up] = createFunctionalReference (functional_preprocessed, upsampled_isotropic_voxelsize)


    checkFiles ([functional_preprocessed])

    isotropic = str(upsampled_isotropic_voxelsize)
    functional_reference = output_path + "functional_reference.nii.gz"
    functional_reference_up = output_path + "functional_reference_up" + isotropic + ".nii.gz"

    jobs = []
    jobs.append("fslroi <functional> <temp/>functional_first_10_volumes.nii.gz 0 10")
    jobs.append("3dresample -dxyz 0.3 0.3 0.3 -prefix <temp/>functional_first_10_volumes_up.nii.gz -inset <temp/>functional_first_10_volumes.nii.gz")
    jobs.append("fslmaths <temp/>functional_first_10_volumes_up.nii.gz -Tmean <functional_reference_up>")
    jobs.append("fslmaths <temp/>functional_first_10_volumes.nii.gz -Tmean <functional_reference>")

    for job in jobs:

        job = job.replace("<functional_reference>", functional_reference)
        job = job.replace("<functional_reference_up>", functional_reference_up)
        job = job.replace("<functional>", functional)
        job = job.replace("<mm>", isotropic)
        job = job.replace("<temp/>", temp_path)

        run(job)


    return [functional_reference, functional_reference_up]



def create_ROI_labels(regions, cortexMasks, maskWM, maskGM, maskCSF):
    # creates a list of identiers (e.g. ctx-rh-lingual, ctx-rh-postcentral-sphere), and filenames of readily prepared masks. If possible, spheres are added to the list.
    # typical use: dict_label_filenames = create_ROI_labels(regions, cortexMasks, maskWM, maskGM, maskCSF)
    # mask_filenames can be used as input for LN2_LAYERS

    checking = [maskWM, maskGM, maskCSF]
    for ctx in cortexMasks:
        checking.append(ctx)
    checkFiles(checking)


    dict_label_filenames = {}


    for ROI in regions:

        ctx_mask = ""               # filename of the cortex masks as created by WM_GM_CSF_CTX_masks.
        for ctx in cortexMasks:      # find matching filename for region.
            if ROI in ctx:
                ctx_mask = ctx
        if ctx_mask == "":
            print ("Error: ctx mask not found. WM_GM_CSF_CTX_masks did not produce the results needed by createROIlabels.")
            quit()

        label_file = mask_path + "label_" + ROI + ".nii.gz"
        dict_label_filenames[ROI] = label_file # used as an output


        # output files
        layer_profile = output_path + "layerProfile.txt"
        layer_map = output_path + "layerMap.nii.gz"

        jobs = []
        jobs.append("fslmaths <ctx_mask> -add 0 <temp/>mask_gm_ROI.nii.gz")
        jobs.append("fslmaths <temp/>mask_gm_ROI.nii.gz -dilF -kernel sphere 2 -fmean <temp/>ROI_extended.nii.gz -odt float")
        jobs.append("fslmaths <temp/>ROI_extended.nii.gz -abs -bin <temp/>ROI_extended_bin.nii.gz")
        jobs.append("fslmaths <maskWM> -mas <temp/>ROI_extended_bin.nii.gz <temp/>mask_wm_ROI.nii.gz")
        jobs.append("fslmaths <maskCSF> -mas <temp/>ROI_extended_bin.nii.gz <temp/>mask_csf_ROI.nii.gz")
        jobs.append("fslmaths <temp/>mask_gm_ROI.nii.gz -abs -bin <temp/>mask_gm_ROI_bin.nii.gz")
        jobs.append("fslmaths <temp/>mask_gm_ROI_bin.nii.gz -mul 2 <temp/>mask_gm_ROI_bin_2.nii.gz")
        jobs.append("fslmaths <temp/>mask_wm_ROI.nii.gz -mul 3 <temp/>mask_wm_ROI_3.nii.gz")
        jobs.append("fslmaths <temp/>mask_csf_ROI.nii.gz -add <temp/>mask_gm_ROI_bin_2.nii.gz -add <temp/>mask_wm_ROI_3.nii.gz <label_file>")



        for job in jobs:
            job = job.replace("<ctx_mask>", ctx_mask)
            job = job.replace("<temp/>", temp_path)
            job = job.replace("<output/>", output_path)
            job = job.replace("<maskWM>", maskWM)
            job = job.replace("<maskGM>", maskGM)
            job = job.replace("<maskCSF>", maskCSF)
            job = job.replace("<label_file>", label_file)


            run(job)


    return dict_label_filenames


def add_ROI_spheres(dict_label_filenames, functional_reference, matrixMNItoFunctional):
    # typical use dict_label_filenames2 = add_ROI_spheres(dict_label_filenames, functional_reference, matrixMNItoFunctional)
    # if possible (e.g., for the postcentral gyrus), spherical ROI masks are added to the existing set of masks.

    result = {} # dictionary: {ROIname1:labelfile1, ROIname2:labelfile2}
    jobs = []
    checking = []


    #sphereDefinitions = {"ctx-rh-postcentral":[155, 119, 125, 8], "ctx-lh-postcentral":[37, 11, 12, 8]} # list contains x, y, z, radius in MNI space

    sphereDefinitions = {"ctx-rh-postcentral": [145, 110, 117, 15],
                         "ctx-lh-postcentral": [45, 110, 117, 15]}  # list contains x, y, z, radius in MNI space. Values refer to "MNI152_2009_template.nii.gz"

    sphereImage = temp_path + "sphere.nii.gz" # temporary, will be reused for different ROIs
    sphere_functionalspace = temp_path + "sphere_functionalspace.nii.gz"

    for key, value in dict_label_filenames.items():
        checking.append(value)
        result[key] = value  # add original

        if key in sphereDefinitions:

            [x, y, z, radius] = sphereDefinitions[key]

            labelfile_original = value
            labelfile_sphere = labelfile_original.replace(".nii.gz", "_ROIsphere.nii.gz")

            result[key + "_sphere"] = labelfile_sphere

            buffer = []
            buffer.append("fslmaths <MNI> -mul 0 -add 1 -roi <x> 1 <y> 1 <z> 1 0 1 <sphere> -odt float")
            sum = 0
            while sum < radius: # use multiple steps: sphere is less fine but calculation becomes faster
                buffer.append("fslmaths <sphere> -kernel sphere 5 -fmean <sphere> -odt float")
                sum += 5
            buffer.append("fslmaths <sphere> -abs -bin <sphere>")
            buffer.append("antsApplyTransforms -d 3 -i <sphere> -r <labelfile_original> -o <sphere_functionalspace> -n NearestNeighbor -t <matrixMNItoFunctional>") # ! reference must be labelfile_original (is in functional space and has the right dimensions and voxel size)
            buffer.append("fslmaths <labelfile_original> -mas <sphere_functionalspace> <labelfile_sphere>")


            for job in buffer:
                job = job.replace("<x>", str(x))
                job = job.replace("<y>", str(y))
                job = job.replace("<z>", str(z))
                job = job.replace("<radius>", str(radius))
                job = job.replace("<sphere>", sphereImage)
                job = job.replace("<labelfile_original>", labelfile_original)
                job = job.replace("<labelfile_sphere>", labelfile_sphere)
                job = job.replace("<MNI>", constant_MNItemplate)
                job = job.replace("<functional_reference>", functional_reference)
                job = job.replace("<matrixMNItoFunctional>", matrixMNItoFunctional)
                job = job.replace("<sphere_functionalspace>", sphere_functionalspace)

                jobs.append(job)


    checkFiles(checking)

    for job in jobs:
        run (job)

    return result



def calculate_layer_profile(dict_LN2results, activationMap):
    # takes # layerbin_equidist from file list in dictionary LN2results[ROI]
    # typical use: dict_LN2results = calculate_layer_profile(dict_LN2results, activationMap) # adds .txt file for layer profile and corresponding PNG plot.

    jobs = []
    for ROIname, layerfiles in dict_LN2results.items():
        for ROIfile in layerfiles:
            if "layerbins" in ROIfile and "equivol" in ROIfile: # different files are associated with an ROIname
                print("\n\n" + ROIfile + "\n")
                buffer = []
                buffer.append("echo ROI: " + ROIname)
                buffer.append("rm <layer_profile.txt>")
                buffer.append("LN2_PROFILE -input <activationMap.nii.gz> -layers <layerbins_equivol.nii.gz> -plot -output <layer_profile.txt>") # to extract layer-profiles from any activity map based on layer masks

                txt_output = layniiresult_path + ROIname + ".txt"
                layerfiles.append(txt_output)
                dict_LN2results[ROIname] = layerfiles               # Add txt file to list. Modified dictionary will be returned to provide an overview.


                for job in buffer:
                    job = job.replace("<activationMap.nii.gz>", activationMap)
                    job = job.replace("<layerbins_equivol.nii.gz>", ROIfile)
                    job = job.replace("<layer_profile.txt>", txt_output)
                    jobs.append(job)


    for job in jobs:
        run(job)

    return dict_LN2results




def create_activation_map(functional_preprocessed, functional_reference):    # currently just a placeholder, creates uniform noise
    # typical use: activation_map = create_activation_map(functional_preprocessed, functional_reference_up)

    activation_map = naming(functional_preprocessed, "activation", output_path)
    jobs = []
    jobs.append("fslmaths <functional_reference> -mul 0 -rand <activation_map>") # uniform noise in functional space

    for job in jobs:
        job = job.replace("<functional_reference>", functional_reference)
        job = job.replace("<activation_map>", activation_map)
        run (job)


    return activation_map





for job in jobs:
    [subjectID, anatomicalT1, anatomicalT2, fieldmap_mag, fieldmap_phase, functional] = job

    with open("run.txt", "a+") as file:
        file.write("\n\n\n" + subjectID + "\n")



    ################################################
    # check if all files and commands are available that are needed in the pipeline
    checkfiles = []

    for file in job[1:]:
        checkfiles.append(file)
        if "func/" in file:
            checkfiles.append(file.replace(".nii.gz", ".json"))

    checkFiles(checkfiles)
    checkCommands()
    print("All checks passed: Required input files and functions available.")
    ################################################


    global_simulated = True

    functional_current = functional

    # 1. undistortion (i.e., correct warping from bias field)
    functional_current =  preprocessing_undistortion (functional_current, fieldmap_mag, fieldmap_phase)

    # 3. bias field correction (Bias field correction is performed before motion correction because the bias field is independent of motion.)
    functional_current = preprocessing_biasfield4D(functional_current)

    # 4. motion correction to distortion corrrected functional images (with additional output of transformation matrix):
    functional_current = preprocessing_motioncorrection(functional_current)




    # Preprocessing done.
    functional_preprocessed = functional_current





    # 5. create functional reference image (intention: Functional image includes a time series of 300 samples. create average of the first 10 to obtain one image and remove noise. Input must be motion corrected.)
    upsampled_isotropic_voxelsize = 0.3
    [functional_reference, functional_reference_up] = createFunctionalReference (functional_preprocessed, upsampled_isotropic_voxelsize)


    # 6. coregistrations
    # example:   input:
    #               - functional_reference
    #               - functional_reference_up (previously upsampled)
    #               - anatomical T1
    #               - MNI template
    #               - Glasser atlas
    #
    #           output:
    #               - anatomical in functional space (upsampled)
    #               - Glasser atlas in functional space (upsampled)
    #               - anatomical_bet_functional (upsampled)
    #               - matrix_atlas_to_functional
    #               - matrix_anatomical_to_functional
    [T1w_functionalspace, atlas_functionalspace, T1w_functionalspace_bet, matrix_atlas_to_functional, matrix_anatomical_to_functional] = coregistrations (functional_preprocessed, functional_reference, functional_reference_up, anatomicalT1, constant_MNItemplate, constant_atlas)



    # 7. brain extraction: uses structural image (already in functional space) for brain extraction of functional images. Result is preprocessed fMRI (undistortion, motion correction, bias field removed, brain extracted)
    functional_preprocessed_bet = bet_functional (functional_preprocessed, T1w_functionalspace_bet)
 
    # 8. obtain intersection-free WM, GM, CSF masks (freesurfer_results  is folder, e.g., "/home/user/freesurfer/myid/")
    [freesurfer_results, maskWM, maskGM, maskCSF, cortexMasks] = WM_GM_CSF_CTX_masks(subjectID, T1w_functionalspace_bet, matrix_anatomical_to_functional) # uses recon-all

    ROIs = ["ctx-rh-lingual", "ctx-rh-cuneus", "ctx-rh-pericalcarine", "ctx-rh-lateraloccipital", "ctx-rh-postcentral", "ctx-rh-superiortemporal",
            "ctx-lh-lingual", "ctx-lh-cuneus", "ctx-lh-pericalcarine", "ctx-lh-lateraloccipital", "ctx-lh-postcentral", "ctx-lh-superiortemporal"]


    dict_label_filenames = create_ROI_labels(ROIs, cortexMasks, maskWM, maskGM, maskCSF) # creates a list of idenfiers (e.g. ctx-rh-lingual, ctx-rh-postcentral-sphere), and filenames of readily prepared masks.
    dict_label_filenames2 = add_ROI_spheres(dict_label_filenames, functional_reference, matrix_atlas_to_functional) # if possible (e.g., for the postcentral gyrus), spherical ROI masks are added to the existing set of masks.


    # 9. create activation map (FSL FEAT?)
    activation_map = create_activation_map(functional_preprocessed, functional_reference_up)  # just a placeholder, creates only white noise currently.


    # 10. generate layers
    numberLayers = 7
    LN2results = createLayers (dict_label_filenames2, numberLayers) # uses LN2_LAYERS
    # dict_layer_filenames: {"ROI1": [path_equivol_metric, path_equidist_metric, path_equivol_layers, path_equidist_layers], "ROI2": "path2", [...], ...}



    # 11. calculate layer profile: "LN2_PROFILE -input <activationMap_3D.nii.gz> -layers <layerbins_equivol.nii.gz> -plot -output <layer_profile.txt>"
    LN2results = calculate_layer_profile(LN2results, activation_map) # generates .txt and .png files associated with layer profile for each ROI in the dictionary in addition to the .nii.gz files created by LN2_LAYERS in a previous step.

