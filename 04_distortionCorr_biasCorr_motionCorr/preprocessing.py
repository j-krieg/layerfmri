import subprocess
import os
import shutil
import nibabel as nib
import numpy as np
import json
import time

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
    subprocess.Popen(task.split(), stdout=subprocess.PIPE)

    

anatomical_path     = "data/sub-s574/anat/"
fieldmap_path       = "data/sub-s574/fmap/"
functional_path     = "data/sub-s574/func/"



anatomicalT1         = anatomical_path  + "8-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_UNI_DEN.nii"
anatomicalT2         = anatomical_path  + "sub-s574_acq-grappatest_T2w.nii.gz"
jobs = []

if True:
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
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()


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


# phase unwrapping with prelude
def run_prelude(input_img, output_img):
    if False:
        task = "prelude -c <input_img> -o <output_img> -m temp/mask.nii.gz"
        task = task.replace("<input_img>", input_img)
        task = task.replace("<output_img>", output_img)


        process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()

    else:
        subprocess.check_output(['prelude', '-c', input_img, '-o', output_img])
    # usage: prelude -c <rawphase> -o <unwrappedphase> [options]
    #print ("Done: Unwrapped fieldmap_phase with FSL prelude.")





def run_mask(input_img_bet, threshold, output_mask):
    # creates a mask from an imput image (typically, bet has already been performed).

    task = "fslmaths <input_img> -thr <threshold> -bin <output_mask>"
    task = task.replace("<input_img>",          input_img_bet)
    task = task.replace("<threshold> ",         str(threshold))
    task = task.replace("<output_mask>",        output_mask)

    process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()



def preprocessing_prepareFieldmapMag (fieldmap_mag, input_img, fieldmap_mag_prepared):
    # this comprises the steps of
    #   - bias correcting the magnitude image for better skullstripping results
    #   - skullstripping
    #   - eroding
    #   - registration to the input image
    # outout: fieldmap_mag_prepared)
    print("- fieldmap magnitude preparation")


    #  1. bias correcting the magnitude image for better skullstripping results
    fieldmap_mag_step1 = temp_path + os.path.basename(fieldmap_mag.replace(".nii", "_step1.nii"))
    task= "fsl_anat --strongbias --nocrop --noreg --nosubcortseg --noseg -i <input> -o <output>"
    task = task.replace("<input>", fieldmap_mag)
    task = task.replace("<output>", fieldmap_mag_step1)

    process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    print("   - mag img: bias correction done.")


    outputfile = fieldmap_mag_step1 + ".anat/T1_biascorr.nii.gz"
    if os.path.exists(outputfile):
        shutil.copy2(outputfile, fieldmap_mag_step1) # output file is temp/[17-gre_field_mapping_HCP_0pp8mm3_e2_step1.nii].anat/T1_biascorr.nii.gz
    else:
        print("Unable to perform bias correction on magnitude image.")
        shutil.copy2(fieldmap_mag, fieldmap_mag_step1)



    # 2. skullstripping
    fieldmap_mag_step2 = fieldmap_mag_step1.replace("_step1.nii", "_step2.nii")
    task = "bet <input> <output> -R"
    task = task.replace("<input>", fieldmap_mag_step1)
    task = task.replace("<output>", fieldmap_mag_step2)

    process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    print("   - mag img: skullstripping done.")


    # 2b create mask
    run_mask(fieldmap_mag_step2, 0.5, "temp/mask.nii")

    # 3. eroding
    fieldmap_mag_step3 = fieldmap_mag_step2.replace("_step2.nii", "_step3.nii")
    task = "fslmaths <input> -ero <output>"
    task = task.replace("<input>", fieldmap_mag_step2)
    task = task.replace("<output>", fieldmap_mag_step3)

    process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    print("   - mag img: eroding done.")


    # 4. registration to input image
    run_reg(fieldmap_mag_step3, input_img, fieldmap_mag_prepared)
    print("   - mag img: registration done.")



def preprocessing_prepareFieldmapPhase (fieldmap_phase, input_img, fieldmap_phase_prepared):
    # this comprises the steps of
    #   - convert phase to radians
    #   - unwrapping the phase image
    #   - registration to the input image
    # output: fieldmap_phase_prepared
    print("- fieldmap phase preparation")

    # 0 brain extraction to make phase unwrapping work better
    # "The reason is that they are not brain extracted, and the non-brain material makes all the unwrapping massively more complicated.  So if you brain extract your images to start with (apply BET to the magnitude and then use this mask with fslmaths on the phase images). " (source: https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;b7eac096.1310)
    if False:
        fieldmap_phase_bet = temp_path + os.path.basename(fieldmap_phase.replace(".nii", "_bet.nii"))
        task = "bet <input> <output> -R"
        task = task.replace("<input>", fieldmap_phase)
        task = task.replace("<output>", fieldmap_phase_bet)

        process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()






    # 1. Convert phasemap to radians
    # -4096 to + 4096: output percentils with fslstats img.nii -p 0 and -p 100




    fieldmap_phase_rad = temp_path + os.path.basename(fieldmap_phase.replace(".nii", "_rad.nii"))
    #task = "fslmaths <input> -div 2048 -sub 1 -mul 3.14159 <output> -odt float" # 0 to 2048 fslstats # -4096 to +4096

    #Siemens: -4096 to +4096, needs this factor to bringt to -2Pi to +2Pi: 651.8986

    # task = "fslmaths <input> -div 651.8986 <output> -odt float"

    task = "fslmaths <input> -div 4096 -mul 3.14159 <output> -odt float"  # -4096 to 4096 to -Pi to +Pi
    task = task.replace("<input>", fieldmap_phase)
    task = task.replace("<output>", fieldmap_phase_rad)

    process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

 #   print ("\n in: " + fieldmap_phase)
 #   print("\n out: " + fieldmap_phase_rad)
    print("   - ph img: phase converted to radians.")



    # 2. phase unwrapping for fieldmap
    fieldmap_phase_unwrapped = temp_path + os.path.basename(fieldmap_phase.replace(".nii", "_unwrapped.nii"))
    run_prelude(fieldmap_phase_rad, fieldmap_phase_unwrapped)
    print("   - ph img: unwrapping done.")

    # 3. registration to input image
    run_reg(fieldmap_phase_unwrapped, input_img, fieldmap_phase_prepared)
    print("   - ph img: registration done.")



    # create a fieldmap with registration to the input image as a distortion estimation
def preprocessing_createFieldmap(input_img, fieldmap_mag, fieldmap_phase, fieldmap_distortion): # uses fsl_prepare_fieldmap
    # create a fieldmap with registration to the input image
    # input_img:            EPI
    # fieldmap_mag:         magnitude image
    # fieldmap_phase        phase image
    # fieldmap_distortion   output that is generated based on the three images above with registration to input_img



    # 1. prepare fieldmap_mag:
    fieldmap_mag_prepared = temp_path + os.path.basename(fieldmap_mag.replace(".nii", "_prepared.nii"))

 #!!!!!!!!!!!!!   preprocessing_prepareFieldmapMag (fieldmap_mag, input_img, fieldmap_mag_prepared)
    # this comprises the steps of
    #   - bias correcting the magnitude image for better skullstripping results
    #   - skullstripping
    #   - eroding
    #   - registration to the input image



    # 2. prepare fieldmap_phase:
    if False:
        fieldmap_phase_prepared = temp_path + os.path.basename(fieldmap_phase.replace(".nii", "_prepared.nii"))
        preprocessing_prepareFieldmapPhase (fieldmap_phase, input_img, fieldmap_phase_prepared)
        # this comprises the steps of
        #   - convert to radians
        #   - unwrapping the phase image
        #   - registration to the input image
    else:
        fieldmap_phase_prepared = fieldmap_phase # actually none of the above steps required and somt part of preprocessing_prepareFieldmapPhase produces an error (manually review images)
        print ("skipped phasemap preparation.")




    # 3. Generate a fieldmap based on the prepared magnitude and phase images:
    print("- distortion estimation...")
    task = "fsl_prepare_fieldmap <scanner> <phase_image> <magnitude_image> <out_image> <deltaTE(in ms)>"
    task += " --nocheck" # currently problem: "Phase image values do not have expected range. Expecting at least 90% of 0 to 4096, but found -4.808575 to 11.701118"
    task = task.replace("<scanner>", "SIEMENS")
#!!!!    task = task.replace("<phase_image>", fieldmap_phase_prepared)
    task = task.replace("<phase_image>", fieldmap_phase) #  page 31 says, i automatically does theunwrapping for Siemens scanners: https://gate.nmr.mgh.harvard.edu/wiki/whynhow/images/8/8f/WhyNhow_FSL_final-WEB.pdf
    task = task.replace("<magnitude_image>", fieldmap_mag_prepared)
    task = task.replace("<out_image>", fieldmap_distortion)
    task = task.replace("<deltaTE(in ms)>", "1.02") # (defaults are *usually* 2.46ms on SIEMENS)



    process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    #print(" ok")

    #print("Done: Fieldmap generation based on phase image and magnitude image with registration to input EPI.")



def fieldmapSimplified (fieldmap_mag, fieldmap_phase, fieldmap_out):
    tasks = ["bet <fieldmap_mag> <fieldmap_mag_bet>",
            "fslmaths <fieldmap_mag_bet> -ero <fieldmap_mag_bet_ero>",
            "fsl_prepare_fieldmap SIEMENS <fieldmap_phase> <fieldmap_mag_bet_ero> <fieldmap_out> 1.02"]

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




def preprocessing_applyFieldmap(input_img, fieldmap_distortion, output_corrected, matrix_distortionCorr):
    # uses fugue(FMRIB's Utility for Geometrically Unwarping EPIs) to perform an unwarping of an EPI image based on fieldmap data.

    dwell_time = get_dwell_time(input_img) # taken from json file (alternatively, DICOM header)

    print("- fieldmap application for undistortion...")
    task = "fugue -i <input_epi> --dwell=<dwell_time> --loadfmap=<fieldmap> -u <result>" # "-savematrix <matrix_distortionCorr>"
    #task += " --smooth3=1" # Gaussian smoothing with sigma in mm (optional)
    task = task.replace("<input_epi>",                  input_img)
    task = task.replace("<dwell_time>",                 str(dwell_time))
    task = task.replace("<fieldmap>",                   fieldmap_distortion)
    task = task.replace("<result>",                     output_corrected)
    #task = task.replace("<matrix_distortionCorr>",      matrix_distortionCorr)
    task += " --verbose"


    process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    #print(" ok")

    #print("Done: Application of fieldmap to epi for distortion correction.")


def preprocessing_motionCorrection (epi_input, epi_motionCorrected, matrix_output):
# saved the motion-corrected .nii image and the matrix that was used to perfrom the transformation
    #task = "3dvolreg -twopass -1Dmatrix_save <matrix_output> -prefix <epi_motionCorrected> <epi_input>"
    task = "3dvolreg -twopass -prefix <epi_motionCorrected> <epi_input>"

    task = task.replace("<epi_input>", epi_input)
    task = task.replace("<matrix_output>", matrix_output)
    task = task.replace("<epi_motionCorrected>", epi_motionCorrected)
    process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()



def biasfieldCorrection (input, output):
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

    # alternative 0: not working, does not provide a bias field, just segmentations
    # option fast -N returns the bias field and does not calculate a correction (needed if working with a time series because the same field should be applied to all images for correction because we are interested in functional data)
    #tasks = ["fast -N -o <biasfield> <input>", "fslmaths <input> -div <biasfield> <output>"]

    # alternative 1:
    # estimates the bias filed but only provides the corrected image as an output (problem: removes everything but the first image from a time series in functional imaging)
    #restore = biasfield.replace(".nii", "_restore.nii")    #   automatically gests this prefix from fast.
    #tasks = ["fast -B -o <biasfield> <input>", "rm <output>", "mv <restore> <output>"]


    # alternative 2:
    # problem with this is, that the output file ends with "bold_biasfield_restore.nii.gz" and option "B" already outputs corrected image
    #tasks = ["fast -B -o <biasfield> <input>", "fslmaths <input> -div <biasfield> <output>"]
    # Warning: An input intended to be a single 3D volume has multiple timepoints. Input will be truncated to first volume, but this functionality is deprecated and will be removed in a future release.
    # proposal: derive bias field from first and apply to all
    #tasks = ["fast -B -o <biasfield> <input>", "fslmaths <input> -div <biasfield> <output>"]


    # alternative 3:
    #tasks = ["fsl_anat --strongbias --nocrop --noreg --nosubcortseg --noseg -i <input> -o <output>"]
    # works only on a 3D volume, not on a 4D volume.


    for task in tasks:
        task = task.replace("<input>", input)
        task = task.replace("<output>", output)
        task = task.replace("<biasfield>", biasfield)
        task = task.replace("<biasfield_smooth>", biasfield_smooth)
 #       task = task.replace("<restore>", restore)
        task = task.replace("<firstimg>", firstimg)
        task = task.replace("<firstbiascorr>", firstbiascorr)
        task = task.replace("<firstbiascorr_restore>", firstbiascorr_restore)

        #print (task)
        run(task)






def naming (source_file, suffix, target_directory):
    # input:    source_file: directory + filename to original file
    #           suffix: new suffix e.g., "biascorr" to name the new file "_biascorr.nii"
    #           target_directory: replace directory in path

    if not target_directory.endswith("/"):
        target_directory += "/"

    return target_directory + os.path.basename(source_file.replace(".nii", ("_" + suffix + ".nii")))


def biasfieldUndistortion(input, fieldmap_mag, fieldmap_phase, output):
    #   input: functional image
    #   output: undistorted functional image

    fieldmap_distortion = naming (input, "biasfieldEstimation", temp_path)
    fieldmapSimplified(fieldmap_mag, fieldmap_phase, fieldmap_distortion)
    # uses fsl_prepare_fieldmap (and bet for skullstripping, fslmaths -ero for eroding, fsl_anat --strongbias for bias correction, epi-reg for registration, fslmaths -div/sub/mul for radian conversion)

    # 2. apply the fieldmap to undistort the input image
    matrix_distortionCorr = temp_path + os.path.basename(
        (output.replace(".nii", "_matrixDistortion.1D")).replace(".gz", ""))  # does not produce matrix as output

    preprocessing_applyFieldmap(input, fieldmap_distortion, output, matrix_distortionCorr)
    # uses fugue



for job in jobs:
    [anatomicalT1, anatomicalT2, fieldmap_mag, fieldmap_phase, functional] = job

    # 1. bias field undistortion
    functional_undistorted = naming(functional, "undistorted", bias_path)
    biasfieldUndistortion (functional, fieldmap_mag, fieldmap_phase, functional_undistorted)



    # 2. bias field correction
    T1_biascorr = naming(anatomicalT1, "biascorr", bias_path)
    T2_biascorr = naming(anatomicalT2, "biascorr", bias_path)
    functional_biascorr = naming(functional_undistorted, "biascorr", bias_path)



    biasfieldCorrection(functional_undistorted, functional_biascorr)
    #biasfieldCorrection(anatomicalT1, T1_biascorr)
    #biasfieldCorrection(anatomicalT2, T2_biascorr)



    # 3. motion correction to distortion corrrected functional images (with additional output of transformation matrix):
    matrix_motion = (naming(functional_biascorr, "matrixMotion", temp_path)).replace(".nii", ".1D").replace(".gz", "")
    functional_motioncorr = naming (functional_biascorr, "motioncorr", output_path)
    preprocessing_motionCorrection(functional_biascorr, functional_motioncorr, matrix_motion)
    # uses 3dvolreg

