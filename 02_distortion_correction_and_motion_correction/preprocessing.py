import subprocess
import os
import shutil
import nibabel as nib
import numpy as np
import json

output_path       = "output/"
temp_path         = "temp/"
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


anatomical          = fieldmap_path + "sub-s574_t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_UNI_DEN.nii"
jobs = []


fieldmap_mag        = fieldmap_path + "sub-s574_acq-0p8mm_magnitude2.nii.gz"
fieldmap_phase      = fieldmap_path + "sub-s574_acq-0p8mm_phasediff.nii.gz"
epi                 = functional_path +"sub-s574_task-rest_acq-0p8mm_bold.nii.gz"
jobs.append([anatomical, fieldmap_mag, fieldmap_phase, epi])


fieldmap_mag        = fieldmap_path + "sub-s574_acq-1p3mm_magnitude2.nii.gz"
fieldmap_phase      = fieldmap_path + "sub-s574_acq-1p3mm_phasediff.nii.gz"
epi                 = functional_path +"sub-s574_task-rest_acq-1p3mm_bold.nii.gz"
jobs.append([anatomical, fieldmap_mag, fieldmap_phase, epi])









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
        quit()

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
    subprocess.check_output(['prelude', '-c', input_img, '-o', output_img])
    # usage: prelude -c <rawphase> -o <unwrappedphase> [options]
    #print ("Done: Unwrapped fieldmap_phase with FSL prelude.")





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

    # 1. Convert phasemap to radians
    fieldmap_phase_rad = temp_path + os.path.basename(fieldmap_phase.replace(".nii", "_rad.nii"))
    task = "fslmaths <input> -div 2048 -sub 1 -mul 3.14159 <output> -odt float"
    task = task.replace("<input>", fieldmap_phase)
    task = task.replace("<output>", fieldmap_phase_rad)

    process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
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
    preprocessing_prepareFieldmapMag (fieldmap_mag, input_img, fieldmap_mag_prepared)
    # this comprises the steps of
    #   - bias correcting the magnitude image for better skullstripping results
    #   - skullstripping
    #   - eroding
    #   - registration to the input image



    # 2. prepare fieldmap_phase:
    if True:
        fieldmap_phase_prepared = temp_path + os.path.basename(fieldmap_phase.replace(".nii", "_prepared.nii"))
        preprocessing_prepareFieldmapPhase (fieldmap_phase, input_img, fieldmap_phase_prepared)
        # this comprises the steps of
        #   - convert to radians
        #   - unwrapping the phase image
        #   - registration to the input image
    else:
        fieldmap_phase_prepared = fieldmap_phase # actually none of the above steps required and somt part of preprocessing_prepareFieldmapPhase produces an error (manually review images)




    # 3. Generate a fieldmap based on the prepared magnitude and phase images:
    print("- distortion estimation...")
    task = "fsl_prepare_fieldmap <scanner> <phase_image> <magnitude_image> <out_image> <deltaTE(in ms)>"
    task += " --nocheck" # currently problem: "Phase image values do not have expected range. Expecting at least 90% of 0 to 4096, but found -4.808575 to 11.701118"
    task = task.replace("<scanner>", "SIEMENS")
    task = task.replace("<phase_image>", fieldmap_phase_prepared)
    task = task.replace("<magnitude_image>", fieldmap_mag_prepared)
    task = task.replace("<out_image>", fieldmap_distortion)
    task = task.replace("<deltaTE(in ms)>", "1.02") # (defaults are *usually* 2.46ms on SIEMENS)



    process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    #print(" ok")

    #print("Done: Fieldmap generation based on phase image and magnitude image with registration to input EPI.")




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
    task = "3dvolreg -twopass -1Dmatrix_save <matrix_output> -prefix <epi_motionCorrected> <epi_input>"

    task = task.replace("<epi_input>", epi_input)
    task = task.replace("<matrix_output>", matrix_output)
    task = task.replace("<epi_motionCorrected>", epi_motionCorrected)
    process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()




def undistort_field_and_motion (anatomical, epi_input, matrix_field, matrix_motion, epi_output): # uses cat_matvec and 3dAllineate

    if not os.path.isfile(matrix_motion):
        print(f"The file '{matrix_motion}' not found.")
        quit()

    if not os.path.isfile(matrix_field):
        print(f"The file '{matrix_field}' not found.")
        quit()

    concatenated_matrix = temp_path + ((os.path.basename(epi_input)).replace(".nii", "_concatenatedMatrix.1D")).replace(".gz", "")
    task = "cat_matvec -ONELINE anat_ns+tlrc::WARP_DATA -I <matrix_distortion> -I <matrix_motion> > <concatenated_matrix>"
    task = task.replace("<concatenated_matrix>", concatenated_matrix)
    task = task.replace("<matrix_distortion>", matrix_field)
    task = task.replace("<matrix_motion>", matrix_motion)
    process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    task = "3dAllineate -base <anatomical> -1Dmatrix_apply <concentrated_matrix> -mast_dxyz 3 -float -prefix <epi_output> <epi_input>"
    task = task.replace("<concatenated_matrix>", concatenated_matrix)
    task = task.replace("<epi_input>", epi_input)
    task = task.replace("<epi_output>", epi_output)
    task = task.replace("<anatomical>", anatomical)
    process = subprocess.Popen(task.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

for job in jobs:
    [anatomical, fieldmap_mag, fieldmap_phase, epi] = job

    epi_corrected        = output_path + os.path.basename(epi.replace(".nii", "_distortionCorrected.nii"))

    print("\n\n" + epi + ":")

    # 1. create a fieldmap with registration to the input image
    fieldmap_distortion = output_path + os.path.basename(epi.replace(".nii", "_distortion.nii"))
    preprocessing_createFieldmap(epi, fieldmap_mag, fieldmap_phase, fieldmap_distortion)
    # uses fsl_prepare_fieldmap (and bet for skullstripping, fslmaths -ero for eroding, fsl_anat --strongbias for bias correction, epi-reg for registration, fslmaths -div/sub/mul for radian conversion)

    # 2. apply the fieldmap to undistort the input image
    matrix_distortionCorr = temp_path + os.path.basename((epi_corrected.replace(".nii", "_matrixDistortion.1D")).replace(".gz", "")) # does not produce matrix as output
    preprocessing_applyFieldmap(epi, fieldmap_distortion, epi_corrected, matrix_distortionCorr)
    # uses fugue

    # 3. motion correction to distortion corrrected functional images (with additional output of transformation matrix):
    matrix_motion = temp_path + os.path.basename((epi_corrected.replace(".gz", "")).replace(".nii", "_matrixMotion.1D"))
    epi_motion_corrected = epi_corrected.replace(".nii", "_motionCorr.nii")
    preprocessing_motionCorrection(epi_corrected, epi_motion_corrected, matrix_motion)
    # uses 3dvolreg

    # 4. Undistort with concatenated matrix to perform only one transformation instead of two
    epi_motion_and_distortion_corrected = output_path + os.path.basename(epi.replace(".nii", "_motion_and_distortion_corrected_one_step.nii"))
    undistort_field_and_motion (anatomical, epi, matrix_distortionCorr , matrix_motion, epi_motion_and_distortion_corrected)
    # uses cat_matvec and 3dAllineate
