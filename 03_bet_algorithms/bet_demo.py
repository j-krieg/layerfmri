import subprocess
import os
import shutil
import nibabel as nib
import numpy as np
import json
import time



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




betfolder = "bet_demo/"
run ("rm -rf " + betfolder)
run("mkdir " + betfolder)

T1 = "data/sub-s574/anat/8-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_UNI_DEN.nii"
T2 = "data/sub-s574/anat/sub-s574_acq-grappatest_T2w.nii.gz"
bet_demo (T1, T2, betfolder+"nobiascorrT1-bet.nii.gz")


T1biascorr = "processed/T1_biascorr.nii.gz"
if False: # already performed before, time intensive
    task= "fsl_anat --strongbias --nocrop --noreg --nosubcortseg --noseg -i <input> -o <output>"    # subfolder: bet_demo/t1_biascorrected.nii.gz.anat/T1_biascorr.nii.gz
    task = task.replace ("<input>", T1)
    task = task.replace ("<output>", T1biascorr)
    run(task)

bet_demo (T1biascorr, T2, "data/biascorrectedT1-bet.nii.gz")