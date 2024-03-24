import os
import re
import sys
import io
import math
import time
import datetime
import shutil
import numpy as np
import matplotlib
matplotlib.use("pgf") # use LaTeX fonts for plots, must occur before importing "plt"
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import json
import psutil # for memory usage
import signal
import colorsys
from matplotlib.colors import ListedColormap
from bs4 import BeautifulSoup
from surfplot import Plot
from brainspace.datasets import load_parcellation
from brainspace.mesh.mesh_io import read_surface
from neuromaps.datasets import fetch_fslr
from PIL import Image, ImageDraw, ImageFont
import nibabel as nib
from nilearn.plotting import plot_stat_map  # pip install nilearn
from nilearn import datasets
from nilearn import image
from nilearn import plotting






global_simulated = False        # simulate runs to re-use file names without running the whole pipeline
number_cores = 3                # used to speed-up antsRegistrationSyN.sh if possible (careful when multiple subjects are processed in parallel)



exec(open("./load_data.py").read()) # now, the dictionary 'datasets' is available.
short_names = list(datasets.keys())
#print ("available datasets: " + str(short_names))






try:
    task = sys.argv[1]
    subject = sys.argv[2]
except IndexError:
    print("usage: 'python3 main.py layerfmri 599', 'python3 main.py recon 574rescan', 'python3 main.py dc 618'")
    display_available_jobs()
    quit()


tasks = ["recon", "dc", "task", "recon-upsampled"]
if not task in tasks:
    print("### Error: Task not found '" + task + "'. Available tasks: "  + str(tasks))
    quit()

short_names.append("summary") # e.g., dc has a function "summary": "python3 main.py dc summary"
if not subject in short_names:
    print("### Error: Subject not found '" + subject + "'. Available subjects: "  + str(short_names))
    quit()
else:
    if subject != "summary":
        dataset = datasets[subject]
        label = dataset["label"] # e.g., "s599-7T"
    else:
        label = "summary_" + task



exec(open("./auxiliary/prepare_directories.py").read())     # initialize directories
exec(open("./auxiliary/auxiliary.py").read())               # load auxiliary methods (e.g., running commands and time keeping)
exec(open("./auxiliary/plotting_auxiliary.py").read())
exec(open("./auxiliary/file_management.py").read())         # upload results to XNAT
exec(open("./auxiliary/color.py").read())
exec(open("./auxiliary/plotting.py").read())
exec(open("./processing_fsl-feat.py").read())
exec(open("./processing.py").read())                        # main methods for fmri processing (e.g., preprocessing, coregistrations)
exec(open("./main_functions.py").read())                    # main functions: layerfmri, dc, recon-all
simulation_switch(0)










# perform one task selected from: recon, task, dc
print(subject + ": " + task)

if task == "dc" and subject == "summary": # "python3 main.py dc summary"
    dc_summary()
    quit()
    

if task == "inverse":
    inverseDCmatrix(subject)

if task == "recon-upsampled":
    layerfmri_reconall_after_upsampling(subject)

if task == "recon":
    anatomicalT1        = dataset["anatomicalT1"]
    reconAnatomical(subject, anatomicalT1)
    
    
if task == "task":
    anatomicalT1        = dataset["anatomicalT1"]
    fieldmap_mag        = dataset["fieldmap_mag"]
    fieldmap_phase      = dataset["fieldmap_phase"]
    functional          = dataset["functionalTap"]


if task in ["task"]:  # layerfmri can be spilt into a first part until recon-all and a second part
    anatomicalT1 = dataset["anatomicalT1"]
    fieldmap_mag = dataset["fieldmap_mag"]
    fieldmap_phase = dataset["fieldmap_phase"]


    activity = ""

    if dataset["functionalTap"] != "N/A":
        functional = dataset["functionalTap"]
        activity = "rft"

    if dataset["functionalChecker"] != "N/A":
        functional = dataset["functionalChecker"]
        activity = "vistim"

    if activity == "":
        print("No activity found.")
        quit()


    if task == "task":
        layerfmri_task(subject, anatomicalT1, fieldmap_mag, fieldmap_phase, functional, activity)


    
    
if task == "dc":
    DC_anatomical           = dataset["DC_anatomical"]
    DC_dcw                  = dataset["DC_dcw"]
    DC_mean_bold            = dataset["DC_mean_bold"]
    DC_mean_bold_GM_mask    = dataset["DC_mean_bold_GM_mask"]
    DC_GM_mask              = dataset["DC_GM_mask"]
    DC_WM_mask              = dataset["DC_WM_mask"]
    DC_CSF_mask             = dataset["DC_CSF_mask"]
    
    dc_layers(subject, DC_anatomical, DC_dcw, DC_mean_bold, DC_mean_bold_GM_mask, DC_GM_mask, DC_WM_mask, DC_CSF_mask)



        








    

    
    
    
    
    




