# auxiliary functions for running commands and measuring time



def simulation_switch(input): # simulate runs to re-use file names without running the whole pipeline
    # returns True or False
    # input: -1: do not change, 0: not simulated, 1: simulated

    simfile = global_subject_folder + "simulation.txt"

    if input == -1:
        return os.path.exists(simfile)

    if input == 0:
        if os.path.exists(simfile):
            os.remove(simfile)
        return False

    if input == 1:
        with open(simfile, "a+") as file:
            file.write("simulated")
        return True

    return False



def sec2HHMMSS (duration): # e.g., 72 seconds => 01:12
    hours = seconds // 3600
    minutes = (seconds % 3600) // 60
    seconds = seconds % 60

    hh = str(ohurs)
    if hours < 10:
        hh = "0" + hh

    mm = str(minutes)
    if minutes < 10:
        mm = "0" + mm

    ss = str(seconds)
    if seconds < 10:
        ss = "0" + ss

    time_format = hh + ":" + mm + ":" + ss  # f"{minutes:02d}:{seconds:02d}"
    return time_format


def sec2HHMMSS (duration): # e.g., 72 seconds => 00:01:12
    hours = duration // 3600
    minutes = (duration % 3600) // 60
    seconds = duration % 60

    hh = str(hours)
    if hours < 10:
        hh = "0" + hh

    mm = str(minutes)
    if minutes < 10:
        mm = "0" + mm

    ss = str(seconds)
    if seconds < 10:
        ss = "0" + ss

    time_format = hh + ":" + mm + ":" + ss  # f"{minutes:02d}:{seconds:02d}"
    return time_format




def run_log (text): # writes log to run.txt
    with open(global_runTXT, "a+") as file:  # global_subject_folder is defined in main.py
        file.write(text + "\n")

def run(cmd, no_sim=False):
    if no_sim or (not simulation_switch(-1)): # simulate runs to re-use file names without running the whole pipeline

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
                    print("\nrunning time-demanding operation: " + cmd)


        if not large: # otherwise already performed
            silent = False
            if "rm " in cmd or "mkdir" in cmd:
                result = os.popen(cmd + " > /dev/null 2>&1").read() # prevent the displaying of an output
                silent = True
            else:
                result = os.popen(cmd).read()

        end_time = time.time()
        duration = int(end_time - start_time)
        time_format = sec2HHMMSS (duration)

        if not silent:
            print("\n" + str(time_format) + "  " + cmd)

        run_log(time_format + " " + cmd)

        error_words = ["dumped", "Abort", "abnormal", "no such file"]
        for word in error_words:
            if word in result:
                print("\n\n# Detected an error, pipeline quits.")
                quit()


    else:
        print("\nsimulated: " + cmd)


def runTimer(cmd, title): # same as run, but for larger jobs to show a timer.
    # returns maximal RAM usage in GB

    start_time = time.time()
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
    print("Press enter to continue...")
    input()
    print("Script continued!")


def now():
    current_datetime = datetime.datetime.now()
    return str(current_datetime.strftime("%Y-%m-%d_%H.%M"))


def nice_round(number, direction="ceil"): # only keep the first two relevant digits and round up
    if number < 100:
        if number > 1:
            return int(number+1)
        else:
            return number
    else:
        number = int(number)

    significant = int(str(number)[:2])
    digits = len(str(number))

    if direction == "ceil":
        difference = 1
    else: # "floor"
        difference = -1

    return int(math.pow(10, digits-2)*(significant + difference))







def checkFiles (filelist, response=False): # checks if all files in list exist, if cales with second parameter respose, it does not quit and, instead gives an answer
    error = False
    for file in filelist:
        if not os.path.isfile(file):
            error = True
            print ("\n checkFiles - file not found: " + file)

    if error:
        if not response:
            if simulation_switch(-1) == 0: # only quit if this is not a simulated run (e.g., large files in the processing pipeline might get deleted to save sotrage space).
                quit()
            else:
                print ("Error, but will continue because this is a simulated run.")

        return False
    return True




        
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
                "LN2_PROFILE",
                "LN_INFO"]



    error = False
    for command in commands:
        testing = "command -v <cmd>"
        testing = testing.replace("<cmd>", command)
        probe = os.popen(testing).read()


        if len(probe) < 2:
            error = True
            print("### error: command not found: " + command)

    if error:
        print ("An error occured. Please make sure that the following are installed: FSL, ANTs, AFNI, freesurfer, laynii")
        quit() # an error occured: command not found, cannot run pipeline or part thereof


def copy_overview_files(subfolder, overview_files, no_replace=False): # copy main files to separate folders
    # example use:
    #       overview_files = []
    #       overview_files.append("output.nii.gz")
    #       copy_overview_files("DC_" + subject, overview_files)
    # no_replace: delete previous folder if true

    folder_base = "results/overview"
    folder_target = folder_base + "/" + subfolder

    jobs = []
    jobs.append("mkdir "        + folder_base)

    if not no_replace:
        jobs.append("rm -rf -f "    + folder_target)
    jobs.append("mkdir "        + folder_target)

    for file in overview_files:
        jobs.append("cp " + file + " " + folder_target + "/")

    for job in jobs:
        run(job, no_sim=True)

def fsl_info (file): # returns dict with info about dimensions
    # returns a dcit with image information, e.g., x_dim (number slices), x_pixdim (mm per slice), x_min (min absolute position in mm), x_max (max absolute position in mm)

    img = nib.load(file)

    # Get the header of the file
    header = img.header

    # definitions in nifti:                                     typical in mm
    #  x = dim1: left -> right         (same in DICOM)          -65 .. +65
    #  y = dim2: posterior -> anterior (opposite in DICOM)      -75 .. +80
    #  z = dim3: caudal -> cranial     (opposite in DICOM)      -60 .. +35


    # 3 dimensions (sic!)
    dim = list(header.get_data_shape())   #   number slices
    pixdim = list(header.get_zooms())     #   mm per slice
    qoffsets = list(img.affine[:3, 3])    #   offsets in mm (! typically first is positive, and second/third are negative)

    qoffsets_pos = [abs(x) for x in qoffsets]


    # get image origin

    origin = list(header.get_sform()[:3, 3])
    #print("origin: " + str(origin))





    #print(str(dim))
    #print(str(pixdim))
    #print(str(qoffsets))

    info = {}
    xyz=["x", "y", "z"]
    for i in range (0, 2+1):
        info[xyz[i]+"_dim"]         = int(dim[i])                                                                   # number slices
        info[xyz[i]+"_pixdim"]      = round(float(pixdim[i]), 2)                                                    # mm per slice
        info[xyz[i]+"_min"]         = round(float(0.0-qoffsets_pos[i]),2)                                           # min abs position in mm
        info[xyz[i] + "_max"]       = round(float(0.0-qoffsets_pos[i]) + (float(dim[i]) * float(pixdim[i])), 2)     # max abs position in mm
        info[xyz[i] + "_origin"]    = round(float(origin[i]), 2)

        # print(str(info))
    return info



def color_within_margin(color1, color2, margin_percent):
    threshold = 255 * margin_percent / 100
    return all(abs(c1 - c2) <= threshold for c1, c2 in zip(color1, color2))


def change_color(myimg, target_color, new_color):
    width, height = myimg.size
    pixels = myimg.load()

    for x in range(width):
        for y in range(height):
            current_color = pixels[x, y]

            if current_color == target_color:
                pixels[x, y] = new_color

    return myimg




def nifti_edit_dimensions(original_image, new_image, x_dim=None, y_dim=None, z_dim=None):
    # changes the dimension values in the header file (no direct methon possible, thence, this function)
    # can be used to extend or crop the image space
    # Slices can be added or removed only at the positive end, but keeps the origin of the image and, therefore, the validity of the header file.

    # definitions in nifti:                                     typical in mm
    #  x = dim1: left -> right         (same in DICOM)          -65 .. +65
    #  y = dim2: posterior -> anterior (opposite in DICOM)      -75 .. +80
    #  z = dim3: caudal -> cranial     (opposite in DICOM)      -60 .. +35

    header_file = temp_path + "header.nii.gz"
    xml_file = temp_path + "header.xml"
    xml_text = os.popen("fslhd -x " + original_image).read()
    # example:
    #   <nifti_image
    #   image_offset = '2880'
    #   ndim = '3'
    #   nx = '312'
    #   ny = '376'
    #   nz = '211'
    #   nt = '1'
    #   dx = '0.5'
    #   dy = '0.5'
    #   dz = '0.5'
    #   dt = '0'
    #   (...)



    if x_dim != None:
        replacement = "nx = '" + str(x_dim) + "'"
        xml_text = re.sub(r"nx = '[0-9]+'", replacement, xml_text)

    if y_dim != None:
        replacement = "ny = '" + str(y_dim) + "'"
        xml_text = re.sub(r"ny = '[0-9]+'", replacement, xml_text)

    if z_dim != None:
        replacement = "nz = '" + str(z_dim) + "'"
        xml_text = re.sub(r"nz = '[0-9]+'", replacement, xml_text)

    with open(xml_file, "w") as file:
        file.write(xml_text)

    jobs = []
    jobs.append ("fslcreatehd <xml_file> <header_file>") # create a .nii.gz file that contains only the header but no image data
    jobs.append("rm <new_image>")
    jobs.append("3dresample -master <header_file> -rmode Li -prefix <new_image> -input <original_image>") # merge header file and image data into an new file

    for job in jobs:
        job = job.replace("<original_image>", original_image)
        job = job.replace("<new_image>", new_image)
        job = job.replace("<header_file>", header_file)
        job = job.replace("<header_file>", header_file)
        job = job.replace("<xml_file>", xml_file)

        run (job)

    return new_image


def nifti_translate_image (image, new_image, x_mm=0, y_mm=0, z_mm=0):
    # translates image in a given direction, preserves the original space (i.e., use "nifti_extend_header" to make sure nothing will be cut off at the image boarder)

    # definitions in nifti:                                     typical in mm
    #  x = dim1: left -> right         (same in DICOM)          -65 .. +65
    #  y = dim2: posterior -> anterior (opposite in DICOM)      -75 .. +80
    #  z = dim3: caudal -> cranial     (opposite in DICOM)      -60 .. +35


    matrix = "1 0 0 <x_mm>\n0 1 0 <y_mm>\n0 0 1 <z_mm>\n0 0 0 1"
    matrix = matrix.replace("<x_mm>", str(x_mm)) # positive direction: to the right
    matrix = matrix.replace("<y_mm>", str(y_mm)) # positive direction: to anterior
    matrix = matrix.replace("<z_mm>", str(z_mm)) # positive direction: to cranial



    matrix_file = temp_path + "shift_matrix.mat"
    with open(matrix_file, "w") as file:
        file.write(matrix)

    job = "flirt -in <image> -ref <image> -applyxfm -init <matrix> -out <output>"

    job = job.replace("<image>", image)
    job = job.replace("<matrix>", matrix_file)
    job = job.replace("<output>", new_image)
    run (job)








def activity_thresholding_for_plotting (activity_map): # Provides min and max values for plotting the activity such that it looks acceptable and only relevant areas are shown. Must be robust for different tasks.
    # returns thresholds: [min_intensity, max_intensity]
    # max intensity: 1mm resampling, bin width 30, at least 5 in bin, ceil to next 100
    # min intensity: see below.

    #   e.g.,
    #           574rescan (task-rft):   [350.0, 1500]
    #           640 (task-vistim):      [250.0, 1100]
    #           641 (task-vistim):      [200.0, 1100]

    threshold           = 50
    bin_size            = 30
    min_content_upper   = 5     # determins max intensity
    downsample_mm       = 1

    [image_small] = resampling(downsample_mm, [activity_map], method="NN", folder=temp_path)  # speed up
    img = nib.load(image_small)
    data = img.get_fdata()

    # Flatten the data array and calculate the intensity distribution
    intensity_values, intensity_counts = np.unique(data.flatten(), return_counts=True)

 

    # Create bins of width 10 for intensity values above the threshold
    bins = np.arange(threshold, np.max(intensity_values) + 10, 10)
    binned_counts, _ = np.histogram(intensity_values, bins=bins, weights=intensity_counts)

    max_bin = bins[np.where(binned_counts >= min_content_upper)[0][-1]] # value of highest relevant bin (contains at least 5 samples)
    max_intensity = (int(max_bin) // 100) * 100 # floor down to next 100



    # min_intensity is calculated by taking the total number of voxels above half of max_intensity, then lowering the threshold by it by 50, until at least 4.5x of the first total number is reached.
    # aim: include all relevant voxels but not too many (cannot determine a size relative to the whole image because area varies with task)
    num_voxels_above_half_max = np.sum(data.flatten() > max_intensity / 2)
    min_intensity = max_intensity / 2
    while (np.sum(data.flatten() > min_intensity) < 4.5 * num_voxels_above_half_max):
        min_intensity -= 50


    result = [min_intensity, max_intensity]
    #print(str(result))
    return result


def cluster_filtering (input_image, min_intensity=None, min_rel_clustersize=None, output_image=None):
    # filters, typically an activity map, such that all clusters that do not contain values of a given minimum intensity are removed
    # min_intensity will be determiend automatically if not given
    # min_rel_clustersize is optional. it refers to the minimum size relative to the largest cluster.
    #   i.e., if given, at least one of the two criteria needs to be fulfilled: 1. minimum intensity of a cluster's maximum; OR 2. minimum size of a cluster relative to the largest cluster

    if output_image == None:
        output_image = input_image.replace(".nii.gz", "_cluster-filtered.nii.gz")

    if min_intensity == None: # automatic mode
        [thr_min, thr_max] = activity_thresholding_for_plotting(input_image) # determins threshold that are useful for plotting (i.e., not the total min/max intensity, but a recommendation for setting the color bar)

        min_intensity = 0.2*thr_min + 0.8*thr_max # only permit clusters with a relatively high intensity

    jobs = []
    #jobs.append("cluster --in=<input_image> --thresh=<min_intensity> --mm --osize=output_cluster_size.txt --oindex=output_cluster_index.txt")
    #jobs.append("3dclust -1Dformat -1dindex 1 -1tindex 1 -savemask cluster_mask.nii.gz -1thresh 3.1 26 <input_image>")
    #jobs.append("3dClusterize -inset <input_image> -ithr 1 -NN 2 -clust_nvox 10 -bisided -pref_map cluster_map -pref_dat cluster_data")
    #jobs.append("fslmaths <input_image> -mas output_cluster_index.nii.gz <output_image>")


    #jobs.append("3dclust -1clip 0.3  5  3000 <input_image>")
    #jobs.append("3dclust -1thresh 0.1 -NN2 -savemask cluster_mask.nii.gz <input_image>  >> clust.txt")
    jobs.append ("rm cluster_mask.nii.gz")
    jobs.append("rm clust.txt")
    #jobs.append("3dclust -1clip 200 5 0 <input_image> -savemask cluster_mask.nii.gz >> clust.txt")
    jobs.append("3dClusterize -inset <input_image> -ithr 0 -NN 1 -clust_nvox 10 -1sided RIGHT_TAIL 3.313 -clust_nvox 157 -out_mask cluster_mask.nii.gz -binary >> clust.txt")

    for job in jobs:
        job = job.replace("<input_image>", input_image)
        job = job.replace("<min_intensity>", str(min_intensity))
        job = job.replace("<output_image>", output_image)
        run (job)

    return output_image



def aligned_sort(list_unsorted, list_criteria): # sorts first list by ascending order of a second list
    zipped_lists = zip(list_criteria, list_unsorted)
    sorted_lists = sorted(zipped_lists)
    return [element for _, element in sorted_lists]



def calculate_average_error(float_list):
    sum_values = sum(float_list)
    avg_value = sum_values / len(float_list)

    sum_squares_diff = sum((val - avg_value) ** 2 for val in float_list)
    variance = sum_squares_diff / (len(float_list) - 1)
    error = math.sqrt(variance) # std_deviation


    return avg_value, error
