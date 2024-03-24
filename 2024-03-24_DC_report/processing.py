

def sample_and_delete (functional, nodelete=False): # avoid storting large intermediate files: create a sample and delete original file
    folder = output_path + "intermediary_samples"
    jobs = []
    jobs.append ("mkdir " + folder)
    folder += "/"



    sample_name = naming(functional, "sample", folder)
    jobs.append("fslroi <functional> <sample_name> 0 5")

    # jobs.append("rm <sample_name>")
    jobs.append("mkdir results/to_delete")
    jobs.append("mv <sample_name> results/to_delete/")


    if not nodelete:                        # i.e., create a sample, without deleting
        jobs.append("rm <functional>")

    for job in jobs:
        job = job.replace("<sample_name>", sample_name)
        job = job.replace("<functional>", functional)
        run (job)


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








def preprocessing_estimate_fieldmap (fieldmap_mag, fieldmap_phase):
    # estimates a fieldmap from the magnitude and phase images

    checkFiles([fieldmap_mag, fieldmap_phase])

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

    checkFiles([fieldmap_out])
    return fieldmap_out





def preprocessing_applyfieldmap(input_img, fieldmap_distortion, output_corrected):
    # uses fugue(FMRIB's Utility for Geometrically Unwarping EPIs) to perform an unwarping of an EPI image based on fieldmap data.
    checkFiles([input_img, fieldmap_distortion])
    dwell_time = get_dwell_time(input_img) # taken from json file (alternatively, DICOM header)

    task = "fugue -i <input_epi> --dwell=<dwell_time> --unwarpdir=y --loadfmap=<fieldmap> -u <result>" # "-savematrix <matrix_distortionCorr>" (--unwarpdir: "Use x, y, z, x-, y- or z- only.")
    task = task.replace("<input_epi>",                  input_img)
    task = task.replace("<dwell_time>",                 str(dwell_time))
    task = task.replace("<fieldmap>",                   fieldmap_distortion)
    task = task.replace("<result>",                     output_corrected)
    task += " --verbose"


    run(task)
    checkFiles([output_corrected])
    return output_corrected


def preprocessing_motioncorrection (input):
    checkFiles([input])
        
    output = naming (input, "motioncorr", output_path)
    
    task = "3dvolreg -twopass -prefix <epi_motionCorrected> <epi_input>"
    task = task.replace("<epi_input>", input)
    task = task.replace("<epi_motionCorrected>", output)

    run (task)
    checkFiles([output])
    sample_and_delete(input) # delete large intermediate files and keep only a sample
    sample_and_delete(output, nodelete=True)
        
    return output



def preprocessing_biasfield4D (input):
    # input is a functional 4D file
    # estimates bias field based on first entry and removes same bias field from all other elements in the time series
    # assumtion: For small distances, the bias field does not change much, in particular, if the functional images are already motion corrected.
    # Please Note: The bias field is not estimated for every element of the time series individually because this is about fMRI data.

    checkFiles([input])


    output = naming (input, "biasFieldRemoval", output_path)
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

    checkFiles([output])

    sample_and_delete(input)  # delete large intermediate files and keep only a sample
    return output













def preprocessing_undistortion(input, fieldmap_mag, fieldmap_phase, output):
    #   input: functional image, fieldmap_mag, fieldmap_phase, output_filename
    #   returns list of output filenames (i.e., the undistorted functional 4D image file and the estimated bias field)
    #   output: undistorted functional image, biasfield

    checkFiles([input, fieldmap_mag, fieldmap_phase])


    # 1. fieldmap estimation
    fieldmap_distortion = preprocessing_estimate_fieldmap(fieldmap_mag, fieldmap_phase)  # uses fsl_prepare_fieldmap (and bet for skullstripping, fslmaths -ero for eroding, fsl_anat --strongbias for bias correction, epi-reg for registration, fslmaths -div/sub/mul for radian conversion)
   
    # 2. apply the fieldmap to undistort the input image
    output = preprocessing_applyfieldmap(input, fieldmap_distortion, output) # uses fugue (set output file name to shorten file pahts)

    checkFiles([output, fieldmap_distortion])
    return [output, fieldmap_distortion]





def coregistrations(functional_reference, anatomical, MNI, atlas, functional_first):
    # typical use: [anatomical_bet_functionalspace, anatomical_bet_functionalaligned, atlas_functionalspace, matrix_anatomical_to_functional] = coregistrations (functional_reference, anatomicalT1, constant_MNItemplate, constant_atlas)

    checkFiles([functional_reference, anatomical, constant_MNItemplate, constant_atlas])

    anatomical_bet                      = temp_path + "anatomical_bet.nii.gz"
    anatomical_bet_functionalspace      = output_path + "anatomical_bet_functionalspace.nii.gz"
    anatomical_bet_functionalaligned    = output_path + "anatomical_bet_functionalaligned.nii.gz"
    MNI_anatomicalspace                 = temp_path + "MNI_anatomicalspace.nii.gz"                    # temp_path + "MNI_anatomicalspace.nii.gz" # gets too long otherwise
    MNI_functionalspace                 = naming(MNI, "functionalspace", output_path)                 # output_path + "MNI_functionalspace.nii.gz"
    atlas_functionalspace               = naming(atlas, "functionalspace", output_path)               # output_path + "atlas_functionalspace.nii.gz"
    functional_reference_fastbet        = naming(functional_reference, "fastbet", output_path)
    matrix_anatomical_to_functional     = output_path + "matrix_anatomical_to_functional.mat"
    matrix_MNI_to_functional            = output_path + "matrix_MNI_to_functional.mat"

    if False and checkFiles([anatomical_bet_functionalspace, anatomical_bet_functionalaligned, atlas_functionalspace, matrix_anatomical_to_functional], True): # skip time-demanding operations if possible
        print("\n\nSkipped coregistration step, all files are already available.\n\n")
        return [anatomical_bet_functionalspace, anatomical_bet_functionalaligned, atlas_functionalspace, matrix_anatomical_to_functional]

    jobs = []

    # Step 1: Register anatomical to functional reference with rigid transformation
    jobs.append("mri_synthstrip -i <functional_reference> --no-csf -o <functional_reference_fastbet>") # only used to match the anatomical image. proper bet needs to performed later with masks.
    jobs.append("mri_synthstrip -i <anatomical> --no-csf -o <anatomical_bet>")

    if False:
        jobs.append("flirt -in <anatomical_bet> -ref <functional_reference_fastbet> -noresample -omat <matrix_anatomical_to_functional> -out <anatomical_bet_functionalspace>")
        jobs.append("flirt -in <anatomical_bet> -ref <functional_reference_fastbet> -noresample -omat <matrix_anatomical_to_functional> -out <anatomical_bet_functionalaligned>")
        # works in most cases, but, for example, for 549rescan the upper few centimeters are missing
        # mind that image dimensions stay the same during transforms, this is not the error source
        # Therefore, two different versions need to be generated, <anatomical_bet_functionalspace> and <anatomical_bet_functionalaligned>. See explanation in main pipeline when "coregistrations" is called.
    else:
        # more robust version: align and then apply matrix (only in second step, the upper part will be preserved)

        if True: #not checkFiles([matrix_anatomical_to_functional], response=True):
            jobs.append("bash antsRegistrationSyN.sh -d 3 -f <functional_reference_fastbet> -m <anatomical_bet> -n <number_cores> -t r -o <temp/>ANTsOutput")
            jobs.append("rm <matrix_anatomical_to_functional>")
            jobs.append("mv <temp/>ANTsOutput0GenericAffine.mat <matrix_anatomical_to_functional>")
        else:
            jobs.append("echo skipped: coregistrations anatomical to functional, pre-computed result found: " + matrix_anatomical_to_functional)

        jobs.append("antsApplyTransforms -d 3 -i <anatomical_bet> -r <functional_first> -o <anatomical_bet_functionalspace> -t <matrix_anatomical_to_functional>")
        jobs.append("antsApplyTransforms -d 3 -i <anatomical_bet> -r <anatomical> -o <anatomical_bet_functionalaligned> -t <matrix_anatomical_to_functional>")






        # Step 2: Register MNI.nii.gz to anatomical_bet_functionalspace with a nonlinear transformation and get matrix
    if True: # not checkFiles([matrix_MNI_to_functional], response=True):
        jobs.append("bash antsRegistrationSyN.sh -d 3 -f <anatomical_bet_functionalaligned> -m <MNI> -n <number_cores> -t s -o <temp/>ANTsOutput") # deformable transformation # mostly interested in the transformation matrix, therefore, not using the uplampled image.
        # parameters:
        #               f: fixed image
        #               t: 'r' (rigid), 'a': rigid + affine (2 stages), 's': rigid + affine + deformable syn (3 stages)
        #               d: dimension (2 or 3, not 4)
        #               m: moving image (that needs transform)
        #               o: outputfile
        # output:
        #       [outputname]Affine.mat
        #       [outputname]Warped.nii.gz

        jobs.append("rm <MNI_functionalspace>")
        jobs.append("rm <matrix_MNI_to_functional>")
        jobs.append("mv <temp/>ANTsOutputWarped.nii.gz <MNI_functionalspace>")
        jobs.append("mv <temp/>ANTsOutput0GenericAffine.mat <matrix_MNI_to_functional>")
    else:
        jobs.append("echo skipped: coregistrations MNI to functional, pre-computed result found: " + matrix_MNI_to_functional)






    # Step 4: Apply matrix to atlas
    jobs.append("antsApplyTransforms -d 3 -i <atlas> -r <anatomical_bet_functionalaligned> -o <atlas_functionalspace> -n NearestNeighbor -t <matrix_MNI_to_functional>")
    # parameters:
    #               -t transformation matrix (can use multiple, each preceded by a -t)
    #               -i input
    #               -o output
    #               -r reference image: "For warping input images, the reference image defines the spacing, origin, size, and direction of the output warped image."
    #





    for job in jobs:
        job = job.replace("<anatomical>",                       anatomical)
        job = job.replace("<anatomical_bet>",                   anatomical_bet)
        job = job.replace("<functional_reference>",             functional_reference)
        job = job.replace("<functional_reference_fastbet>",     functional_reference_fastbet)
        job = job.replace("<matrix_anatomical_to_functional>",  matrix_anatomical_to_functional)
        job = job.replace("<matrix_MNI_to_functional>",         matrix_MNI_to_functional)
        job = job.replace("<anatomical_bet_functionalspace>",   anatomical_bet_functionalspace)
        job = job.replace("<anatomical_bet_functionalaligned>", anatomical_bet_functionalaligned)
        job = job.replace("<MNI_functionalspace>",              MNI_functionalspace)
        job = job.replace("<atlas>",                            atlas)
        job = job.replace("<MNI>",                              atlas)
        job = job.replace("<atlas_functionalspace>",            atlas_functionalspace)
        job = job.replace("<temp/>",                            temp_path)
        job = job.replace("<number_cores>",                     str(number_cores))  # defined in main
        job = job.replace("<functional_first>",                 functional_first)

        run(job)



    result = [anatomical_bet_functionalspace, anatomical_bet_functionalaligned, atlas_functionalspace, matrix_anatomical_to_functional]
    checkFiles(result)
    return result







def perform_bet_with_masks (functional, anatomical_bet_functionalspace):
    # typical use: functional_bet = perform_bet_with_masks (functional_upsampled, anatomical_bet_functionalspace_upsampled)
    # !! both input images must be in the same same and have the same voxel size


    checkFiles([functional, anatomical_bet_functionalspace])

    anatomical_bet_functionalspace_mask = naming(anatomical_bet_functionalspace, "mask", temp_path)
    functional_bet = naming(functional, "bet", output_path)


    jobs = []
    jobs.append("fslmaths <anatomical_bet_functionalspace> -abs -bin <anatomical_bet_functionalspace_mask>")
    jobs.append("fslmaths <functional> -mas <anatomical_bet_functionalspace_mask> <functional_bet>")


    for job in jobs:
        job = job.replace("<functional>", functional)
        job = job.replace("<functional_bet>", functional_bet)
        job = job.replace("<anatomical_bet_functionalspace>", anatomical_bet_functionalspace)
        job = job.replace("<anatomical_bet_functionalspace_mask>", anatomical_bet_functionalspace_mask)

        run(job)

    checkFiles([functional_bet])
    return functional_bet




def resampling (voxelsize_mm, list_of_images, method=None, folder=None):
    # returns a list of filenames of the resampled images
    # typical use: [anatomical_bet_functionalspace_upsampled, functional_preprocessed_upsampled, maskWM_upsampled, maskGM_upsampled, maskCSF_upsampled] = upsampling (0.5, [anatomical_bet_functionalspace, functional_preprocessed, maskWM, maskGM, maskCSF])
    # method = "NN", "Li", None. If no sampling method is given, linear interpolation will be performed unless the filename suggests it is a mask or an atlas

    # about AFNI 3dresample: https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dresample.html
    #   -rmode RESAM     : use this resampling method
    #           e.g.  -rmode Linear
    #           default is NN (nearest neighbor)
    #
    #           The resampling method string RESAM should come
    #           from the set {'NN', 'Li', 'Cu', 'Bk'}.  These
    #           are for 'Nearest Neighbor', 'Linear', 'Cubic'
    #           and 'Blocky' interpolation, respectively.
    #   => use without quotation marks, e.g.,  3dresample -orient asl -rmode NN -prefix asl.dset -input in+orig


    checkFiles(list_of_images)

    isotropic = "{:.2f}".format(voxelsize_mm)

    list_upsampled = []

    for image in list_of_images:
        jobs = []
        jobs.append ("rm <image_up>")  # 3dresample does not replace existing files and states "** ERROR: dataset NOT written to disk! failure: cannot write dataset, exiting..."
        jobs.append("3dresample -dxyz <mm> <mm> <mm> -rmode <sampling_method> -prefix <image_up> -inset <image>")

        image_up = image.replace(".nii.gz", "_resampled_" + isotropic + "mm.nii.gz")
        if folder != None:
            image_up = folder + os.path.basename(image_up)

        if ("data/" in image_up) or ("atlas/" in image_up): # resampling an original file - do not change content of folder data/
            image_up = output_path + os.path.basename(image_up)


        if method != None:
            sampling_method = method
        else:
            if ("mask" in image) or ("atlas" in image) or ("glasser" in image) or ("schaefer" in image) or ("parcel" in image) or ("WM" in image) or ("GM" in image) or ("CSF" in image): # use nearest neighbor
                sampling_method = "NN" # nearest neighbor
            else:
                sampling_method = "Li" # linear

        for job in jobs:
            job = job.replace("<image>", image)
            job = job.replace("<image_up>", image_up)
            job = job.replace("<mm>", isotropic)
            job = job.replace("<sampling_method>", sampling_method)
            run (job)

        list_upsampled.append(image_up)

    checkFiles(list_upsampled)
    return list_upsampled






def create_LN2_layers (title, LN2_input_labels, numberLayers, anatomical=None):
    # typical use: layer_results = createLayers (dict_label_filenames, numberLayers)
    #              selected_layers = layer_results["layers_equivol"]
        # LN2_input_labels is the label file with 3 intensities for GM, WM, and CSF.
        # anatomical is an optional parameter to needed if the label file (3 layers), and the output layers file should be shown; anatomical_functionalspace
        # returns list of six output files ([layers, metrix, layerbins] x [equidist, equivol])

    # title is only ofr file management e.g., "DC". Not a plot title
    if " " in title:
        print("Error: create_LN2_layers was called with a title that was only meant for file management. It contains unallowed characters (e.g., spaces).")
        quit() # Prevent resulting script from failing.


    checkFiles([LN2_input_labels])

    output_files = {} # dict of six output files, e.g. "layers_equivol"  # ([layers, metrix, layerbins] x [equidist, equivol])
    jobs = []

    LN2_results_folder = "<output/>LN2_results"

    jobs.append("mkdir <LN2_results_folder>")
    jobs.append("rm <temp/>LN2temp.nii.gz")
    jobs.append ("cp <LN2_input_labels> <temp/>LN2temp.nii.gz")                                            # needed because LN2_LAYERS has many output files that should not all be in the main output folder
    jobs.append("LN2_LAYERS -rim <temp/>LN2temp.nii.gz -nr_layers <nr_layers> -equal_counts -equivol")      # equivol creates equivol layers in addition to equidist. resulting file is, for example, "[LN2temp]_metric_equidist.nii.gz"

    for type1 in ["layers", "metric", "layerbins"]:
        for type2 in ["equidist", "equivol"]:  # layerbins needed for profile

            output_file = "<LN2_results_folder>/" + title + "_result_" + str(numberLayers) + "-layers_" + "<type1>_<type2>.nii.gz"          # output: e.g., output/layer_results/layniiLabels_layerbins_equivol.nii.gz
            output_file = output_file.replace("<type1>", type1)
            output_file = output_file.replace("<type2>", type2)
            output_file = output_file.replace("<layerresults/>", layniiresult_path)
            tag = type1 + "_" + type2

            output_file = output_file.replace("<LN2_results_folder>", LN2_results_folder)
            output_file = output_file.replace("<output/>", output_path)
            output_files[tag] = output_file


            job = "rm <output_file>"
            job = job.replace("<output_file>", output_file)
            jobs.append(job)

            job = "cp <temp/>LN2temp_<type1>_<type2>.nii.gz <output_file>"
            job = job.replace("<output_file>", output_file)
            job = job.replace("<type1>", type1)
            job = job.replace("<type2>", type2)
            jobs.append(job)



    for job in jobs:
        job = job.replace("<LN2_results_folder>", LN2_results_folder)
        job = job.replace("<temp/>", temp_path)
        job = job.replace("<output/>", output_path)
        job = job.replace("<nr_layers>", str(numberLayers))  # e.g., 9
        job = job.replace("<LN2_input_labels>", LN2_input_labels)
        run (job)



    if anatomical != None and not "temp" in title:
        plot_type = "3layers"




        plot_overlay(anatomical, LN2_input_labels, LN2_input_labels.replace(".nii.gz", "_ortho.png"),       title=None, type="3layers", view="ortho", upsampling = 0.4)


        for selection in ["layerbins_equivol"]: #["layerbins_equivol", "layerbins_equidist"]:
            selected_layers = output_files[selection]
            plot_type = str(numberLayers) + "layers"
            plot_title = "LN2_layers output: " + selection
            plot_overlay(anatomical, selected_layers, selected_layers.replace(".nii.gz", "_ortho.png"),     title=None, type=plot_type, view="ortho", upsampling = 0.4)


    #output_files_list =  [entry for entry in output_files]
    #checkFiles(output_files_list)

    return output_files





def create_functional_reference (functional, functional_preprocessed):
    # create functional reference image (intention: Functional image includes a time series of 300 samples. create average of the first 10 to obtain one image and remove noise. Input must be motion corrected.)
    # typical use:
    #       upsampled_isotropic_voxelsize = 0.3
    #       [functional_reference, functional_first] = create_functional_reference (functional, functional_preprocessed)


    checkFiles ([functional, functional_preprocessed])
    functional_reference = output_path + "functional_reference.nii.gz"
    first_original = output_path + "functional_original_first.nii.gz"

    jobs = []
    jobs.append("fslroi <functional_preprocessed> <temp/>functional_first_10_volumes.nii.gz 0 10")
    jobs.append("fslmaths <temp/>functional_first_10_volumes.nii.gz -Tmean <functional_reference>")
    jobs.append("fslroi <functional> <first_original> 0 1")

    for job in jobs:

        job = job.replace("<functional_reference>", functional_reference)
        job = job.replace("<functional>", functional)
        job = job.replace("<functional_preprocessed>", functional_preprocessed)
        job = job.replace("<first_original>", first_original)
        job = job.replace("<temp/>", temp_path)

        run(job)

    checkFiles([functional_reference])
    return [functional_reference, first_original]






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



def create_LN2_labels(voxelsize_mm, title, maskWM, maskGM, maskCSF): # creates label file (intensities 1, 2, 3) as an input for LN2layers
    # typical use: LN2_input_labels = create_LN2_labels (0.3mm, "main", maskWM, maskGM, maskCSF)
    # performes the following tasks:
    #   - upsampling to target sirotropic voxel size (freesurfer output is always 1mm)
    #   - smoothing where possible  (i.e., if CSF and WM are sufficiently spaced from each other, smoothing is applied. If not possible, keep the original NN upsampling)
    #   - creation of a 1-voxel rim around the GM volume tha tis either pial or superficial white matter
    #   - creation of a label file as input for LN2_layers (i.e., 1: pial, 2: superfixial WM, 3: GM)

    isotropic = "{:.2f}".format(voxelsize_mm)
    LN2_label_file = output_path + title + "_LN2_labels.nii.gz" #  in functional space
    partially_smoothed = output_path + "labels_partially_smoothed_" + isotropic + "mm.nii.gz"   # label file with GM and full WM and CSF (not jsut rim), smoothed and upsampled where applicable
    # intensities
    #   1: CSF / pial
    #   2: WM / superficial WM
    #   3: GM


    jobs = []

    jobs.append("rm <temp/>labels_joined_upsampled_NN.nii.gz")
    jobs.append("rm <temp/>labels_joined_upsampled_Li.nii.gz")
    jobs.append("rm <temp/>LN2_labels_layer1_1mm.nii.gz")
    jobs.append("rm <temp/>LN2_labels_layer3_1mm.nii.gz")
    jobs.append("rm <temp/>LN2_labels_layer3.nii.gz")
    jobs.append("rm <temp/>labels_joined_1mm.nii.gz")
    jobs.append("rm <temp/>labels_smoothed_CSF.nii.gz")
    jobs.append("rm <temp/>labels_smoothed_CSF_mask.nii.gz")
    jobs.append("rm <temp/>labels_smoothed_WM.nii.gz")
    jobs.append("rm <temp/>labels_smoothed_WM_mask.nii.gz")
    jobs.append("rm <temp/>labels_smoothed_GM.nii.gz ")
    jobs.append("rm <temp/>labels_smoothed_GM_mask.nii.gz")
    jobs.append("rm <temp/>labels_smoothed_CSF_mask_1.nii.gz")
    jobs.append("rm <temp/>labels_smoothed_WM_mask_2.nii.gz")
    jobs.append("rm <temp/>labels_smoothed_GM_mask_3.nii.gz")
    jobs.append("rm <temp/>labels_smoothed_joined.nii.gz")
    jobs.append("rm <temp/>labels_CSF_selection.nii.gz")
    jobs.append("rm <temp/>labels_WM_selection.nii.gz")
    jobs.append("rm <temp/>labels_CSF_selection_smear.nii.gz")
    jobs.append("rm <temp/>labels_WM_selection_smear.nii.gz")
    jobs.append("rm <temp/>labels_no_smoothing_mask.nii.gz")
    jobs.append("rm <temp/>labels_smoothing_mask.nii.gz")
    jobs.append("rm <temp/>labels_smoothed_part.nii.gz ")
    jobs.append("rm <temp/>labels_not_smoothed_part.nii.gz")
    jobs.append("rm <partially_smoothed>")


    jobs.append("rm <temp/>LN2_labels_layer1.nii.gz")
    jobs.append("rm <temp/>LN2_labels_layer2.nii.gz")
    jobs.append("rm <temp/>LN2_labels_layer3.nii.gz")
    jobs.append("rm <LN2_labels>")
    jobs.append("fslmaths <maskCSF> -mul 1 <temp/>LN2_labels_layer1_1mm.nii.gz") # this is only for smoothing and different from the final label file. GM and WM intensities are swapped here
    jobs.append("fslmaths <maskWM> -mul 3 <temp/>LN2_labels_layer3_1mm.nii.gz")
    jobs.append("fslmaths <maskGM> -mul 2 <temp/>LN2_labels_layer2_1mm.nii.gz")
    jobs.append("fslmaths <temp/>LN2_labels_layer1_1mm.nii.gz -add <temp/>LN2_labels_layer2_1mm.nii.gz -add <temp/>LN2_labels_layer3_1mm.nii.gz <temp/>labels_joined_1mm.nii.gz")

    jobs.append("3dresample -dxyz <mm> <mm> <mm> -rmode NN -prefix <temp/>labels_joined_upsampled_NN.nii.gz -inset <temp/>labels_joined_1mm.nii.gz") # nearest neighbor (where smoothing is not possible)
    if False:
        # not working, this introduces additionas WM lines between GM and CSF
        jobs.append("rm <temp/>labels_joined_upsampled_Li.nii.gz")
        jobs.append("3dresample -dxyz <mm> <mm> <mm> -rmode Li -prefix <temp/>labels_joined_upsampled_Li.nii.gz -inset <temp/>labels_joined_1mm.nii.gz") # linear (where smoothing is possible)
    else:
        # alternative: smoothing with spherical kernel and mean values
        jobs.append("fslmaths <temp/>labels_joined_upsampled_NN.nii.gz -kernel sphere 1.0 -fmean <temp/>labels_joined_upsampled_Li.nii.gz -odt double") # use sphere-shaped kernel with radius 1mm (corresponds to recon-all output) for smearing




    # find out where WM and CSF have enough distance, such that smoothing is applicable

    jobs.append("fslmaths <temp/>labels_joined_upsampled_Li.nii.gz -thr 0.5 -uthr 1.4999 <temp/>labels_smoothed_CSF.nii.gz -odt double") # only temporary numbers: WM = 3, GM=2, CSF=1, otherwise different in final label file.
    jobs.append("fslmaths <temp/>labels_smoothed_CSF.nii.gz -abs -bin <temp/>labels_smoothed_CSF_mask.nii.gz -odt int")

    jobs.append("fslmaths <temp/>labels_joined_upsampled_Li.nii.gz -thr 2.5 -uthr 3.4999 <temp/>labels_smoothed_WM.nii.gz -odt double")
    jobs.append("fslmaths <temp/>labels_smoothed_WM.nii.gz -abs -bin <temp/>labels_smoothed_WM_mask.nii.gz -odt int")

    jobs.append("fslmaths <temp/>labels_joined_upsampled_Li.nii.gz -thr 1.5 -uthr 2.4999 <temp/>labels_smoothed_GM.nii.gz -odt double")
    jobs.append("fslmaths <temp/>labels_smoothed_GM.nii.gz -abs -bin <temp/>labels_smoothed_GM_mask.nii.gz -odt int")

    jobs.append("fslmaths <temp/>labels_smoothed_CSF_mask.nii.gz -mul 1 <temp/>labels_smoothed_CSF_mask_1.nii.gz -odt int") # numbers now correspond to LN2_layers inout
    jobs.append("fslmaths <temp/>labels_smoothed_WM_mask.nii.gz -mul 2 <temp/>labels_smoothed_WM_mask_2.nii.gz -odt int")
    jobs.append("fslmaths <temp/>labels_smoothed_GM_mask.nii.gz -mul 3 <temp/>labels_smoothed_GM_mask_3.nii.gz -odt int")
    jobs.append("fslmaths <temp/>labels_smoothed_CSF_mask_1.nii.gz -add <temp/>labels_smoothed_WM_mask_2.nii.gz -add <temp/>labels_smoothed_GM_mask_3.nii.gz <temp/>labels_smoothed_joined.nii.gz")  # but not applicable everywhere

    if False: # not neccasary, smoothing works well without corrections
        # determine how to mix smoothed and non-smoothed labels
        jobs.append("fslmaths <temp/>labels_joined_upsampled_NN.nii.gz -thr 0.9 -uthr 1.1 <temp/>labels_CSF_selection.nii.gz")
        jobs.append("fslmaths <temp/>labels_joined_upsampled_NN.nii.gz -thr 1.9 -uthr 2.1 <temp/>labels_WM_selection.nii.gz")

        jobs.append("fslmaths <temp/>labels_CSF_selection.nii.gz -kernel box 1 -dilF <temp/>labels_CSF_selection_smear.nii.gz -odt int")  # 1 mm (not voxel) around target CSF layer (because recon-all yields 1 mm precision)
        jobs.append("fslmaths <temp/>labels_WM_selection.nii.gz -kernel box 1 -dilF <temp/>labels_WM_selection_smear.nii.gz -odt int")
        jobs.append("fslmaths <temp/>labels_CSF_selection_smear.nii.gz -mul <temp/>labels_WM_selection_smear.nii.gz <temp/>labels_WM_CSF_overlap.nii.gz -odt int")
        jobs.append("fslmaths <temp/>labels_WM_CSF_overlap.nii.gz -div 6 <temp/>labels_WM_CSF_overlap.nii.gz -odt double") # get between 0 and 1
        jobs.append("fslmaths <temp/>labels_WM_CSF_overlap.nii.gz -abs -bin <temp/>labels_no_smoothing_mask.nii.gz -odt int")  # no smoothing possible WM and CSF too close
        jobs.append("fslmaths <temp/>labels_no_smoothing_mask.nii.gz -sub 1 -mul -1 <temp/>labels_smoothing_mask.nii.gz")  # smoothing possible here, WM and CSF distanced more than 1mm

        # mix smoothed and non-smoothed labels => <output/>labels_partially_smoothed.nii.gz
        jobs.append("fslmaths <temp/>labels_smoothed_joined.nii.gz -mas <temp/>labels_smoothing_mask.nii.gz <temp/>labels_smoothed_part.nii.gz -odt int")
        jobs.append("fslmaths <temp/>labels_joined_upsampled_NN.nii.gz -mas <temp/>labels_no_smoothing_mask.nii.gz <temp/>labels_not_smoothed_part.nii.gz -odt int")
        jobs.append("fslmaths <temp/>labels_smoothed_part.nii.gz -add <temp/>labels_not_smoothed_part.nii.gz <partially_smoothed> -odt int")
    else:
        jobs.append("cp <temp/>labels_smoothed_joined.nii.gz <partially_smoothed>") # jsut use the smoothed version.




    # now create rim file
    jobs.append("fslmaths <partially_smoothed> -thr 0.9 -uthr 1.1 <temp/>labels_CSF_smoothed.nii.gz")
    jobs.append("fslmaths <partially_smoothed> -thr 2.9 -uthr 3.1 <temp/>labels_GM_smoothed.nii.gz")
    jobs.append("fslmaths <partially_smoothed> -thr 1.9 -uthr 2.1 <temp/>labels_WM_smoothed.nii.gz")

    jobs.append("fslmaths <partially_smoothed> -kernel sphere <mm> -fmean <temp/>labels_smoothed_smeared.nii.gz -odt double")                   # sphere of 1 voxel radius
    jobs.append("fslmaths <temp/>labels_smoothed_smeared.nii.gz -thr 1.1 -uthr 1.9 <temp/>labels_CSF_candidate.nii.gz -odt double")             # mix between 1 (CSF) and 3 (GM)
    jobs.append("fslmaths <temp/>labels_smoothed_smeared.nii.gz -thr 2.1 -uthr 2.9 <temp/>labels_WM_candidate.nii.gz -odt double")              # mix betwen 2 (WM) and 3 (GM)
    jobs.append("fslmaths <temp/>labels_CSF_smoothed.nii.gz -mas <temp/>labels_CSF_candidate.nii.gz <temp/>labels_CSF_rim.nii.gz -odt double")  # only take part outside of GM
    jobs.append("fslmaths <temp/>labels_WM_smoothed.nii.gz -mas <temp/>labels_WM_candidate.nii.gz <temp/>labels_WM_rim.nii.gz -odt double")     # only take part outside of GM
    jobs.append("fslmaths <temp/>labels_GM_smoothed.nii.gz -add <temp/>labels_CSF_rim.nii.gz -add <temp/>labels_WM_rim.nii.gz <LN2_label_file> -odt int")


    

    for job in jobs:
        job = job.replace("<temp/>", temp_path)
        job = job.replace("<output/>", output_path)
        job = job.replace("<maskWM>", maskWM)
        job = job.replace("<maskGM>", maskGM)
        job = job.replace("<maskCSF>", maskCSF)
        job = job.replace("<mm>", isotropic)
        job = job.replace("<partially_smoothed>", partially_smoothed)
        job = job.replace("<LN2_label_file>", LN2_label_file)

        run (job)

    checkFiles([LN2_label_file])
    return LN2_label_file


def preprocessing_insights (functional_preprocessed, functional_original, biasfield):
    # provides insights in folder "preprocessing_insights"
    # typical use: reprocessing_insights (functional_current, functional, biasfield)



    jobs = []

    jobs.append("rm -rf -f <output/>../preprocessing_insights")
    jobs.append("mkdir <output/>../preprocessing_insights")
    jobs.append("fslroi " + functional_original     + " <temp/>functional_original_first_10_volumes.nii.gz 0 10")
    jobs.append("fslroi " + functional_preprocessed + " <temp/>functional_preprocessed_first_10_volumes.nii.gz 0 10")
    jobs.append("fslmaths <temp/>functional_original_first_10_volumes.nii.gz -Tmean <output/>../preprocessing_insights/functional_original.nii.gz")
    jobs.append("fslmaths <temp/>functional_preprocessed_first_10_volumes.nii.gz -Tmean <output/>../preprocessing_insights/functional_preprocessed.nii.gz")
    jobs.append ("cp " + biasfield + " <output/>../preprocessing_insights/bias_field.nii.gz")

    for job in jobs:
        job = job.replace ("<output/>",         output_path)
        job = job.replace("<temp/>",            temp_path)
        run(job)



def WM_GM_CSF_CTX_masks (subjectID, anatomical_functionalspace_bet, DC=None):
    # typical use: [maskWM, maskGM, maskCSF, cortex_mask_files] = WM_GM_CSF_CTX_masks(subjectID, T1w_functionalspace_bet)
        # cortex_mask_files is a dict, e.g., cortex_mask_files["ctx-rh-postcentral"] yields the respective file name
        # freesurfer_results is folder from recon-all, e.g., "/home/user/freesurfer/subject-574/"

    checkFiles ([anatomical_functionalspace_bet])

    
    current_datetime = datetime.datetime.now()
    myid = "subj_" + subjectID # + "_" + str(current_datetime.strftime("%Y-%m-%d_%H.%M"))



    if DC == None:
        maskWM = mask_path + "mask_wm.nii.gz"
        maskGM = mask_path + "mask_gm.nii.gz"
        maskOuter = mask_path + "mask_csf.nii.gz"
    else:
        maskWM = mask_path + "DC_mask_wm.nii.gz"
        maskGM = mask_path + "DC_mask_gm.nii.gz"
        maskOuter = mask_path + "DC_mask_csf.nii.gz"
    freesurfer_results = "~/freesurfer/" + myid + "/"







    jobs = []
    #if not checkFiles([freesurfer_results + "mri/aparc+aseg.mgz"], True):
    if False:
        # get intersection-free masks in functional space for: WM, GM, CSF (outer)
        jobs.append("mri_convert <anatomical_functionalspace_bet> <temp/>anatomical_bet_functionalspace.mgh")
        jobs.append("mkdir ~/freesurfer") # location for results of recon-all (avoid problems with access rights by locating this within the home directory)
        jobs.append("rm -rf -f ~/freesurfer/<myid>") # neccessary otherwise recon-all will stop if directory exists
        jobs.append("recon-all -i <temp/>anatomical_bet_functionalspace.mgh -sd ~/freesurfer -subjid <myid> -all")
    else:
        print("skipping recon-all: found previous results")

    #if not checkFiles([temp_path + "aparc+aseg_functionalspace.nii.gz"], True):
    if True:
        jobs.append("mri_convert ~/freesurfer/<myid>/mri/aparc+aseg.mgz <temp/>aparc+aseg_talairach.nii.gz")
        jobs.append("mri_convert ~/freesurfer/<myid>/mri/brainmask.mgz <temp/>brainmask_talairach.nii.gz")      # need to match this with anatomical-bet image


        # convert from talairach_space to functional_space by aligning brainmask.mgz with anatomical_bet_functionalspace
        jobs.append("bash antsRegistrationSyN.sh -d 3 -f <anatomical_functionalspace_bet> -m <temp/>brainmask_talairach.nii.gz -n <number_cores> -t s -o <temp/>ANTsOutput_talairach") # takes 90 minutes (if not unneccesarily upsampled)
        jobs.append("rm <output/>matrix_talairach_to_functional.mat")
        jobs.append("mv <temp/>ANTsOutput_talairach0GenericAffine.mat <output/>matrix_talairach_to_functional.mat")
        jobs.append("antsApplyTransforms -d 3 -i <temp/>aparc+aseg_talairach.nii.gz -r <anatomical_functionalspace_bet> -o <temp/>aparc+aseg_functionalspace.nii.gz -n NearestNeighbor -t <output/>matrix_talairach_to_functional.mat")
    else:
        print("skipping another part of recon. already got results.")





    # remove cerebellum and brainstem: currently not needed (main pipeline: only selected cortex areas, DC: use schaefer100_functoinalspace to mask layers and remove cerebellum)
    if False:
        intensities_cerebellum = [7, 8, 15, 16, 46, 47]
        list_cerebellum_files=[]
        for intensity in intensities_cerebellum:
            file = "<temp/>aparc+aseg_functionalspace_" + str(intensity) + ".nii.gz"
            jobs.append ("rm " + file)
            list_cerebellum_files.append(file)
            job = ("fslmaths <temp/>aparc+aseg_functionalspace.nii.gz -thr <intensity> -uthr <intensity> <file> -odt double")
            job = job.replace("<intensity>", str(intensity))
            job = job.replace("<file>", file)
            jobs.append(job)


        job = ""
        cerebellum_mask = "<temp/>cerebellum_mask.nii.gz"
        for n, file in enumerate(list_cerebellum_files):
            if n == 0:
                jobs.append ("rm " + file)
                job = "fslmaths " + file
            else:
                job = " -add " + file
        job += " " + cerebellum_mask
        jobs.append(job)

        segementation_no_cerebellum = "<temp/>aparc+aseg_functionalspace_no-cerebellum.nii.gz"
        jobs.append("rm " + segementation_no_cerebellum)
        jobs.append("fslmaths <temp/>aparc+aseg_functionalspace.nii.gz -mas " + cerebellum_mask + " " + segementation_no_cerebellum)






    # process masks
    if True: # keep cerebellum and brainstem
        jobs.append("mri_binarize --i <temp/>aparc+aseg_functionalspace.nii.gz --gm --o <maskGM>")
        jobs.append("mri_binarize --i <temp/>aparc+aseg_functionalspace.nii.gz --all-wm --o <maskWM>")
    else: # remove cerebellum and brainstem
        jobs.append("mri_binarize --i " + segementation_no_cerebellum + " --gm --o <maskGM>")
        jobs.append("mri_binarize --i " + segementation_no_cerebellum + " --all-wm --o <maskWM>")


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
    locations["ctx-rh-postcentral"]         = 2022  # somatosensory cortex
    locations["ctx-rh-precentral"]          = 2024  # primary motor
    locations["ctx-rh-superiortemporal"]    = 2030  # primary auditory

    for locR in list(locations.keys()): # create entries for other hemisphere
        locL = locR.replace("ctx-rh", "ctx-lh")
        locations[locL] = locations[locR] - 1000

    cortex_mask_files = {}

    for loc, intensity in locations.items():
        job = "fslmaths <temp/>aparc+aseg_functionalspace.nii.gz -thr <intensity> -uthr <intensity> <CTXmaskFile> -odt double"
        CTXmaskFile = mask_path + "mask_" + loc + ".nii.gz"
        job = job.replace("<CTXmaskFile>", CTXmaskFile)
        job = job.replace ("<loc>", loc)
        job = job.replace ("<intensity>", str(locations[loc]))
        jobs.append(job)
        cortex_mask_files[loc] = CTXmaskFile




    for job in jobs:
        job = job.replace("<anatomical_functionalspace_bet>", anatomical_functionalspace_bet)
        job = job.replace("<myid>", myid)
        job = job.replace("<maskGM>", maskGM)
        job = job.replace("<maskWM>", maskWM)
        job = job.replace("<maskOuter>", maskOuter)
        job = job
        job = job.replace("<temp/>", temp_path)
        job = job.replace("<output/>", output_path)
        job = job.replace("<number_cores>", str(number_cores)) # defined in main

        run(job)



    return [maskWM, maskGM, maskOuter, cortex_mask_files]




def dc_layers_prepare_masks(DC_mean_bold, DC_mean_bold_GM_mask, DC_anat_GM_mask, DC_anat_WM_mask, DC_anat_CSF_mask, DC_anatomical):
    # typical use: [DC_bold_space_GM_mask, DC_bold_space_WM_mask, DC_bold_space_CSF_mask, DC_anatomical_functional_space] = dc_layers_prepare_masks(DC_mean_bold, DC_mean_bold_GM_mask, DC_GM_mask, DC_WM_mask, DC_CSF_mask, DC_anatomical)
    # DC-specific function: create Gm, WM and CSF masks in DC_bold_space such that they match the CPAC DCW file.

    # steps:
    #   1. create transformation from DC_anatomical_space to DC_bold_space by rigid transform (rotion, translation) between given GM masks in both spaces
    #   2. create GM, WM, CSF masks in DC_bold_space from matrix resulting from (1) and given masks in DC_anatomical_space

    output_WM                           = output_path + "DC_functionalspace_WM.nii.gz"
    output_GM                           = output_path + "DC_functionalspace_GM.nii.gz"
    output_CSF                          = output_path + "DC_functionalspace_CSF.nii.gz"
    DC_anatomical_functional_space      = output_path + "DC_anatomical_functionalspace.nii.gz"


    jobs = []

    jobs.append("rm <temp/>DC_ANTsOutput0GenericAffine.mat")
    jobs.append("bash antsRegistrationSyN.sh -d 3 -f <cpac_bold_GM.nii.gz> -m <cpac_anat_GM.nii.gz> -n <number_cores> -o <temp/>DC_ANTsOutput -t r") # rigid transformation: within same subject
    # parameters:
    #               f: fixed image
    #               t: 'r' (rigid), 'a': rigid + affine (2 stages), 's': rigid + affine + deformable syn (3 stages)
    #               d: dimension (2 or 3, not 4)
    #               m: moving image (that needs transform)
    #               o: outputfile
    # output:
    # [outputname]Affine.mat
    # [outputname]Warped.nii.gz

    jobs.append("rm <output/>matrixDC-anatomical_to_functional.mat")
    jobs.append("mv <temp/>DC_ANTsOutput0GenericAffine.mat <output/>matrixDC-anatomical_to_functional.mat")

    jobs.append("antsApplyTransforms -d 3 -i <cpac_anat_GM.nii.gz> -r <cpac_bold.nii.gz> -o <output_GM.nii.gz> -n NearestNeighbor -t <output/>matrixDC-anatomical_to_functional.mat")
    jobs.append("antsApplyTransforms -d 3 -i <cpac_anat_WM.nii.gz> -r <cpac_bold.nii.gz> -o <output_WM.nii.gz> -n NearestNeighbor -t <output/>matrixDC-anatomical_to_functional.mat")
    jobs.append("antsApplyTransforms -d 3 -i <cpac_anat_CSF.nii.gz> -r <cpac_bold.nii.gz> -o <output_CSF.nii.gz> -n NearestNeighbor -t <output/>matrixDC-anatomical_to_functional.mat")
    jobs.append("antsApplyTransforms -d 3 -i <DC_anatomical> -r <cpac_bold.nii.gz> -o <DC_anatomical_functional_space> -n Linear -t <output/>matrixDC-anatomical_to_functional.mat")







    for job in jobs:
        job = job.replace("<temp/>", temp_path)
        job = job.replace("<output/>", output_path)
        job = job.replace("<cpac_bold_GM.nii.gz>", DC_mean_bold_GM_mask)
        job = job.replace("<cpac_anat_GM.nii.gz>", DC_anat_GM_mask)
        job = job.replace("<cpac_anat_WM.nii.gz>", DC_anat_WM_mask)
        job = job.replace("<cpac_anat_CSF.nii.gz>", DC_anat_CSF_mask)
        job = job.replace("<cpac_bold.nii.gz>", DC_mean_bold)
        job = job.replace("<output_GM.nii.gz>", output_GM)
        job = job.replace("<output_WM.nii.gz>", output_WM)
        job = job.replace("<output_CSF.nii.gz>", output_CSF)
        job = job.replace("<DC_anatomical>", DC_anatomical)
        job = job.replace("<DC_anatomical_functional_space>", DC_anatomical_functional_space)
        job = job.replace("<number_cores>", str(number_cores))  # defined in mainf



        run(job)

    checkFiles([output_GM, output_WM, output_CSF, DC_anatomical_functional_space])
    return [output_GM, output_WM, output_CSF, DC_anatomical_functional_space]

def DC_atlas_to_functional_space(atlasfile, DC_anatomical_functional_space, DC_functional, DC_DCW):
    # typical use: [schaefer100_DC_functionalspace, DC_DCW_MNIspace, matrix_MNI_to_DCfunctional] = DC_atlas_to_functional_space(constant_schaefer100, DC_anatomical_functional_space, DC_functional, DC_DCW)
    # output: schaefer100 atlas in DC-functional-space



    DC_anatomical_functional_space_bet      = temp_path + "DC_anatomical_functional_space_bet.nii.gz"
    matrix_MNI_to_DC_functional             = output_path + "DC_matrix_MNI_to_DC_functional.mat"
    atlasfile_DC_functionalspace            = output_path + "DC_schaefer100_functionalspace.nii.gz"
    matrix_DC_functional_to_MNI             = output_path + "DC_matrix_DC_functional_to_MNI.mat"
    DC_DCW_MNIspace                         = output_path + "DC_DCW_MNIspace.nii.gz"


    jobs = []

    # DC-anatomical_functionalspace: bet
    jobs.append("mri_synthstrip -i <DC_anatomical_functional_space> --no-csf -o <DC_anatomical_functional_space_bet>")

    # matrix MNI -> DC-functional space
    jobs.append("bash antsRegistrationSyN.sh -d 3 -f <DC_anatomical_functional_space_bet> -m <MNI> -n <number_cores> -t s -o <temp/>DC_ANTsOutput") # deformable transformation: between subject and template
    # parameters:
    #               f: fixed image
    #               t: 'r' (rigid), 'a': rigid + affine (2 stages), 's': rigid + affine + deformable syn (3 stages)
    #               d: dimension (2 or 3, not 4)
    #               m: moving image (that needs transform)
    #               o: outputfile
    # output:
    # [outputname]Affine.mat
    # [outputname]Warped.nii.gz
    jobs.append("rm <matrix_MNI_to_DC_functional>")
    jobs.append("mv <temp/>DC_ANTsOutput0GenericAffine.mat <matrix_MNI_to_DC_functional>")

    # schaefer-100 to  -> DC-functional space
    jobs.append("antsApplyTransforms -d 3 -i <atlasfile> -r <DC_functional> -o <atlasfile_DC_functionalspace> -n NearestNeighbor -t <matrix_MNI_to_DC_functional>") # NN
    # parameters:
    #               -t transformation matrix (can use multiple, each preceded by a -t)
    #               -i input
    #               -o output
    #               -r reference image: "For warping input images, the reference image defines the spacing, origin, size, and direction of the output warped image."
    #



    # functional to MNI (DC_anatomical_functionalspace to MNI) and then DC_DCW to MNI (for later averaging between subjects: "python3 main.py dc summary")
    jobs.append("bash antsRegistrationSyN.sh -d 3 -f <MNI> -m <DC_anatomical_functional_space_bet> -n <number_cores> -o <temp/>DC_ANTsOutput -t s")  # deformable transformation: between subject and template
    jobs.append("rm <matrix_DC_functional_to_MNI>")
    jobs.append("mv <temp/>DC_ANTsOutput0GenericAffine.mat <matrix_DC_functional_to_MNI>")
    jobs.append("antsApplyTransforms -d 3 -i <DC_DCW> -r <MNI> -o <DC_DCW_MNIspace> -n Linear -t <matrix_DC_functional_to_MNI>") # Linear transformation





    
    
    for job in jobs:
        job = job.replace("<temp/>",                                temp_path)
        job = job.replace("<output/>",                              output_path)
        job = job.replace("<atlasfile>",                            atlasfile)
        job = job.replace("<DC_anatomical_functional_space>",       DC_anatomical_functional_space)
        job = job.replace("<DC_anatomical_functional_space_bet>",   DC_anatomical_functional_space_bet)
        job = job.replace("<matrix_MNI_to_DC_functional>",          matrix_MNI_to_DC_functional)
        job = job.replace("<matrix_DC_functional_to_MNI>",          matrix_DC_functional_to_MNI)
        job = job.replace("<MNI>",                                  constant_MNItemplate)
        job = job.replace("<DC_functional>",                        DC_functional)
        job = job.replace("<atlasfile_DC_functionalspace>",         atlasfile_DC_functionalspace)
        job = job.replace("<DC_DCW>",                               DC_DCW)
        job = job.replace("<DC_DCW_MNIspace>",                      DC_DCW_MNIspace)
        job = job.replace("<number_cores>",                         str(number_cores))  # defined in main

        run (job)


    checkFiles([atlasfile_DC_functionalspace, DC_DCW_MNIspace, matrix_MNI_to_DC_functional])
    return [atlasfile_DC_functionalspace, DC_DCW_MNIspace, matrix_MNI_to_DC_functional]


def create_LN2_profile (activity_map, labels, config_title, number_layers, short_name=None):
    # typical use: layer_profile = create_LN2_profile (title, DC_dcw, "7layers_equist", 7)
        # output layer_profile is a dict: {"PNG":"<output>.png", "table":[intensity_layer0, intensity_layer1, ....]}

    jobs = []

    LN2_output_folder = output_path + "LN2_PROFILES"
    jobs.append("mkdir " + LN2_output_folder)
    LN2_output_folder += "/"
    title = os.path.basename(activity_map)
    if short_name == None:
        LN2_output_folder += title.replace(".nii.gz", "")   # if no short_name is given, output will be named after activity map
    else:
        LN2_output_folder += short_name         # e.g., "main" or "motor_left"
    jobs.append("mkdir " + LN2_output_folder)
    LN2_output_folder += "/"

    output_txt = LN2_output_folder + (os.path.basename(activity_map)).replace(".nii.gz", "") + "_output_" + config_title + ".txt"
    #                  - Column 1 is the layer number.
    #                  - Column 2 is the mean signal in this layer.
    #                  - Column 3 is the STDEV of the signal variance across all voxels in this layer.
    #                  - Column 4 is the number of voxels per layer.



    jobs.append("rm <output_txt>")
    jobs.append("LN2_PROFILE -input <activity_map> -layers <labels> -plot -output <output_txt>")

    for job in jobs:
        job = job.replace("<activity_map>", activity_map)
        job = job.replace("<labels>", labels)
        job = job.replace("<output_txt>", output_txt)

        run(job)

    ###########################
    # read out the resulting text file (intensities in second column, layer number starting at 1 in first column; cortical depth is increasing with layer number: CSF -> WM)

    checkFiles([output_txt])    # resulting file always exists, but may be empty in case LN2_PROFILE has found no sginal (e.g., when using small ROIs). Also, one or more lines may be missing.
    #mean_signal = np.loadtxt(output_txt, usecols=(1))




    ######################################
    # read in the LN2_PROFILE text output, make it robust, because different cases may occur
    #   - a value is "-nan"
    #   - one or all lines are missing
    #   - the file is empty
    # add dummy lines to the file, such that it can be read, even though it is empty or some lines (i.e., layer information) is missing. This occurs for small ROIs when some layer bins are empty.

    output_txt_padding = output_txt.replace(".txt", "_padding.txt")
    run ("rm " + output_txt_padding)
    run ("cp " + output_txt + " " + output_txt_padding)
    with open(output_txt_padding, "a+") as file:
        for n in range(number_layers):
            file.write("-1\t-1\t-1\n")

    mean_signal_unprocessed = np.loadtxt(output_txt_padding, usecols=(0, 1)) # e.g. col1, line 3 = mean_signal[2][0]
    mean_signal_processing = {}
    mean_signal = []
    for line in mean_signal_unprocessed:
        layer = line[0]
        signal = line[1]
        if not ("nan" in str(signal) or layer == -1): # two cases: a value may be "-nan", or it layer number is "-". The first case is a result of LN2_PROFILE. The second case is the result of the previous padding and missing lines.
            mean_signal_processing[layer] = signal # add to dict and check later which values are avilable

    for n in range(1, number_layers+1): # first line is layer 1
        if n in mean_signal_processing:
            mean_signal.append(mean_signal_processing[n]) # entry was available in LN2_PROFILE output
        else:
            mean_signal.append(0) # entry was unavailable
    #################################################################### # mean_signal is now a list of number_layers values (with 0.0 in case of errors or empty ROIs)










    if config_title == "temp": # do not output any files
        return {"table": mean_signal}
    else:

        # Plot the data
        cortical_depth = np.arange(1, len(mean_signal) + 1)     #   xvals
        list_yvals = [mean_signal]                              #   yvals (here: only one set)
        title = os.path.basename(activity_map)
        output_png = output_txt.replace(".txt", ".png")


        dict_xticks = {}
        dict_xticks['sWM'] = min(cortical_depth)
        dict_xticks['pial'] = max(cortical_depth)

        plot_multiple_graphs(cortical_depth, list_yvals, output_png,
                             axis_title_x       = None,
                             axis_title_y       = "Mean Signal",
                             list_labels        = None,
                             plot_title=        "layer profile for " + title,
                             dict_xticks=dict_xticks)






        return {"PNG":output_png, "TXT":output_txt, "table":mean_signal}







def LN2_layer_profiles_ROIwise (main_image, label_file, atlas_imagespace, numberROIs, numberLayers, configuration):
    # typical use: layer_profile_ROI_wise = LN2_layer_profiles_ROIwise (main_image, label_file, schaefer100_DC_dunctionalspace, numberROIs=100, numberLayers=3, configuration)
    #   input:
    #       main_image:         main image (e.g., DC_dcw or activity map); atlas and image file must be aligned in the same space
    #       label_file          label file, used as an input for LN2_LAYERS, that includes CSF, WM, and GM masks. E.g., created by "create_LN2_labels"
    #       atlas_imagespace    atlas and image file must be aligned in the same space
    #       numberROIs:         needs to match the atlas
    #       numberLayers:       typically, 3
    #       configuration:      typically, "layerbins_equivol"
    #   output:
    #       list of lists: layer_profile_ROI_wise[0..numberLayers-1][0..numberROIs-1]       typically as input for a subsequent visualisation
    
    layer_profile_ROI_wise = np.zeros((numberLayers, numberROIs), dtype=float)

    # calculate layers with LN2_LAYERS for given number of layers
    layer_results = create_LN2_layers ("DC_layer_profile_ROI_wise", label_file, numberLayers, anatomical=None) # do not set anatomical image because no plot is wanted here
    selected_layers = layer_results[configuration] # typically, configuration="layerbins_equivol"


    for ROI in range(numberROIs):
        text = "\n\nprocessing layers for ROI: " + str(ROI)
        print(text)
        if simulation_switch(-1) == 0:
            run_log(text)

        selected_layers_currentROI = atlas_ROI_masking (selected_layers, atlas_imagespace, ROI, temp_path + "temp_selected_layers_currentROI.nii.gz") # by applying a mask on the LN2_LAYERS output file, a specific region of interest is selected

        layer_profile = create_LN2_profile(main_image, selected_layers_currentROI, "temp", numberLayers) # marked as temporary, such that not output files (PNG, TXT) are created.
        layer_profile_table = layer_profile["table"]  # list of mean signals for different layers as defined in numberLayers (starting from 0, i.e., closer to white matter)


        #print("\n\n ROI: " + str(ROI) + " - " + str(layer_profile_table))

        for layer_number, meanSignal in enumerate(layer_profile_table):
            layer_profile_ROI_wise[layer_number][ROI] = meanSignal




    return layer_profile_ROI_wise


def atlas_select_ROI (atlas, intensity, ROI_output):   # selects a ROI by setting a target intensity in an atlas, ROI_output determins the output file
    job = "fslmaths <atlas> -thr <intensity> -uthr <intensity> <ROI_output>"
    job = job.replace("<atlas>",                atlas)
    job = job.replace("<intensity>",            str(intensity))
    job = job.replace("<ROI_output>",           ROI_output)

    run(job)
    return ROI_output


def atlas_ROI_masking (image, atlas_imagespace, ROI, outputfile):
    # selects a region from an image by using a mask based on an ROI (i.ie, an intensity) and an atlas
    ROI_mask = temp_path + "temp_ROI.nii.gz"

    atlas_select_ROI (atlas_imagespace, ROI, ROI_mask)
    job = "fslmaths <image> -mas <ROI_mask> <outputfile>"
    job = job.replace("<atlas_imagespace>",     atlas_imagespace)
    job = job.replace("<ROI_mask>",             ROI_mask)
    job = job.replace("<outputfile>",           outputfile)
    job = job.replace ("<image>",               image)

    run (job)
    return outputfile


def apply_mask (image, mask, suffix=None, erode_voxels=False):
    # erode_voxels: e.g. 1 to shrink the mask by 1 voxel
    if suffix == None:
        suffix = "mask"

    outputfile = naming (image, suffix, output_path)

    if True: # simple version or with smoothing to preserve original resolution
        mask_imagespace = transform_space(mask, image, sampling="NN")
    else:
        mask_binary = naming(mask, "binary", temp_path)
        run("fslmaths " + mask + " -abs -bin " + mask_binary + "-odt double")

        if erode_voxels:
            run("fslmaths " + mask_binary + " -kernel box 3 -ero " + mask_binary + " -odt double")


        mask_transformed = transform_space(mask_binary, image, sampling="Li")


        mask_imagespace = naming(mask_transformed, "smooth", temp_path)
        run("fslmaths " + mask_transformed + " -thr 0.50 -uthr 1.01 " + mask_imagespace + " -odt double") # re-sampled and smoothed mask
    
    job = "fslmaths <image> -mas <mask> <outputfile> -odt double"
    job = job.replace("<mask>",                 mask_imagespace)
    job = job.replace("<outputfile>",           outputfile)
    job = job.replace ("<image>",               image)

    run (job)
    return outputfile



def apply_mask_with_margin (image, mask, title):
    # typical use: functional_bet_GMmargin = apply_mask_with_margin (functional_bet, maskGM, "GMmargin")
    # Applies a mask with a margin to an image. Typically used to reduce the functional image to the GM area to speed-up the creation of an activity map with fsl-feat.
    output_img = naming(image, title, output_path)
    output_img_insight = naming(output_img, "insight", output_path)
    mask_dil = naming(mask, "dil", temp_path)

    jobs = []
    jobs.append("")

    jobs.append("fslmaths <mask> -dilF -kernel box 3.0 -fmean <mask_dil> -odt float")       # TODO: check the "-fmean" option. moreover, kerlal before filtering. should likely be "fslmaths <mask> -kernel box 3 -dilF <mask_dil> -odt int"
    #jobs.append("fslmaths <mask> -dilF -kernel 3D <mask_dil> -odt float")                   # 3x3x3 box centered on target voxel
    jobs.append("fslmaths <mask_dil> -abs -bin <mask_dil>")                                 # make sure there are no voxels other than 0 or 1.
    jobs.append("fslmaths <image> -mas <mask_dil> <output_img>")
    jobs.append("fslroi <output_img> <output_img_insight> 0 1")

    for job in jobs:
        job = job.replace("<image>", image)
        job = job.replace("<output_img>", output_img)
        job = job.replace("<mask>", mask)
        job = job.replace("<mask_dil>", mask_dil)
        job = job.replace("<output_img_insight>", output_img_insight)

        run(job)

    return output_img


def create_preview(functional):
    # extracts the first image of a 4D image
    preview_file = functional.replace(".nii.gz", "_first.nii.gz")

    jobs = []
    jobs.append("rm <preview_file>")
    jobs.append("fslroi <functional> <preview_file> 0 1")

    for job in jobs:
        job = job.replace("<functional>", functional)
        job = job.replace("<preview_file>", preview_file)
        run(job)

    return preview_file


def apply_mask_cortex_selection(functional, cortex_mask_files, ctx_list):
    # typical use: functional_bet_ctx_selection = apply_mask_cortex_selection(functional_bet, cortex_mask_files, ["ctx-rh-precentral, ctx-lh-precentral"])
    # Creates reduced version of the functional image for GM in a selected number of cortices. Adds a margin of a few voxels.
    # Resulting masked functional image can be used to run fsl-feat faster than with the complete brain.

    ctx_selection = temp_path + "ctx-selection.nii.gz"

    ctx_filenames = []
    for location, file_ctx_mask in cortex_mask_files.items():
        if location in ctx_list:
            ctx_filenames.append(file_ctx_mask)

    jobs = []
    jobs.append("rm <ctx_selection>")

    # create template for mask from functional image such that comon dimensions are guaranteed. This avoids fslmaths: "Attempted to multiply images of different sizes terminate called after throwing an instance of 'std::runtime_error'"
    insight_file = output_path + "insight_functional_ctx-selection.nii.gz"
    jobs.append("fslroi <functional> <ctx_selection> 0 1")
    jobs.append("rm " + insight_file)
    jobs.append("cp <ctx_selection> " + insight_file)
    jobs.append("fslmaths <ctx_selection> -mul 0 <ctx_selection>")


    for file in ctx_filenames:
       jobs.append("fslmaths <ctx_selection> -max " + file + " <ctx_selection>")



    for job in jobs:
        job = job.replace("<ctx_selection>", ctx_selection)
        job = job.replace("<functional>", functional)
        run(job)

    return apply_mask_with_margin (functional, ctx_selection, "ctx-selection")






def dc_dcw_zscore_plots(DC_dcw, DC_anatomical_functional_space, DC_bold_space_GM_mask, result_folder, title=None):

    ###############################################################
    # plot z-score maps with DC map projection

    image = DC_anatomical_functional_space
    projection = DC_dcw
    output_png = result_folder + "DC_overlay.png"


    plot_overlay(None, projection, output_png, title=title, type="activity", view="zscore")



    #######################################################
    # now only cortex:

    img_file = DC_dcw
    img_file_ctx = result_folder + "DC_functionalspace_GM.nii.gz"
    ctx_file = DC_bold_space_GM_mask
    outputCTX_png = result_folder + "DC_overlay_GMmask.png"

    img_file_correct_space = transform_space(img_file, ctx_file, sampling="Li", suffix="zcor-correctspace")

    job = "fslmaths <img_file> -mas <ctx_file> <img_file_ctx>"
    job = job.replace("<img_file>", img_file_correct_space)
    job = job.replace("<ctx_file>", ctx_file)
    job = job.replace("<img_file_ctx>", img_file_ctx)
    run(job)


    image = DC_anatomical_functional_space
    plot_overlay(image, img_file_ctx, outputCTX_png, title=title, type="activity", view="zscore", alpha=1.0)

    return [output_png, outputCTX_png]




def transform_space(input_image, reference_image, sampling="Li", suffix=None):
    # This function transforms an image to another space. Typically the input image (e.g., LN2_labels, cortex masks) is already aligend to the reference image (e.g. functional_reference), but has larger dimensions.
    # This is needed because the functional impage sometimes misses a few centimeters of the cortex (e.g., s548rescan, s618) and the pipeline needs the aligned, but complete anatomical image.
    # Only in later steps, "fslmaths -mas" required identical spaces and, hence, this function.

    if suffix == None:
        description = "transformedspace"
    else:
        description = "transformedspace-" + suffix


    output_image = naming(input_image, description, temp_path)
    
    # sampling can be "Li" (default) or "NN"
    jobs = []
    jobs.append("rm <output_image>")
    jobs.append("3dresample -master <reference_image> -rmode <NN_Li> -prefix <output_image> -input <input_image>")

    for job in jobs:
        job = job.replace("<input_image>", input_image)
        job = job.replace("<reference_image>", reference_image)
        job = job.replace("<output_image>", output_image)
        job = job.replace("<NN_Li>", sampling)

        run (job)

    checkFiles([output_image])
    if ("label" in input_image or "mask" in input_image) and sampling != "NN":
        print("\n\n Warming: You may be trying to use a sampling method other than nearest neigbor on a label or mask file.")
        wait()

    return output_image


def mean_image (list_images, output_image, base_space=None):
    # calculates the mean image from a set of images (in same image space)
    # base_space is an optional input image to determine the image space (e.g. MNI template) whose intentsity will be multiplied by zero

    checkFiles(list_images)

    run ("rm " + output_image)

    cmd = "fslmaths"

    if base_space == None:
        for n, image in enumerate (list_images):
            if (n > 0):
                cmd += " -add"
            cmd += " " + image
    else:
        cmd += " " + base_space + " -mul 0"
        for image in list_images:
            cmd += " -add " + image

    cmd += " -div " + str(len(list_images))  + ".0"
    cmd += " " + output_image
    cmd += " -odt double"
    run(cmd)

    checkFiles([output_image])
    return output_image