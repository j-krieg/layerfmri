# main functions of the pipeline: layerfmri, dc, recon-all


def dc_layers(subject, DC_anatomical, DC_dcw, DC_mean_bold, DC_mean_bold_GM_mask, DC_GM_mask, DC_WM_mask, DC_CSF_mask):
    checkFiles([DC_anatomical, DC_dcw, DC_mean_bold, DC_mean_bold_GM_mask, DC_GM_mask, DC_WM_mask, DC_CSF_mask])

    # input:
    #   - GM, WM and CSF masks in DC_anatomical_space
    #   - only GM mask in DC_bold_space
    #   - DC_dcw in DC_bold_space

    # output:
    #   - layer profile for DC_dcw

    # steps:
    #   1a. create transformation from DC_anatomical_space to DC_bold_space by rigid transform (rotion, translation) between given GM masks in both spaces
    #   1b. create GM, WM, CSF masks in DC_bold_space from matrix resulting from (1) and given masks in DC_anatomical_space
    #   2.  schaefer-100 parcellation in DC-functional space
    #   3.  create label file from the three masks (Gm , WM, CSF)) as an input for laynii
    #   4.  create layers with laynii
    #   5.  calculate layer profile for DC_dcw (in DC_bold_space)
    #   6.  generate LN2 profile data parcellation-wise
    #   7.  visualisation: generate brainspace plots

    simulation_switch(0)

    overview_files = []

    # step 1
    run_log("\n\n" + now() + "\n1. DC: generate WM, GM, CSF masks in DC functional space")

    # if ("574" in label) or ("618" in label): # better results with recon-method (todo: check why). (GM mask method: 548re, 549re, 5888re, 599re; recon method: 574re, 618, 599re)
    if False:
        if False:
            # DC_anatomical_functional_space = ... output_path + "DC_anatomical_functionalspace.nii.gz"
            [maskWM, maskGM, maskCSF, cortex_mask_files] = WM_GM_CSF_CTX_masks("DC_" + subject, DC_anatomical_functional_space, DC="DC")  # uses recon-all
        else:
            # re-use from previous runs
            DC_anatomical_functional_space = output_path + "DC_anatomical_functionalspace.nii.gz"
            maskGM = output_path + "masks/DC_mask_gm.nii.gz"
            maskWM = output_path + "masks/DC_mask_wm.nii.gz"
            maskCSF = output_path + "masks/DC_mask_csf.nii.gz"
            checkFiles([DC_anatomical_functional_space, maskGM, maskWM, maskCSF])

    else:
        [maskGM, maskWM, maskCSF, DC_anatomical_functional_space] = dc_layers_prepare_masks(DC_mean_bold, DC_mean_bold_GM_mask, DC_GM_mask, DC_WM_mask, DC_CSF_mask, DC_anatomical)

    simulation_switch(1)
    run_log(
        "\n\n" + now() + "\n2. DC: registration of MNI template to DC_anatomical_functionalspace_bet and transformation of schaefer100 atlas into DC-functional space")
    [schaefer100_DC_functionalspace, DC_dcw_MNIspace, matrix_MNI_to_DC_functional] = DC_atlas_to_functional_space(
        constant_schaefer100, DC_anatomical_functional_space, DC_mean_bold, DC_dcw)
    # "DC_dcw_MNIspace" will be used later by "DC_summary(subjects)" to create an average DC map from all subjects in MNI space.
    # matrix MNI_to_DCfunctional neeeded to use freesurfer results later.
    simulation_switch(0)

    voxelsize_mm = 0.50  # for label file and images/masks

    # step 3
    run_log("\n\n" + now() + "\n3. DC: create LN2 label file from GM, WM, and CSF masks")
    title = datasets[subject]["label"] + "DC_CPAC-dcw"
    label_file = create_LN2_labels(voxelsize_mm, title, maskWM, maskGM, maskCSF)
    overview_files.append(label_file)

    # 2 purposes: 1. upsampling for better results, 2. ensure that masks have same resulutio as the images
    run_log("\n\n" + now() + "\n2c. DC: upsampling")
    [DC_dcw, schaefer100_DC_functionalspace, DC_anatomical_functional_space, maskWM, maskGM, maskCSF] = resampling(voxelsize_mm, [DC_dcw, schaefer100_DC_functionalspace, DC_anatomical_functional_space, maskWM, maskGM, maskCSF])
    # from now on, the upsampled files are used.
    overview_files.append(schaefer100_DC_functionalspace)
    overview_files.append(DC_dcw)
    overview_files.append(DC_anatomical_functional_space)

    run("fslmaths " + schaefer100_DC_functionalspace + " -kernel box 1 -dilF " + temp_path + "schaefer100_functional_dil.nii.gz -odt int")  # dilate by 1mm
    run("fslmaths " + maskGM + " -kernel box 1 -dilF " + temp_path + "maskGM_1mmdil.nii.gz -odt int")
    cerebral_GM = apply_mask(temp_path + "maskGM_1mmdil.nii.gz", temp_path + "schaefer100_functional_dil.nii.gz", suffix="cerebral")  # typically: remove cerebellum and brainstem with schaefer100_functionalspace mask (same done with layer files before calculating the profile)
    label_file_cerebral = apply_mask(label_file, cerebral_GM, suffix="cerebral-GM")

    # step 3b: create plots
    [plot_DC, plot_DC_GM] = dc_dcw_zscore_plots(DC_dcw, DC_anatomical_functional_space, cerebral_GM, output_path, label)
    overview_files.append(plot_DC)
    overview_files.append(plot_DC_GM)
    print("zscore done.\n\n")

    # step 4
    run_log("\n\n" + now() + "\n4. DC: calculating the different layers with LN2_layers")
    number_layers = 7
    configuration = "layerbins_equivol"
    layer_results = create_LN2_layers(title, label_file_cerebral, number_layers, DC_anatomical_functional_space)
    selected_layers = layer_results[configuration]

    overview_files.append(mricrogl_plot_labelfile(DC_anatomical_functional_space, label_file, orientation='x'))
    overview_files.append(mricrogl_plot_labelfile(DC_anatomical_functional_space, label_file, orientation='y'))
    overview_files.append(mricrogl_plot_labelfile(DC_anatomical_functional_space, label_file, orientation='z'))

    if True:  # demo of different files for the report

        # demo of layer creation pipeline: anatomical, 1mm freesurfer segmentation, upsamling and smoothing, 1-voxel rim, LN2_LAYERS result
        GM_unprocessed = temp_path + "LN2_labels_layer2_1mm.nii.gz"  # has intensity 2
        GM_unprocessed3 = temp_path + "LN2_labels_layer2_1mm_demo.nii.gz"
        run("fslmaths " + GM_unprocessed + " -abs -bin " + GM_unprocessed3, no_sim=True)
        run("fslmaths " + GM_unprocessed3 + " -mul 3 " + GM_unprocessed3, no_sim=True)

        cuts = [0.66]
        orientation = "coronal"
        demo_files = []
        # demo_files.append(mricrogl_plot_labelfile(DC_anatomical_functional_space, None, output_filename=output_path + "segmentation01_anat.png",            cuts=cuts, orientation=orientation))
        # demo_files.append(mricrogl_plot_labelfile(DC_anatomical_functional_space, GM_unprocessed3, output_filename=output_path + "segmentation02_unsmoothed-1mm.png",  cuts=cuts, orientation=orientation))
        demo_files.append(mricrogl_plot_labelfile(DC_anatomical_functional_space, temp_path + "labels_smoothed_joined.nii.gz", output_filename=output_path + "segmentation03_smoothed.png", cuts=cuts, orientation=orientation))
        demo_files.append(mricrogl_plot_labelfile(DC_anatomical_functional_space, label_file, output_filename=output_path + "segmentation04_label.png", cuts=cuts, orientation=orientation, skip_GM=True))
        demo_files.append(mricrogl_plot_layerfile(DC_anatomical_functional_space, selected_layers, number_layers, output_filename=output_path + "segmentation05_layer.png", cuts=cuts, orientation=orientation))

        demo_files.append(stack_images(demo_files, output_path + "segmentation_demo.png", direction="horizontal", overlap=0.1))

        overview_files.extend(demo_files)

    overview_files.append(mricrogl_plot_layerfile(DC_anatomical_functional_space, selected_layers, number_layers))

    copy_overview_files("DC_" + subject, overview_files, no_replace=True)

    overview_files.append(label_file)
    overview_files.append(label_file.replace(".nii.gz", "_ortho.png"))
    overview_files.append(selected_layers)
    overview_files.append(selected_layers.replace(".nii.gz", "_ortho.png"))

    # step 5
    run_log("\n\n" + now() + "\n5. DC: calculating the resulting layer profile")
    config_title = str(number_layers) + "_" + configuration
    output_PNG = output_path + "LN2_PROFILES/DC_" + str(number_layers) + "_layers_" + configuration + "_s" + subject + ".png"
    output_TXT = output_PNG.replace(".png", ".txt")
    layer_profile = create_LN2_profile(DC_dcw, selected_layers, config_title, number_layers)
    layer_profile_PNG = layer_profile["PNG"]
    layer_profile_TXT = layer_profile["TXT"]
    layer_profile_table = layer_profile["table"]  # list of mean signals for different layers as defined in number_layers (starting from 0, i.e., closer to white matter)
    run("mv " + layer_profile_PNG + " " + output_PNG)
    run("mv " + layer_profile_TXT + " " + output_TXT)
    np.savetxt(output_path + "DC_profile_full_" + str(number_layers) + "layers.txt", layer_profile_table, delimiter=" ")

    overview_files.append(output_PNG)
    overview_files.append(output_TXT)
    copy_overview_files("DC_" + subject, overview_files)

    ###############################################
    # create region-wise brainspace plots
    # uses "schaefer100_DC_functionalspace" that has already been created before the upsampling step

    # step 6
    run_log("\n\n" + now() + "\n6. DC: calculating the resulting layer profile")

    for number_layers in [3]:
        configuration = "layerbins_equivol"
        number_ROIs = 100  # schaefer-100

        layer_profile_ROI_wise = LN2_layer_profiles_ROIwise(DC_dcw, label_file_cerebral, schaefer100_DC_functionalspace, number_ROIs, number_layers, configuration)  # generates a list of lists: layer_profile_ROI_wise[0..number_layers-1][0..numberROIs-1]
        if simulation_switch(-1) == 0:
            np.savetxt(output_path + "layer_profile_ROI_wise.txt", layer_profile_ROI_wise, delimiter=" ")
        else:
            layer_profile_ROI_wise = np.loadtxt(output_path + "layer_profile_ROI_wise.txt")

        # step 7: visualisation
        output_PNG = output_path + "LN2_PROFILES/DC_" + str(
            number_layers) + "_layers_" + configuration + "_schaefer100_s" + subject + ".png"
        brainspaces_schaefer100_stacked(layer_profile_ROI_wise, output_PNG)
        overview_files.append(output_PNG)

    copy_overview_files("DC_" + subject, overview_files)  # copy main files to separate folders

    print("DC done: " + subject)
    quit()


def dc_summary():
    simulation_switch(0)

    ##########
    overview_files = []

    scans = ['548rescan', '549rescan', '574rescan', '588rescan', '599rescan', '618']
    # output_path is automatically set to "summary_DC"

    labels = []
    for scan in scans:
        labels.append(datasets[scan]["label"])

    # step 1: read-in text files with individually calculated DC layer profiles and calculate average and error to create a plot
    result_list = []
    N = 0
    L = 7
    for label in labels:  # different label files
        file = "results/<subject>/output/DC_profile_full_7layers.txt"
        file = file.replace("<subject>", label)

        checkFiles([file])

        profile = np.loadtxt(file)  # list of 7 values
        result_list.append(profile)
        N += 1

    average = []
    error = []
    for l in range(L):
        values = []
        for n in range(N):
            values.append(result_list[n][l])

        avg, err = calculate_average_error(values)
        average.append(avg)
        error.append(err)

    all_y = [average]
    all_err = [error]
    list_linestyles = [None]
    list_labels = ["average", None]  # do not show label for standard deviation
    list_colors = [None]  # use standard color for the first, averaged graph
    for n in range(N):
        all_y.append(result_list[n])
        all_err.append(None)
        list_linestyles.append(":")
        list_labels.append(labels[n])

        mapping = [1, 2, 3, 4, 6, 8]  # avoid blue tones (reserved for averaged graph)
        color = col_seq[mapping[n]]
        list_colors.append(color)  # ("0.5")

    label_sorting = [1, 2, 2, 2, 2, 2, 2]  # label for average value plotted first, than all others in descending order

    cortical_depth = [1, 2, 3, 4, 5, 6, 7]
    dict_xticks = {}
    dict_xticks['sWM'] = min(cortical_depth)
    dict_xticks['pial'] = max(cortical_depth)

    # Plot the data
    title = "DC layer profile averaged over " + str(N) + " scans with 7 cortical layers"
    output_png = output_path + "DC_averaged.png"
    overview_files.append(output_png)

    plot_multiple_graphs(cortical_depth, all_y, output_png, list_y_errs=all_err, axis_title_x=None,
                         axis_title_y="Mean Signal", list_labels=list_labels,
                         plot_title=title, dict_xticks=dict_xticks, list_linestyles=list_linestyles,
                         list_colors=list_colors, label_sorting=label_sorting)

    # step 2: create mean image from set of DC maps in MNI space and calculate layer profile
    voxelsize_mm = 0.5  # everything will be upsampled or created with this voxel size (e.g., mean DC image, label file)

    list_DCW = []
    for label in labels:
        list_DCW.append("results/" + label + "/output/DC_DCW_MNIspace.nii.gz")

    list_DCW_upsampled = resampling(voxelsize_mm, list_DCW, method="Li")
    DC_mean = mean_image(list_DCW_upsampled, output_path + "DC_mean.nii.gz")
    overview_files.append(DC_mean)

    # create layers and calculate layer profile:

    maskWM = constant_MNItemplate_WM
    maskGM = constant_MNItemplate_GM
    maskCSF = constant_MNItemplate_CSF

    title = "average_DC_CPAC-dcw"
    label_file = create_LN2_labels(voxelsize_mm, title, maskWM, maskGM, maskCSF)  # the resulting label file is already upsampled and smoothed (maskGM etc. still have an isotropic voxel size of 1mm though)
    label_file = transform_space(label_file, DC_mean, sampling="NN")  # the label file was alreay created with target voxelsize and is, therefore, smoothed. Now make sure that space dimensions match DC_mean
    DC_anatomical_functionalspace = transform_space(constant_MNItemplate, DC_mean, sampling="Li")
    schaefer100_functionalspace = transform_space(constant_schaefer100, DC_mean, sampling="NN")

    # need upsampled GM of only the cerebral cortex: apply schaefer100_functionalspace as a mask to remove cerebellum and brainstem (careful: 1) schaefer100 is sometimes more than GM; 2) do not loose outer pixels due to slighly imperfect matches)
    run("fslmaths " + schaefer100_functionalspace + " -kernel box 1 -dilF " + temp_path + "schaefer100_functional_dil.nii.gz -odt int")  # dilate by 1mm
    cerebral_GM = apply_mask(label_file, temp_path + "schaefer100_functional_dil.nii.gz", suffix="cerebral-GM")  # use label file as GM source because it is already upsampled and smoothed
    label_file_cerebral = apply_mask(label_file, cerebral_GM, suffix="cerebral-GM")

    overview_files.append(DC_mean)
    overview_files.append(schaefer100_functionalspace)
    overview_files.append(DC_anatomical_functionalspace)
    overview_files.append(label_file_cerebral)

    # plot overlay maps after transforms and resampling
    title = "DC averaged from " + str(len(scans)) + " subjects"
    [plot_DC, plot_DC_GM] = dc_dcw_zscore_plots(DC_mean, DC_anatomical_functionalspace, cerebral_GM, output_path, title=title)
    overview_files.append(plot_DC)
    overview_files.append(plot_DC_GM)

    number_layers = 7
    configuration = "layerbins_equivol"
    layer_results = create_LN2_layers("DC-summary", label_file_cerebral, number_layers, DC_anatomical_functionalspace)  # mask is needed to only select the cerebral cortex (i.e., layer profile without brain stem and cerebellum)
    selected_layers = layer_results[configuration]
    overview_files.append(selected_layers)

    overview_files.append(mricrogl_plot_labelfile(DC_anatomical_functionalspace, label_file, orientation='x'))
    overview_files.append(mricrogl_plot_labelfile(DC_anatomical_functionalspace, label_file, orientation='y'))
    overview_files.append(mricrogl_plot_labelfile(DC_anatomical_functionalspace, label_file, orientation='z'))
    overview_files.append(mricrogl_plot_layerfile(DC_anatomical_functionalspace, selected_layers, number_layers))

    overview_files.append(mricrogl_plot_labelfile(DC_anatomical_functionalspace, label_file, orientation='x'))
    overview_files.append(mricrogl_plot_labelfile(DC_anatomical_functionalspace, label_file, orientation='y'))
    overview_files.append(mricrogl_plot_labelfile(DC_anatomical_functionalspace, label_file, orientation='z'))

    config_title = str(number_layers) + "_" + configuration
    output_PNG = output_path + "DC_" + str(number_layers) + "_layers_" + configuration + "_DC_average.png"
    layer_profile = create_LN2_profile(DC_mean, selected_layers, config_title, number_layers)
    layer_profile_PNG = layer_profile["PNG"]
    layer_profile_TXT = layer_profile["TXT"]
    layer_profile_table = layer_profile["table"]  # list of mean signals for different layers as defined in number_layers (starting from 0, i.e., closer to white matter)
    run("mv " + layer_profile_PNG + " " + output_PNG)
    np.savetxt(output_path + "DC_profile_full_" + str(number_layers) + "layers.txt", layer_profile_table, delimiter=" ")
    overview_files.append(output_PNG)

    copy_overview_files("DC_summary", overview_files)  # copy main files to separate folders

    ###############################################
    # create region-wise brainspace plots
    # uses "schaefer100_DC_functionalspace" that has already been created before the upsampling step

    # step 6
    run_log("\n\n" + now() + "\n6. DC: calculating the resulting layer profile")

    for number_layers in [3, 7]:
        configuration = "layerbins_equivol"
        number_ROIs = 100  # schaefer-100

        layer_profile_ROI_wise = LN2_layer_profiles_ROIwise(DC_mean, label_file_cerebral, schaefer100_functionalspace, number_ROIs, number_layers, configuration)  # generates a list of lists: layer_profile_ROI_wise[0..number_layers-1][0..numberROIs-1]
        if simulation_switch(-1) == 0:
            np.savetxt(output_path + "layer_profile_ROI_wise.txt", layer_profile_ROI_wise, delimiter=" ")
        else:
            layer_profile_ROI_wise = np.loadtxt(output_path + "layer_profile_ROI_wise.txt")

        # step 7: visualisation
        output_PNG = output_path + "LN2_PROFILES/DC_" + str(number_layers) + "_layers_" + configuration + "_schaefer100_s" + subject + ".png"
        overview_files.append(output_PNG)
        brainspaces_schaefer100_stacked(layer_profile_ROI_wise, output_PNG)

    copy_overview_files("DC_summary", overview_files)  # copy main files to separate folders
    quit("DC summary done.")



def reconAnatomical(subject, anatomicalT1):
    # runs recon-all on anatomical image

    checkFiles([anatomicalT1])

    jobs = []

    jobs.append("mri_convert <anatomical> <temp/>anatomical.mgh")
    jobs.append("mkdir ~/freesurfer")  # location for results of recon-all (avoid problems with access rights by locating this within the home directory)
    jobs.append("rm -rf -f ~/freesurfer/<myid>_anatomical")  # neccessary otherwise recon-all will stop if directory exists
    jobs.append("recon-all -i <temp/>anatomical.mgh -sd ~/freesurfer -subjid <myid>_anatomical -all")


    for job in jobs:
        job = job.replace("<anatomical>", anatomicalT1)
        job = job.replace("<myid>", subject)
        job = job.replace("<temp/>", temp_path)

        run(job)

    print ("recon-all on anatomical image finished for " + subject + ".")



def testing():
    # testing without changing subject data

    quit()












def layerfmri_task(subjectID, anatomical, fieldmap_mag, fieldmap_phase, functional, activity):
    # process task-related functional images: all steps from preprocessing to layer profile creation.
    checkFiles([anatomicalT1, fieldmap_mag, fieldmap_phase, functional])

    with open(global_runTXT, "a+") as file:
        file.write("\n\n\n" + subjectID + "\n")





    ################################################
    # check if all files and commands are available that are needed in the pipeline
    checkFiles([anatomicalT1, fieldmap_mag, fieldmap_phase, functional])
    checkCommands()
    print("All checks passed: Required input files and functions available.")
    ################################################
    simulation_switch(1)

    functional_current = functional

    run_log("\n\n" + now() + "\n1.   preprocessing")
    # 1. undistortion (i.e., correct warping from bias field)
    run_log("\n\n" + now() + "\n1.1 undistortion")
    output_undistortion = output_path + "sub-s" + subjectID + "_" + activity + "_undistorted.nii.gz" # simplify filename (activity selected from ["rft", "vistim"])
    [functional_current, biasfield] =  preprocessing_undistortion (functional_current, fieldmap_mag, fieldmap_phase, output_undistortion)

    # 2. bias field correction (Bias field correction is performed before motion correction because the bias field is independent of motion.)
    run_log("\n\n" + now() + "\n1.2 bias field correction")
    functional_current = preprocessing_biasfield4D(functional_current)


    # 3. motion correction to distortion corrrected functional images (with additional output of transformation matrix):
    run_log("\n\n" + now() + "\n1.3 motion correction")
    functional_current = preprocessing_motioncorrection(functional_current)





    # Preprocessing done.
    functional_preprocessed = functional_current
    preprocessing_insights(functional_preprocessed, functional, biasfield)








    # 4. create functional reference image (intention: Functional image includes a time series of 300 samples. create average of the first 10 to obtain one image and remove noise. Input must be motion corrected.)
    run_log("\n\n" + now() + "\n2.  reference image creation")
    [functional_first, functional_reference] = create_functional_reference (functional_preprocessed, functional)



    # 5. coregistrations
    run_log("\n\n" + now() + "\n3.   coregistrations")
    [anatomical_bet_functionalspace, anatomical_bet_functionalaligned, atlas_functionalspace, matrix_anatomical_to_functional] = coregistrations (functional_reference, anatomicalT1, constant_MNItemplate, constant_atlas, functional_first)
    # anatomical_bet_functionalspace:       pro: can be used to apply masks on funcitonal image
    #                                       con: some part will be missing if functional image is badly aligned (e.g., s618, s548re) and freeserufer segmentation would fail because brainmask.mgz cannnot be alligned properly
    #                                       => used to create a bet-mask for functional image that is more reliable than any bet operations on noisy functional image
    # anatomical_bet_functionalaligned      pro: aligned to functional image by rigid transform, reliable recon-all results
    #                                       con: cannot be directly used as a mask for functional image
    #                                       => used to create GM/WM/CSF masks via segnmentations from freesurfer by transforming recon-all's brainmask.mgz and aparc+aseg in talairach space to anatomical_bet_functionalaligned.
    #                                       In a later step, the LN2-lebel file that was created form said GM/WM/CSF masks, needs to be cut down to match functional space




    # 5. functional bet (You must make sure that input images have the same voxel size, for example, by upsampling them before!)
    run_log("\n\n" + now() + "\n5.   brain extraction in functional space with anatomical image as a reference")
    functional_bet = perform_bet_with_masks(functional_preprocessed, anatomical_bet_functionalspace) # Please Note: this must use anatomical_bet_functionalspace and not anatomical_bet_functionalaligned.
    create_preview(functional_preprocessed)
    create_preview(functional_bet)








    # 7. upsampling: two purposes: 1. better resolution for functional images, 2. ensure that anatomical_bet_functionalspace and functional have the same isotropic voxel size
    run_log("\n\n" + now() + "\n7.   upsampling")
    voxelsize_mm = 0.5

    # 7a) upsample the main images and maps
    # if "mask" or "atlas" is included in the filepath, nearest neighbor resampling will be selected instead of linear.
    [anatomical_bet_functionalspace_upsampled, anatomical_bet_functionalaligned_upsampled, functional_bet_upsampled] = resampling (voxelsize_mm, [anatomical_bet_functionalspace, anatomical_bet_functionalaligned, functional_bet])
    functional_bet_upsampled_preview = create_preview(functional_bet_upsampled)


    print("\n\nupsampling done!")



    # from now on, the above filesname are replaced by their upsampled version


    # 4. obtain intersection-free WM, GM, CSF, CTX masks (freesurfer_results  is folder, e.g., "/home/user/freesurfer/myid/")
    run_log("\n\n" + now() + "\n4.   recon-all for GM, WM and CSF masks from anatomical image in functional space")
    [maskWM, maskGM, maskCSF, cortex_mask_files] = WM_GM_CSF_CTX_masks(subjectID, anatomical_bet_functionalaligned)  # uses recon-all
    # Please Note: this must use anatomical_bet_functionalaligned and not anatomical_bet_functionalspace.


    # compute activity map (for 0.5mm isotropic voxel size: 60 +/- 2 hours and approx. 48 GB RAM)
    activity_map = output_path + "activity_map.nii.gz" # defined in "feat_processing", uses "stats/cope1.nii.gz"
    if False: #not checkFiles([activity_map], True):
        print("\n\nstarting FEAT\n\n")
        feat_processing(functional_bet)
        print("\n\n FEAT DONE: " + subjectID  +"\n\n")
        quit()
    else:
        print ("\n Skipped producing a new activity map. File already found form previous run: " + activity_map)



    # resample to finest possible
    #voxelsize_mm = 0.3
    #[activity_map_upsampled, anatomical_bet_functionalspace_upsampled, functional_bet_upsampled_preview] = resampling (voxelsize_mm, [activity_map, anatomical_bet_functionalspace, functional_bet_upsampled_preview])
    voxelsize_mm = 0.5

    # 8. calculate label file for LN2layers from GM, WM and CSF masks
    run_log("\n\n8..   creation of a label file as an input for laynii")
    LN2_input_labels = create_LN2_labels (voxelsize_mm, "main", maskWM, maskGM, maskCSF)
    print("main pipeline: label file computed")


    simulation_switch(0)
    print("\n\nuntil here ok.\n\”")


    # 9. generate layers
    run_log("\n\n9.   computation of multiple layers with laynii")
    number_layers = 7
    configuration = "layerbins_equivol"
    layer_results = create_LN2_layers("layerfmri", LN2_input_labels, number_layers, anatomical_bet_functionalaligned)
    selected_layers = layer_results[configuration]
    print("main pipeline: layers generated")




    # 10: calculate layer profile
    run_log("\n\n" + now() + "\n10. DC: calculating the resulting layer profile")
    config_title = str(number_layers) + "_" + configuration
    layer_profile = create_LN2_profile(activity_map, transform_space(selected_layers, activity_map, sampling="Li", suffix="main"), config_title, number_layers, short_name="full_activity_map")
    # transform_space(selected_layers, activity_map) above is needed to go from functional_aligned to functional_spcace. The difference is in the image dimensions and is essential for fslmaths -mas to perform correctly.
    layer_profile_PNG = layer_profile["PNG"]
    layer_profile_TXT = layer_profile["TXT"]
    layer_profile_table = layer_profile["table"]  # list of mean signals for different layers as defined in nnumber_layers (starting from 0, i.e., closer to white matter)


    np.savetxt(output_path + "MAINPIPELINE_profile_full_" + str(number_layers) + "layers.txt", layer_profile_table, delimiter=" ")



    simulation_switch(0)



    for num in [7, 15]:
        if not num == number_layers: # otherwise: recycle previous layer file
            number_layers = num
            layer_results = create_LN2_layers("layerfmri", LN2_input_labels, number_layers, anatomical_bet_functionalaligned)

        # create profiles for different ROIs
        configuration = "layerbins_equivol"
        config_title = str(number_layers) + "_" + configuration


        ROIs = {}
        ROIs["visual-rhlh"]             = ["ctx-rh-lingual", "ctx-rh-cuneus", "ctx-rh-pericalcarine", "ctx-rh-lateraloccipital", "ctx-lh-lingual", "ctx-lh-cuneus", "ctx-lh-pericalcarine", "ctx-lh-lateraloccipital"]
        ROIs["motor-lh"]                = ["ctx-lh-precentral"]
        ROIs["motor-rh"]                = ["ctx-rh-precentral"]
        ROIs["somatosensory-lh"]        = ["ctx-lh-postcentral"]


        ROI_label_folder = output_path + "ROI_labels/"
        run ("mkdir " + ROI_label_folder)


        ROI_tables = {}
        for ROI, cortex_list in ROIs.items():
            ROI_mask = ROI_label_folder + "mask_" + ROI + ".nii.gz"
            ROI_layer_results = ROI_label_folder + "layers_" + ROI + ".nii.gz"
            job = ""
            for n, cortex in enumerate(cortex_list):
                if n == 0:
                    job = "fslmaths " + cortex_mask_files[cortex]
                else:
                    job += " -max " + cortex_mask_files[cortex]
            job += " " + ROI_mask
            run(job)
            [ROI_mask] = resampling (voxelsize_mm, [ROI_mask]) # same voxel size as anatomical_functional and layer_results (currently 0.5 mm)
            print("resampling done: " + ROI_mask)

            run("fslmaths " + layer_results[configuration] + " -mas " + ROI_mask + " " + ROI_layer_results)
            layer_profile = create_LN2_profile(activity_map, transform_space(ROI_layer_results, activity_map, sampling="NN", suffix=ROI), config_title, number_layers, short_name=ROI)
            # transform_space(ROI_layer_results, activity_map) above is needed to go from functional_aligned to functional_spcace. The difference is in the image dimensions and is essential for fslmaths -mas to perform correctly.
            layer_profile_table = layer_profile["table"]  # list of mean signals for different layers as defined in nnumber_layers
            ROI_tables[ROI] = layer_profile_table
            np.savetxt(ROI_label_folder + ROI  + "_" + str(number_layers) + "layers.txt", layer_profile_table, delimiter=" ")
            print ("  calculated profile: " + ROI)


        # Plot the data

        cortical_depth = [i for i in range(1, number_layers+1)] # e.g., for number_layers=7: cortical_depth = [1, 2, 3, 4, 5, 6, 7]
        title = subjectID + ": " + str(number_layers) + " layers - " + os.path.basename(functional)
        output_png = output_path + "LN2_PROFILES/comparison_" + subjectID + "_" + str(number_layers) + "_layers.png"


        list_yvals = []
        list_labels = []
        for ROI, meanSignal in ROI_tables.items():
            list_yvals.append(meanSignal) # meanSignal is a list of values that corresponeds to [1 ... N] cortical layers in "cortical_depth"
            list_labels.append(ROI)

        dict_xticks = {}
        dict_xticks['sWM'] = min(cortical_depth)
        dict_xticks['pial'] = max(cortical_depth)


        plot_multiple_graphs(cortical_depth, list_yvals, output_png,
                             axis_title_x       = None,
                             axis_title_y       = "Mean Signal",
                             list_labels        = list_labels,
                             plot_title         = "layer profile for " + title,
                             dict_xticks        = dict_xticks)





    print("Main pipeline: layer profile calculated for " + subjectID)
    quit()











