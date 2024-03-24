datasets = {}



anatomical_path             =   "data/<label>/anat/" # e.g., label = "sub-s588"
fieldmap_path               =   "data/<label>/fmap/"
functional_path             =   "data/<label>/func/"
DC_path                     =   "data/<label>/DC/"

atlas_path = "atlas/"
constant_MNItemplate        = atlas_path + "MNI152_2009_template.nii.gz"
constant_MNItemplate_GM     = atlas_path + "MNI152_2009_GM.nii.gz"
constant_MNItemplate_WM     = atlas_path + "MNI152_2009_WM.nii.gz"
constant_MNItemplate_CSF    = atlas_path + "MNI152_2009_CSF.nii.gz"

constant_atlas              = atlas_path + "HCPMMP1_on_MNI152_ICBM2009a_nlin_hd.nii.gz" # glasser atlas
constant_schaefer100        = atlas_path + "schaefer100MNI.nii.gz"




################################################################################################
label = "s548-7T"
content = {}
content["label"] = label

content["anatomicalT1"]        = anatomical_path           + "sub-s548_acq-mp2rage_T1w.nii.gz"
content["anatomicalT2"]        = anatomical_path           + "sub-s548_T2w.nii.gz"
content["fieldmap_mag"]        = fieldmap_path             + "sub-s548_acq-0p8mm_magnitude2.nii.gz"
content["fieldmap_phase"]      = fieldmap_path             + "sub-s548_acq-0p8mm_phasediff.nii.gz"
content["functionalRest"]      = functional_path           + "sub-s548_task-rest_acq-0p8mm_bold.nii.gz"
content["functionalTap"]       = "N/A"
content["functionalChecker"]   = "N/A"
content["DC_anatomical"]       = DC_path                   + "sub-s548_ses-1_desc-head_T1w.nii.gz"
content["DC_mean_bold"]        = DC_path                   + "sub-s548_ses-1_task-rest_acq-0p8mm_desc-mean_bold.nii.gz"
content["DC_mean_bold_GM_mask"]= DC_path                   + "sub-s548_ses-1_space-bold_label-GM_mask.nii.gz"
content["DC_dcw"]              = DC_path                   + "sub-s548_ses-1_task-rest_acq-0p8mm_desc-sm5_dcw.nii.gz"
content["DC_GM_mask"]          = DC_path                   + "sub-s548_ses-1_label-GM_mask.nii.gz"
content["DC_WM_mask"]          = DC_path                   + "sub-s548_ses-1_label-WM_mask.nii.gz"
content["DC_CSF_mask"]         = DC_path                   + "sub-s548_ses-1_label-CSF_mask.nii.gz"

datasets["548"] = content


################################################################################################
label = "s574-7T"
content = {}
content["label"] = label

content["anatomicalT1"]        = anatomical_path           + "8-t1_mp2rage_sag_HCP_0pp8mm3_T1map_BL5_230717_UNI_DEN.nii"
content["anatomicalT2"]        = anatomical_path           + "sub-s574_acq-grappatest_T2w.nii.gz"
content["fieldmap_mag"]        = fieldmap_path             + "sub-s574_acq-0p8mm_magnitude2.nii.gz"
content["fieldmap_phase"]      = fieldmap_path             + "sub-s574_acq-0p8mm_phasediff.nii.gz"
content["functionalRest"]      = functional_path           + "sub-s574_task-rest_acq-0p8mm_bold.nii.gz"
content["functionalTap"]       = "N/A"
content["functionalChecker"]   = "N/A"
content["DC_anatomical"]       = DC_path                   + "sub-s574_ses-1_desc-head_T1w.nii.gz"
content["DC_mean_bold"]        = DC_path                   + "sub-s574_ses-1_task-rest_acq-1p3mm_desc-mean_bold.nii.gz"
content["DC_mean_bold_GM_mask"]= DC_path                   + "sub-s574_ses-1_space-bold_label-GM_mask.nii.gz"
content["DC_dcw"]              = DC_path                   + "sub-s574_ses-1_task-rest_acq-1p3mm_desc-sm5_dcw.nii.gz"
content["DC_GM_mask"]          = DC_path                   + "sub-s574_ses-1_label-GM_mask.nii.gz"
content["DC_WM_mask"]          = DC_path                   + "sub-s574_ses-1_label-WM_mask.nii.gz"
content["DC_CSF_mask"]         = DC_path                   + "sub-s574_ses-1_label-CSF_mask.nii.gz"

datasets["574"] = content

################################################################################################
label = "s574-7Trescan"
content = {}
content["label"] = label

content["anatomicalT1"]        = anatomical_path                + "sub-s574_acq-mp2rage_T1w.nii.gz"
content["anatomicalT2"]        = anatomical_path                + "sub-s574_T2w.nii.gz"
content["fieldmap_mag"]        = fieldmap_path                  + "sub-s574_acq-0p8mm_magnitude2.nii.gz"
content["fieldmap_phase"]      = fieldmap_path                  + "sub-s574_acq-0p8mm_phasediff.nii.gz"
content["functionalRest"]      = functional_path                + "sub-s574_task-rest_acq-0p8mm_bold.nii.gz"
content["functionalTap"]       = functional_path                + "sub-s574_task-rft_acq-0p8mm_bold.nii.gz"
content["functionalChecker"]   = "N/A"
content["DC_anatomical"]       = DC_path                        + "sub-s574_ses-1_desc-head_T1w.nii.gz"
content["DC_mean_bold"]        = DC_path                        + "sub-s574_ses-1_task-rest_acq-0p8mm_desc-mean_bold.nii.gz"
content["DC_mean_bold_GM_mask"]= DC_path                        + "sub-s574_ses-1_space-bold_label-GM_mask.nii.gz"
content["DC_dcw"]              = DC_path                        + "sub-s574_ses-1_task-rest_acq-0p8mm_desc-sm5_dcw.nii.gz"
content["DC_GM_mask"]          = DC_path                        + "sub-s574_ses-1_label-GM_mask.nii.gz"
content["DC_WM_mask"]          = DC_path                        + "sub-s574_ses-1_label-WM_mask.nii.gz"
content["DC_CSF_mask"]         = DC_path                        + "sub-s574_ses-1_label-CSF_mask.nii.gz"

datasets["574rescan"] = content



################################################################################################
label = "s588-7T"
content = {}
content["label"] = label

content["anatomicalT1"]        = anatomical_path           + "sub-s588_acq-mprage_T1w.nii.gz"
content["anatomicalT2"]        = anatomical_path           + "sub-s588_T2w.nii.gz"
content["fieldmap_mag"]        = fieldmap_path             + "sub-s588_acq-0p8mm_magnitude2.nii.gz"
content["fieldmap_phase"]      = fieldmap_path             + "sub-s588_acq-0p8mm_phasediff.nii.gz"
content["functionalRest"]      = functional_path           + "sub-s588_task-rest_acq-0p8mm_bold.nii.gz"
content["functionalTap"]       = "N/A"
content["functionalChecker"]   = "N/A"
content["DC_anatomical"]       = DC_path                   + "sub-s588_ses-1_desc-head_T1w.nii.gz"
content["DC_mean_bold"]        = DC_path                   + "sub-s588_ses-1_task-rest_acq-0p8mm_desc-mean_bold.nii.gz"
content["DC_mean_bold_GM_mask"]= DC_path                   + "sub-s588_ses-1_space-bold_label-GM_mask.nii.gz"
content["DC_dcw"]              = DC_path                   + "sub-s588_ses-1_task-rest_acq-0p8mm_desc-sm5_dcw.nii.gz"
content["DC_GM_mask"]          = DC_path                   + "sub-s588_ses-1_label-GM_mask.nii.gz"
content["DC_WM_mask"]          = DC_path                   + "sub-s588_ses-1_label-WM_mask.nii.gz"
content["DC_CSF_mask"]         = DC_path                   + "sub-s588_ses-1_label-CSF_mask.nii.gz"

datasets["588"] = content


################################################################################################
label = "s588-7Trescan"
content = {}
content["label"] = label

content["anatomicalT1"]        = anatomical_path           + "sub-s588_acq-mprage_T1w.nii.gz"
content["anatomicalT2"]        = anatomical_path           + "sub-s588_T2w.nii.gz"
content["fieldmap_mag"]        = fieldmap_path             + "sub-s588_acq-0p8mm_magnitude2.nii.gz"
content["fieldmap_phase"]      = fieldmap_path             + "sub-s588_acq-0p8mm_phasediff.nii.gz"
content["functionalRest"]      = functional_path           + "sub-s588_task-rest_acq-0p8mm_bold.nii.gz"
content["functionalTap"]       = "N/A"
content["functionalChecker"]   = "N/A"
content["DC_anatomical"]       = DC_path                   + "sub-s588_ses-1_desc-head_T1w.nii.gz"
content["DC_mean_bold"]        = DC_path                   + "sub-s588_ses-1_task-rest_acq-0p8mm_desc-mean_bold.nii.gz"
content["DC_mean_bold_GM_mask"]= DC_path                   + "sub-s588_ses-1_space-bold_label-GM_mask.nii.gz"
content["DC_dcw"]              = DC_path                   + "sub-s588_ses-1_task-rest_acq-0p8mm_desc-sm5_dcw.nii.gz"
content["DC_GM_mask"]          = DC_path                   + "sub-s588_ses-1_label-GM_mask.nii.gz"
content["DC_WM_mask"]          = DC_path                   + "sub-s588_ses-1_label-WM_mask.nii.gz"
content["DC_CSF_mask"]         = DC_path                   + "sub-s588_ses-1_label-CSF_mask.nii.gz"

datasets["588rescan"] = content




################################################################################################
label = "s597-7T"
content = {}
content["label"] = label

content["anatomicalT1"]        = anatomical_path               + "sub-s597_acq-mp2rage_run-02_T1w.nii.gz"
content["anatomicalT2"]        = anatomical_path               + "sub-s597_T2w.nii.gz"
content["fieldmap_mag"]        = fieldmap_path                 + "sub-s597_acq-0p8mm_magnitude2.nii.gz"
content["fieldmap_phase"]      = fieldmap_path                 + "sub-s597_acq-0p8mm_phasediff.nii.gz"
content["functionalRest"]      = functional_path               + "sub-s597_task-rest_acq-0p8mm_bold.nii.gz" # only task-rest available
content["functionalTap"]       = "N/A"
content["functionalChecker"]   = "N/A"
content["DC_anatomical"]       = "N/A"
content["DC_mean_bold"]        = "N/A"
content["DC_dcw"]              = "N/A"
content["DC_GM_mask"]          = "N/A"
content["DC_WM_mask"]          = "N/A"
content["DC_CSF_mask"]         = "N/A"

datasets["597"] = content



################################################################################################
label = "s599-7T"
content = {}
content["label"] = label

content["anatomicalT1"]        = anatomical_path           + "sub-s599_acq-mprage_T1w.nii.gz"
content["anatomicalT2"]        = anatomical_path           + "sub-s599_T2w.nii.gz"
content["fieldmap_mag"]        = fieldmap_path             + "sub-s599_acq-0p8mm_magnitude2.nii.gz"
content["fieldmap_phase"]      = fieldmap_path             + "sub-s599_acq-0p8mm_phasediff.nii.gz"
content["functionalRest"]      = functional_path           + "sub-s599_task-rest_acq-0p8mm_bold.nii.gz"
content["functionalTap"]       = "N/A"
content["functionalChecker"]   = "N/A"
content["DC_anatomical"]       = DC_path                   + "sub-s599_ses-1_desc-head_T1w.nii.gz"
content["DC_mean_bold"]        = DC_path                   + "sub-s599_ses-1_task-rest_acq-0p8mm_desc-mean_bold.nii.gz"
content["DC_mean_bold_GM_mask"]= DC_path                   + "sub-s599_ses-1_space-bold_label-GM_mask.nii.gz"
content["DC_dcw"]              = DC_path                   + "sub-s599_ses-1_task-rest_acq-0p8mm_desc-sm5_dcw.nii.gz"
content["DC_GM_mask"]          = DC_path                   + "sub-s599_ses-1_label-GM_mask.nii.gz"
content["DC_WM_mask"]          = DC_path                   + "sub-s599_ses-1_label-WM_mask.nii.gz"
content["DC_CSF_mask"]         = DC_path                   + "sub-s599_ses-1_label-CSF_mask.nii.gz"

datasets["599"] = content


################################################################################################
label = "s618-7T"
content = {}
content["label"] = label

content["anatomicalT1"]        = anatomical_path           + "sub-s618_acq-mp2rage_T1w.nii.gz"
content["anatomicalT2"]        = anatomical_path           + "sub-s618_T2w.nii.gz"
content["fieldmap_mag"]        = fieldmap_path             + "sub-s618_acq-0p8mm_magnitude2.nii.gz"
content["fieldmap_phase"]      = fieldmap_path             + "sub-s618_acq-0p8mm_phasediff.nii.gz"
content["functionalRest"]      = functional_path           + "sub-s618_task-rest_acq-0p8mm_bold.nii.gz"
content["functionalTap"]       = functional_path           + "sub-s618_task-rft_acq-0p8mm_bold.nii.gz"
content["functionalChecker"]   = "N/A"
content["DC_anatomical"]       = DC_path                   + "sub-s618_ses-1_desc-head_T1w.nii.gz"
content["DC_mean_bold"]        = DC_path                   + "sub-s618_ses-1_task-rest_acq-0p8mm_desc-mean_bold.nii.gz"
content["DC_mean_bold_GM_mask"]= DC_path                   + "sub-s618_ses-1_space-bold_label-GM_mask.nii.gz"
content["DC_dcw"]              = DC_path                   + "sub-s618_ses-1_task-rest_acq-0p8mm_desc-sm5_dcw.nii.gz"
content["DC_GM_mask"]          = DC_path                   + "sub-s618_ses-1_label-GM_mask.nii.gz"
content["DC_WM_mask"]          = DC_path                   + "sub-s618_ses-1_label-WM_mask.nii.gz"
content["DC_CSF_mask"]         = DC_path                   + "sub-s618_ses-1_label-CSF_mask.nii.gz"

datasets["618"] = content




################################################################################################
label = "s548-7Trescan"
content = {}
content["label"] = label

content["anatomicalT1"]        = anatomical_path           + "sub-s548_acq-mp2rage_T1w.nii.gz"
content["anatomicalT2"]        = anatomical_path           + "sub-s548_T2w.nii.gz"
content["fieldmap_mag"]        = fieldmap_path             + "sub-s548_acq-0p8mm_magnitude2.nii.gz"
content["fieldmap_phase"]      = fieldmap_path             + "sub-s548_acq-0p8mm_phasediff.nii.gz"
content["functionalRest"]      = functional_path           + "sub-s548_task-rest_acq-0p8mm_bold.nii.gz"
content["functionalTap"]       = "N/A"
content["functionalChecker"]   = functional_path           + "sub-s548_task-vistim_acq-0p8mm_bold.nii.gz"
content["DC_anatomical"]       = DC_path                   + "sub-s548_ses-1_desc-head_T1w.nii.gz"
content["DC_mean_bold"]        = DC_path                   + "sub-s548_ses-1_task-rest_acq-0p8mm_desc-mean_bold.nii.gz"
content["DC_mean_bold_GM_mask"]= DC_path                   + "sub-s548_ses-1_space-bold_label-GM_mask.nii.gz"
content["DC_dcw"]              = DC_path                   + "sub-s548_ses-1_task-rest_acq-0p8mm_desc-sm5_dcw.nii.gz"
content["DC_GM_mask"]          = DC_path                   + "sub-s548_ses-1_label-GM_mask.nii.gz"
content["DC_WM_mask"]          = DC_path                   + "sub-s548_ses-1_label-WM_mask.nii.gz"
content["DC_CSF_mask"]         = DC_path                   + "sub-s548_ses-1_label-CSF_mask.nii.gz"

datasets["548rescan"] = content



################################################################################################
label = "s549-7Trescan"
content = {}
content["label"] = label

content["anatomicalT1"]        = anatomical_path           + "sub-s549_acq-mp2rage_T1w.nii.gz"
content["anatomicalT2"]        = anatomical_path           + "sub-s549_T2w.nii.gz"
content["fieldmap_mag"]        = fieldmap_path             + "sub-s549_acq-0p8mm_magnitude2.nii.gz"
content["fieldmap_phase"]      = fieldmap_path             + "sub-s549_acq-0p8mm_phasediff.nii.gz"
content["functionalRest"]      = functional_path           + "sub-s549_task-rest_acq-0p8mm_bold.nii.gz"
content["functionalTap"]       = functional_path           + "sub-s549_task-rft_acq-0p8mm_bold.nii.gz"
content["functionalChecker"]   = "N/A"
content["DC_anatomical"]       = DC_path                   + "sub-s549_ses-1_desc-head_T1w.nii.gz"
content["DC_mean_bold"]        = DC_path                   + "sub-s549_ses-1_task-rest_acq-0p8mm_desc-mean_bold.nii.gz"
content["DC_mean_bold_GM_mask"]= DC_path                   + "sub-s549_ses-1_space-bold_label-GM_mask.nii.gz"
content["DC_dcw"]              = DC_path                   + "sub-s549_ses-1_task-rest_acq-0p8mm_desc-sm5_dcw.nii.gz"
content["DC_GM_mask"]          = DC_path                   + "sub-s549_ses-1_label-GM_mask.nii.gz"
content["DC_WM_mask"]          = DC_path                   + "sub-s549_ses-1_label-WM_mask.nii.gz"
content["DC_CSF_mask"]         = DC_path                   + "sub-s549_ses-1_label-CSF_mask.nii.gz"

datasets["549rescan"] = content


################################################################################################
label = "s599-7Trescan"
content = {}
content["label"] = label

content["anatomicalT1"]        = anatomical_path           + "sub-s599_acq-mp2rage_T1w.nii.gz"
content["anatomicalT2"]        = anatomical_path           + "sub-s599_T2w.nii.gz"
content["fieldmap_mag"]        = fieldmap_path             + "sub-s599_acq-0p8mm_magnitude2.nii.gz"
content["fieldmap_phase"]      = fieldmap_path             + "sub-s599_acq-0p8mm_phasediff.nii.gz"
content["functionalRest"]      = functional_path           + "sub-s599_task-rest_acq-0p8mm_bold.nii.gz"
content["functionalTap"]       = "N/A"
content["functionalChecker"]   = functional_path           + "sub-s599_task-vistim_acq-0p8mm_bold.nii.gz"
content["DC_anatomical"]       = DC_path                   + "sub-s599_ses-1_desc-head_T1w.nii.gz"
content["DC_mean_bold"]        = DC_path                   + "sub-s599_ses-1_task-rest_acq-0p8mm_desc-mean_bold.nii.gz"
content["DC_mean_bold_GM_mask"]= DC_path                   + "sub-s599_ses-1_space-bold_label-GM_mask.nii.gz"
content["DC_dcw"]              = DC_path                   + "sub-s599_ses-1_task-rest_acq-0p8mm_desc-sm5_dcw.nii.gz"
content["DC_GM_mask"]          = DC_path                   + "sub-s599_ses-1_label-GM_mask.nii.gz"
content["DC_WM_mask"]          = DC_path                   + "sub-s599_ses-1_label-WM_mask.nii.gz"
content["DC_CSF_mask"]         = DC_path                   + "sub-s599_ses-1_label-CSF_mask.nii.gz"

datasets["599rescan"] = content


################################################################################################
label = "s640-7T"
content = {}
content["label"] = label

content["anatomicalT1"]        = anatomical_path           + "sub-s640_acq-mp2rage_T1w.nii.gz"
content["anatomicalT2"]        = "N/A"
content["fieldmap_mag"]        = fieldmap_path             + "sub-s640_acq-0p8mm_magnitude2.nii.gz"
content["fieldmap_phase"]      = fieldmap_path             + "sub-s640_acq-0p8mm_phasediff.nii.gz"
content["functionalRest"]      = functional_path           + "sub-s640_task-rest_acq-0p8mm_bold.nii.gz"
content["functionalTap"]       = "N/A"
content["functionalChecker"]   = functional_path           + "sub-s640_task-vistim_acq-0p8mm_bold.nii.gz"
content["DC_anatomical"]       = "N/A"
content["DC_mean_bold"]        = "N/A"
content["DC_mean_bold_GM_mask"]= "N/A"
content["DC_dcw"]              = "N/A"
content["DC_GM_mask"]          = "N/A"
content["DC_WM_mask"]          = "N/A"
content["DC_CSF_mask"]         = "N/A"

datasets["640"] = content


################################################################################################
label = "s641-7T"
content = {}
content["label"] = label

content["anatomicalT1"]        = anatomical_path           + "sub-s641_acq-mp2rage_T1w.nii.gz"
content["anatomicalT2"]        = "N/A"
content["fieldmap_mag"]        = fieldmap_path             + "sub-s641_acq-0p8mm_magnitude2.nii.gz"
content["fieldmap_phase"]      = fieldmap_path             + "sub-s641_acq-0p8mm_phasediff.nii.gz"
content["functionalRest"]      = functional_path           + "sub-s641_task-rest_acq-0p8mm_bold.nii.gz"
content["functionalTap"]       = "N/A"
content["functionalChecker"]   = functional_path           + "sub-s641_task-vistim_acq-0p8mm_run-02_bold.nii.gz"
content["DC_anatomical"]       = "N/A"
content["DC_mean_bold"]        = "N/A"
content["DC_mean_bold_GM_mask"]= "N/A"
content["DC_dcw"]              = "N/A"
content["DC_GM_mask"]          = "N/A"
content["DC_WM_mask"]          = "N/A"
content["DC_CSF_mask"]         = "N/A"

datasets["641"] = content




#################################
# replace foldernames:
for quicklabel, content in datasets.items():
    label = content["label"]
    for entry, filename in content.items():
        filename = filename.replace("<label>", label)
        content[entry] = filename
    datasets[quicklabel] = content











def display_available_jobs():
    # checks dataset and, based on available files, determines, which jobs can be performed
    list_tapping = ["anatomicalT1", "fieldmap_mag", "fieldmap_phase", "functionalRest", "functionalTap"]
    list_visual = ["anatomicalT1", "fieldmap_mag", "fieldmap_phase", "functionalRest", "functionalChecker"]
    list_DC = ["DC_anatomical", "DC_mean_bold", "DC_mean_bold_GM_mask", "DC_dcw", "DC_WM_mask", "DC_CSF_mask"]

    lists={}
    lists["tapping"]    = list_tapping
    lists["visual"]     = list_visual
    lists["DC"]         = list_DC

    commands = {}
    commands["tapping"] = "cd ~/Desktop/pipeline/\npython3 main.py layerfmri <id>"
    commands["visual"] = "cd ~/Desktop/pipeline/\npython3 main.py layerfmri <id>"
    commands["DC"] = "cd ~/Desktop/pipeline/\npython3 main.py dc <id>"


    results = ""
    results_commands = ""

    for analysis, files in lists.items():
        for subject, content in datasets.items():
            possible = True
            for file in files:
                if not file in content:
                    possible = False
                else:
                    if content[file] == "N/A":
                        possible = False

            if possible:
                results += analysis + ": " + subject + "\n"

                command = "\n" + commands[analysis]
                command = command.replace("<id>", subject)
                results_commands += command + "\n"

    print(results)
    print(results_commands)











