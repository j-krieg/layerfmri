# prepare output directories:
#   global_subject_folder
#   output_path
#   temp_path
#   mask_path
#   layniiresult_path
#   atlas_path
#   constant_MNItemplate
#   constant_atlas



global_subject_folder = "results/" + label + "/"
global_runTXT = global_subject_folder + "run.txt" # log commands and runtimes here
output_path = global_subject_folder + "output/"
temp_path = global_subject_folder + "temp/"
mask_path = output_path + "masks/"
layniiresult_path = global_subject_folder + "output/layer_results/"

atlas_path = "atlas/"
constant_MNItemplate = atlas_path + "MNI152_2009_template.nii.gz"
constant_atlas = atlas_path + "HCPMMP1_on_MNI152_ICBM2009a_nlin_hd.nii.gz"  # Glasser


if not os.path.exists("results"):
    os.makedirs("results")

if not os.path.exists(global_subject_folder):
    os.makedirs(global_subject_folder)

if not os.path.exists(output_path):
    os.makedirs(output_path)

if not os.path.exists(temp_path):
    os.makedirs(temp_path)

if not os.path.exists(mask_path):
    os.makedirs(mask_path)

if not os.path.exists(layniiresult_path):
    os.makedirs(layniiresult_path)








