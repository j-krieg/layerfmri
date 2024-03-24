


def stack_images (list_files, output_file, direction="horizontal", overlap=0.0):
    # overlap: percentage of single image dimension in stacking directoin
    N = len(list_files)

    width_new = -1
    height_new = -1
    for file in list_files:
        input = Image.open(file)

        width, height = input.size

        if width_new < width:
            width_new = width

        if height_new < height:
            height_new = height

    if direction == "horizontal":
        width_new = int(width_new*(1.0-overlap))
        width_single = width_new
        width_new *= N
    else:
        height_new = int(height_new*(1.0-overlap))
        height_single = height_new
        height_new *= N




    final_image = Image.new("RGB", (width_new, height_new), "black")

    for n, file in enumerate(list_files):
        input = Image.open(file)
        if direction == "horizontal":
            final_image.paste(input, (n * width_single, 0))
        else:
            final_image.paste(input, (0, n * height_single))



    final_image.save(output_file)
    return output_file




def alternative_color(color_hex):

    r, g, b = int(color_hex[1:3], 16), int(color_hex[3:5], 16), int(color_hex[5:7], 16)
    h, l, s = colorsys.rgb_to_hls(r / 255, g / 255, b / 255)

    l = 1.0
    s = (0.6 * s) + 0.4

    r, g, b = colorsys.hls_to_rgb(h, l, s)

    return color_hex # TO-DO



def calculate_image_margins (input_file):  # finds the minimum dimensions of the image such that the boarders do not contain empty voxels (within a small pre-defined margin)
    # returns [[slices_x_low, slices_x_high], [slices_y_low, slices_z_high], [slices_z_low, slices_z_high]]

    # determine pixdim (mm/slice) of original image
    info = fsl_info(input_file)
    x_pixdim = info["x_pixdim"]  # mm
    y_pixdim = info["y_pixdim"]
    z_pixdim = info["z_pixdim"]

    x_dim = info["x_dim"]  # slices
    y_dim = info["y_dim"]
    z_dim = info["z_dim"]
    original_pixdim = [x_pixdim, y_pixdim, z_pixdim]

    intensity_threshold = 1000  # minimum intensity
    minimum_cluster_size = 5  # shrink image by given number of voxels to remove noise etc. (in downsampled space)
    margin_voxels = 3  # add small margin (refers to pixdim in original image)
    downsampling_mm = 2  # speed-up search for empty voxels

    [image_small] = resampling(downsampling_mm, [input_file], method="NN", folder=temp_path)  # speed up
    image_filtered = temp_path + "img_filtered.nii.gz"

    jobs = []
    jobs.append("fslmaths <image> -thr <intensity_threshold> <image_filtered> -odt int")
    jobs.append("fslmaths <image_filtered> -abs -bin <image_filtered> -odt int")
    jobs.append(
        "fslmaths <image_filtered> -kernel box <minimum_cluster_size> -ero <image_filtered> -odt int")  # remove small clsuters. Amount of voxels will be added to the result later.

    for job in jobs:
        job = job.replace("<image>", image_small)
        job = job.replace("<image_filtered>", image_filtered)
        job = job.replace("<intensity_threshold>", str(intensity_threshold))
        job = job.replace("<minimum_cluster_size>", str(minimum_cluster_size))
        job = job.replace("<margin_voxels>", str(margin_voxels))
        run(job)

    img = nib.load(image_filtered)
    data = img.get_fdata()

    profile_x = [np.mean(data[i, :, :]) for i in range(data.shape[0])]
    profile_y = [np.mean(data[:, i, :]) for i in range(data.shape[1])]
    profile_z = [np.mean(data[:, i, :]) for i in range(data.shape[2])]

    empty_slices = [[0, 0], [0, 0], [0,
                                     0]]  # [[x_leading, x_trailing], [y_leading, y_trailing], [z_leading, z_trailing] (in original space)

    for dim, profile in enumerate([profile_x, profile_y, profile_z]):
        empty_slices_leading = next((i for i, x in enumerate(profile) if x != 0.0), len(profile))
        empty_slices_trailing = next((i for i, x in enumerate(profile[::-1]) if x != 0.0), len(profile))

        # determine the number of slices to be deleted in original space after consideration of the small margin
        empty_slices[dim][0] = max(
            int((empty_slices_leading - minimum_cluster_size) * downsampling_mm / original_pixdim[dim]) - margin_voxels,
            0)  # number of slices to be removed in original image
        empty_slices[dim][1] = max(
            int((empty_slices_trailing - minimum_cluster_size) * downsampling_mm / original_pixdim[
                dim]) - margin_voxels, 0)

    return empty_slices


def strip_image_in_y_direction(image, activity):
    # removes empty spaces from anatomical image
    # simulaneously transforms the activity map such that both still match afterwards

    output_image = temp_path + (os.path.basename(image)).replace(".nii.gz", "_ystripped.nii.gz")
    output_activity = temp_path + (os.path.basename(activity)).replace(".nii.gz", "_ystripped.nii.gz")

    empty_slices = calculate_image_margins(
        image)  # returns [[slices_x_low, slices_x_high], [slices_y_low, slices_z_high], [slices_z_low, slices_z_high]]

    y_low = empty_slices[1][0]
    y_high = empty_slices[1][1]

    info = fsl_info(image)
    y_pixdim = info["y_pixdim"]
    y_dim = info["y_dim"]

    shift_down_mm = y_low / y_pixdim  # mm
    remove_upper_slices = y_low + y_high  # number of slices

    # 1a. translate (i.e. remove posterior margin)
    translated_image = temp_path + (os.path.basename(image)).replace(".nii.gz", "_translated.nii.gz")
    nifti_translate_image(image, translated_image, y_mm=0.0 - shift_down_mm)

    # 1b: perform same shift with activity map
    nifti_translate_image(activity, output_activity, y_mm=0.0 - shift_down_mm)

    # 2. remove upper (i.e., cutoff most anterior slices)
    y_dim_new = y_dim - remove_upper_slices
    nifti_edit_dimensions(translated_image, output_image, y_dim=y_dim_new)

    return [output_image, output_activity]


def fix_shader_problem(file_iteration_main, file_iteration_anatomical, file_iteration_comparison, file_iteration_surface):
    # problem: activity map visible through surface which cannot be fixed by shader settings
    # solution: compare img_iteration_anatomical (no activity overlay) and img_iteration_surface (darker surface with ambient light off) to replace pixels only on the surface (i.e. activity only shown on cut)

    img_iteration_main          = Image.open(file_iteration_main)
    img_iteration_anatomical    = Image.open(file_iteration_anatomical)
    img_iteration_surface       = Image.open(file_iteration_surface)
    img_iteration_comparison    = Image.open(file_iteration_comparison)



    pixels_main                 = img_iteration_main.load()
    pixels_anatomical           = img_iteration_anatomical.load()
    pixels_surface              = img_iteration_surface.load()
    pixels_comparison           = img_iteration_comparison.load()


    width, height = img_iteration_main.size
    for x in range(width):
        for y in range(height):

            if pixels_main[x, y] != (0, 0, 0):
                pixel_anatomical    = pixels_anatomical[x, y]
                pixel_surface       = pixels_surface[x, y]
                pixel_comparison    = pixels_comparison[x, y]

                if pixel_surface[0] < pixel_comparison[0]-2: # Surface area appears darker with lower ambient light. Cut surfaces remain unchanged.
                    pixels_main[x, y] = pixel_anatomical


    img_iteration_main.save(file_iteration_main)
    return file_iteration_main


def image_paste(target, file_source, coords_target, coords_source, file_output=None, resize=None):
    # pastes one image onto another
    #       - target: background image: can be PIL image or string filename => returns filename or image
    #       - file_source: what is pasted
    #       - coords_target: [x_target, y_target]
    #       - coords_source: [[x_min, y_min], [x_max, y_max]]: pastes [x_min, y_min] onto [x_target, y_target]
    #       - resize: float value: resize the source image



    is_file = False
    if isinstance(target, str):
        img_target = Image.open(target)
        is_file = True
    else:
        img_target = target


    img_source = Image.open(file_source)
    cropped_source = img_source.crop((coords_source[0][0], coords_source[0][1], coords_source[1][0], coords_source[1][1]))
    if resize != None:
        cropped_source = cropped_source.resize((int(cropped_source.size[0] * resize), int(cropped_source.size[1] * resize)))


    img_target.paste(cropped_source, (coords_target[0], coords_target[1]))


    if is_file:
        if file_output == None:
            file_output = file_target

        img_target.save(file_output)
        return file_output
    else:
        if file_output != None:
            img_target.save(file_output)
        return img_target



def finalize_sliced_plot (img_joined, anatomical_stripped, activity_masked_stripped, activity_thresholds, cuts, img_joined_finalized):

    image_result = Image.new("RGB", (2700, 1450), "black")

    # paste slices
    coords_target = [0, 0]
    coords_source = [[500, 1700], [3200, 2900]] # coords_source: [[x_min, y_min], [x_max, y_max]]
    image_result = image_paste(image_result, img_joined, coords_target, coords_source)

    # paste color bar
    coords_target = [490, 1200]
    coords_source = [[130, 70], [3840, 450]]
    image_result = image_paste(image_result, img_joined, coords_target, coords_source, resize=0.30)







    # create projection and indication of slice positions
    proj_image      = temp_path + "projection_anatomical.nii.gz"
    proj_activity   = temp_path + "projection_activity.nii.gz"
    proj_PNG        = temp_path + "projection_activity.png"

    if False: # test: width of image is orignal y span.
        proj_area = temp_path + "projection_area.nii.gz"
        PNG_area = temp_path + "projection_area.png"
        run("fslmaths " + anatomical_stripped + " -mul 0 -rand " + proj_area, no_sim=True)
        mricrogl_plot_activity(proj_area, proj_activity, activity_thresholds, output_filename=PNG_area, orientation='sagittal', cuts=[0.5])

    run("fslmaths " + anatomical_stripped       + " -Xmean " + proj_image,      no_sim=True)
    run("fslmaths " + activity_masked_stripped  + " -Xmax " + proj_activity,   no_sim=True)


    mricrogl_plot_activity(proj_image, proj_activity, activity_thresholds, output_filename=proj_PNG, orientation='sagittal-noseleft', cuts=[0.5])

    # width of image is correct and corresponds to original y dimensions (tested!)
    # now: remove black margin from the top

    sagittal_img = Image.open(proj_PNG)
    pixels_sagittal = sagittal_img.load()
    width_sagittal, height_sagittal = sagittal_img.size


    for y in range(height_sagittal):
        if all(sagittal_img.getpixel((x, y)) == (0, 0, 0) for x in range(width_sagittal)):
            sagittal_img = sagittal_img.crop((0, y, width_sagittal, height_sagittal))
        else:
            break


    height_sagittal_new = 1290
    #sagittal_img = sagittal_img.resize((width_sagittal, height_sagittal_new), Image.LANCZOS)
    #sagittal_img = sagittal_img.crop((0, 0, width_sagittal, height_sagittal_new))





    # Draw vertical lines on the image
    draw = ImageDraw.Draw(sagittal_img)
    line_color = col_reds[0] # (255, 0, 0)  # Red color
    line_thickness = 6

    for cut in cuts:
        x = int(width_sagittal * (1.0-cut))
        draw.line((x, 0, x, height_sagittal), fill=line_color, width=line_thickness)

    sagittal_img.save(proj_PNG)

    # paste projection
    coords_target = [2000, 1000]
    coords_source = [[0, 0], [width_sagittal, height_sagittal_new]]
    image_result = image_paste(image_result, proj_PNG, coords_target, coords_source, resize=0.30)
    image_result.save(img_joined_finalized)


    return img_joined_finalized




