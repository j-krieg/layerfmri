




def plot_multiple_graphs (xvals, list_yvals, output_file, list_y_errs=None, axis_title_x=None, axis_title_y = None, list_labels = None, plot_title=None, dict_xticks=None, list_linestyles=None, list_colors=None, label_sorting=None):

    # xvals: [1, 2, 3] # same for all graphs
    # yvals: [[9, 8, 7], [5, 5, 5]] # example for two graphs
    # color: None for default color scheme
    # linestyle: supported values are '-', '--', '-.', ':', 'None', ' ', '', 'solid', 'dashed', 'dashdot', 'dotted'
    # label_sorting: can be None/False, true, or a list of integer numbers with the same length as xvals. In case of True, the labels will be sorted by their graph mean values. If a list of integers is given, mean values are taken as a secondary criterium if some integers are same


    # use LaTeX fonts for plots:
    #   already done:
    #       install latex on system:
    #           sudo apt update
    #           sudo apt install texlive-xetex texlive-latex-recommended texlive-fonts-recommended
    #       import matplotlib
    #       matplotlib.use("pgf")  # use LaTeX fonts for plots
    #       import matplotlib.pyplot as plt

    plt.rcParams.update({
        "pgf.texsystem": "pdflatex",
        'font.family': 'sans-serif',
        'font.size': 12,
        'text.usetex': True,
        'pgf.rcfonts': False,
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'legend.fontsize': 10,
    })




    plt.clf()  # initialize
    plt.figure(figsize=(8, 5))
    
    average_vals_for_eventual_sorting = []

    for n, yvals in enumerate(list_yvals):
        if list_labels != None:
            label = list_labels[n]
        else:
            label = None


        if list_linestyles != None:
            linestyle = list_linestyles[n]
        else:
            linestyle = None


        if list_colors != None:                         # color list given
            probe = list_colors[n]
            if probe != None:                           # color list given AND current color given
                color = probe
            else:                                       # color list given, but current color not given
                color = col_seq[n%col_sec_entries]      # use own color scheme if no color is given.
        else:
            color = col_seq[n%col_sec_entries]          # no color list given


        plt.plot(xvals, yvals, label=label, linestyle=linestyle, color=unify_color(color))

        if list_y_errs != None: # activate/deactivate errors for all graphs
            yerr = list_y_errs[n]

            if yerr != None: # activate/deactivate errors for individual graphs
                yerr_array = np.array(yerr)
                plt.fill_between(xvals,
                                 yvals - yerr_array, yvals + yerr_array,
                                 alpha=0.5, edgecolor=unify_color(color), facecolor=adjust_color(color, 0.1, 0.2), label=None)


        # now calculate average values in case the labels should be sorted:
        avg, err = calculate_average_error(yvals)
        average_vals_for_eventual_sorting.append(avg)






    if plot_title != None:
        plt.title(plot_title, fontsize=8)


    if axis_title_x == None:
        plt.xlabel('')
    else:
        plt.xlabel(axis_title_x)


    if axis_title_y == None:
        plt.ylabel('')
    else:
        plt.ylabel(axis_title_y)


    if dict_xticks != None:                         # name xticks individually: dict with x-value and tile
        # e.g., dict_xticks = {'WM':1, 'pial':9}
        list_names = []
        list_values = []
        for name, value in dict_xticks.items():
            list_names.append(name)
            list_values.append(value)

        plt.xticks(list_values, list_names)
        # e.g., plt.xticks([min(cortical_depth), max(cortical_depth)], ['WM', 'pial'])


    plt.grid(axis='both', which='major', linestyle=':', linewidth=0.5, alpha=0.5)
    plt.minorticks_on()
    plt.grid(axis='x', which='both', linestyle=':', linewidth=0.5, alpha=0.5)

    # Set the spacing of minor ticks
    minor_locator = MultipleLocator(1)
    ax = plt.gca()
    ax.xaxis.set_minor_locator(minor_locator)


    # if labels are provided create a legend (i.e., make plot a 20% smaller and add legend on the right)
    if list_labels != None:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        loc="upper left"
        plt.legend(list_labels) # only a dummy, legend will be replotted after labels are checked.

        # Remove labels/handles with value None from the legend (typically, the None labels created when plotting error margins)
        handles = plt.gca().get_legend().legend_handles
        labels = plt.gca().get_legend().get_texts()

        filtered_handles = []
        filtered_labels = []

        for handle, label in zip(handles, labels):
            if label.get_text() != '':
                filtered_handles.append(handle)
                filtered_labels.append(label.get_text())





        # now: sort the labels if label_sorting is True or a list of integer numbers:
        if label_sorting != False and label_sorting != None:
            if label_sorting == True:
                first_sorting_criteria = [1] * len(xvals) # no other list given, only "True" was given as a parameter
            else:
                first_sorting_criteria = label_sorting # list of integer values

            sorting_criteria=[]
            for n, avg in enumerate(average_vals_for_eventual_sorting):
                sorting_criteria.append (100*max(average_vals_for_eventual_sorting)*first_sorting_criteria[n] - avg) # use both criteria: given integer numer (first priority, ascending), average value (second priority, descending)


            filtered_handles = aligned_sort(filtered_handles, sorting_criteria)
            filtered_labels = aligned_sort(filtered_labels, sorting_criteria)



                

        plt.legend(filtered_handles, filtered_labels, bbox_to_anchor=(1, 1.02), loc=loc) # 1.02 to align upper boarders of diagram and legend box



    plt.savefig(output_file, dpi=300)

    if True: # get printable versions for report
        output_PDF              = output_file.replace(".png", ".pdf")
        output_PDF_notitle      = output_file.replace(".png", "_no-title.pdf")

        plt.savefig(output_PDF)
        plt.title("")
        plt.savefig(output_PDF_notitle)

    plt.clf()
    plt.close()







def plot_overlay (anatomical, projection, output_png, title="", type="activity", view="zscore", alpha=1.0, threshold = None, upsampling = None):
    # plots an overlay (activity map or layers) on anatomical image
    # input:
    #           anatomical: anatomical image
    #           projection: activity map or layers map
    #           type: "activity", "layers", "3layers", "5layers", "7layers":    "activity" (diverging color map) or "layers" (qualitative color map, e.g., tab10, FlatUI, ggplot, default, 538)
    #           view: "zscore", "x", "y", "z", "ortho"
    #           upsampling (default None): e.g., 0.3 if background image should be upsampled to 0.3mm for a better resolution of the resulting image.
    #

    ###############################################################
    # plot z-score maps with DC map projection

    if anatomical == None: # only show overlay (plotting.plot_stat_map accepts no "None", therefore create empty image with the size of the overlay)
        anatomical_black = naming(projection, "black", temp_path)
        run("fslmaths " + projection + " -mul 0 " + anatomical_black, no_sim=True)
        anatomical = anatomical_black

    if upsampling != None:
        [anatomical] = resampling (upsampling, [anatomical])

    bg_img = image.load_img(anatomical)      # anatomical
    img = image.load_img(projection)         # overlay

    # type can be "activity" or "layers". "layers" may have an additional number. This is used to set a diverging color map (for activity) or a qualitative color map (for layers)
    vmin = None
    vmax = None
    cbar_tick_format = None
    colorbar = None
    if (threshold== None):
        threshold = 0.00001  # 0 should be transparent.

    if type == "activity":
        color_map = "inferno"
        cbar_tick_format = "%.2g"
        colorbar = True
    else:
        if not "layers" in type:
            print("\n\nError: type unknown: " + type)
            quit()
        else:
            # plotting layers
            cbar_tick_format = "%i"

            color_map = "tab10"                         # alternative: 'Accent'
            vmin = 1                                    # tab10 has 10 different colors
            vmax = 10
            if "3" in type:                             # e.g., "3layers"
                color_map = "Accent"    # "Set1"
                vmin = 0                                 # Set1 has 9 different colors
                vmax = 8
                colorbar = False

            if "5" in type or "7" in type:              # e.g., "7layers"
                color_map = "Set2"
                vmin = 1
                vmax = 8                                # Set2 has 8 different colors
                colorbar = True

            if "15" in type:
                color_map = "tab20"        # alternative: Set3 has  12
                vmin = 1
                vmax = 20                                # Set2 has 8 different colors
                colorbar = True



    if view not in ["zscore", "ortho", "x", "y", "z"]:
        print("\n\nError: view unknown: " + view)
        quit()

    if view == "zscore":
        display_mode = "z"
        cut_coords = 5

    if view in ["x", "y", "z"]:
        display_mode = view
        cut_coords=1

    if view == "ortho":
        display_mode = "ortho"
        cut_coords = None #[36, -27, -30]


    # documentation: https://nilearn.github.io/dev/modules/generated/nilearn.plotting.plot_stat_map.html
    plotting.plot_stat_map(img, bg_img=bg_img, display_mode=display_mode, cut_coords=cut_coords, alpha=alpha, cmap=color_map, vmin=vmin, vmax=vmax, cbar_tick_format=cbar_tick_format, draw_cross=False, colorbar=colorbar, title=title, output_file=output_png, threshold=threshold)





def brainspaces_schaefer100_single(values, minimum, maximum):
    # example use: image = brainspaces_schaefer100_single ([0, 1, 2, ..., 99], 0, 1000)
    # input: list of 100 values that will be presented on the schaefer100 parcellated surface
    # minimum and maximum are needed for normalization in case plots are generated for multiple layers:
    # returns a figure that can be used like that:
    #       image.save('plot.png')

    surfaces = fetch_fslr()  # surface files included in the neuromaps package
    lh, rh = surfaces['inflated']

    # atlas selection:
    # 		https://github.com/brainspaces/schaefer400
    # 		https://github.com/brainspaces/schaefer100
    # 		https://github.com/brainspaces/glasser360

    atlas = load_parcellation('schaefer', 100, join=True)  # location: /home/user/fsl/lib/python3.11/site-packages/brainspace/datasets/parcellations/
    lh_parc, rh_parc = load_parcellation('schaefer', 100)

    unique = np.unique(atlas)  # get unique values in the 'atlas' array
    unique = unique[1:len(unique)]  # remove the first element from the 'unique' array
    # print(len(atlas)) 			# 64984
    # print (str(atlas))			# [50 21 13 ... 91 91 91]
    # print(str(unique.shape[0])) # 100 (from schaefer-100)

    for i in range(100):
        value = values[i]
        atlas = np.where(atlas == unique[i], value, atlas)

    # Generate plot
    p = Plot(lh, rh, views=['lateral', 'medial'], size=(2000, 500), zoom=1.2, layout='row')
    p.add_layer(atlas, cbar=True, cmap='viridis', color_range=(minimum, maximum))
    p.add_layer({'left': lh_parc, 'right': rh_parc}, cmap='gray', as_outline=True, cbar=False)  # add outline

    fig = p.build()
    # fig.show()
    # fig.savefig('plot.png')

    # convert fig to image
    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', pad_inches=0)
    buf.seek(0)
    img = Image.open(buf)
    # img.save("demo.png")
    return img





def brainspaces_schaefer100_stacked(list_of_sets, output_PNG):
    # input: list of sets (each with 100 intensity values), output, filename
    # example use: rainspaces_schaefer100_stacked ([[0, 1, 2, ..., 99], [0, 1, 2, ..., 99], [0, 1, 2, ..., 99]], "output_stacked.PNG")
    number_layers = len(list_of_sets)  # e.g. 3 sets of 100 intensity values to represent 3 different layers

    ##################
    ## introduce same min and max values in datasets, such that the plots are normalized.
    minimum = -1.0
    maximum = -1.0

    for set in list_of_sets:
        if False:
            for value in set:
                if minimum == -1.0 or minimum > value:
                    minimum = value

                if maximum == -1.0 or maximum < value:
                    maximum = value
        else:
            # avoid outliers

            percentile_low = np.percentile(np.array(set), 5)
            percentile_high = np.percentile(np.array(set), 95)

            if minimum == -1.0 or minimum > percentile_low:
                minimum = percentile_low

            if maximum == -1.0 or maximum < percentile_high:
                maximum = percentile_high



    # make minimum and maxima look nicer:
    if minimum < 0.2 * maximum:
        minimum = 0.0
    else:
        minimum = nice_round(minimum, "floor")

    maximum = nice_round(maximum, "ceil")
    ##################

    P = 100  # 100 parcelations (e.g., schaefer-100)
    overlap = 50  # vertical overlap of 100 pixels (to remove the scale)
    number_sets = len(list_of_sets)

    for n, set in enumerate(list_of_sets):
        if len(set) != P:
            print("### error in brainspaces_schaefer100_stacked: need exactly " + str(P) + " parcellations in each set.")
            exit()

        img = brainspaces_schaefer100_single(set, minimum, maximum)


        if n == 0:
            width, height_single = img.size
            height_total = height_single + ((number_sets - 1) * (height_single - overlap))
            stacked_image = Image.new('RGB', (width, height_total))

        stacked_image.paste(img, (0, (n * (height_single - overlap))))  # n starts at 0

    # add text
    margin_top = 30
    margin_left = 250
    height_new = height_total + margin_top
    width_new = width + margin_left

    final_image = Image.new("RGB", (width_new, height_new), "white")
    final_image.paste(stacked_image, (margin_left, margin_top))

    drawing = ImageDraw.Draw(final_image)
    myFont = ImageFont.truetype('FreeMono.ttf', 36)
    myFill = (0, 0, 0)
    drawing.text((190 + margin_left, 10), "LH", font=myFont, fill=myFill)
    drawing.text((1000 + margin_left, 10), "RH", font=myFont, fill=myFill)

    for n in range(number_sets):
        text = str(n + 1)  # in LN2_PROFILE, the first layer is 1.
        drawing.text((25, 210 + (n * (height_single - overlap))), text, font=myFont, fill=myFill)

    offset = 40
    drawing.text((25, 210 - offset + (0 * (height_single - overlap))), "(sWM)", font=myFont, fill=myFill)  # top: WM (smaller layer number)
    drawing.text((25, 210 + offset + ((number_sets - 1) * (height_single - overlap))), "(pial)", font=myFont, fill=myFill)  # bottom: pial (larger layer number)

    final_image.save(output_PNG)








    
def mricrogl_plot_labelfile (anatomical, label_file, output_filename=None, orientation='z', cuts=[0.20, 0.35, 0.5, 0.65, 0.80], skip_GM=False):

    # cuts = [0.20, 0.35, 0.5, 0.65, 0.80] # relative values between min and max on z axis

    if output_filename == None:
        output_filename = label_file.replace(".nii.gz", "_" + orientation + "_MRIcroGL.png")

    info = fsl_info(anatomical)
    x_min = info["x_min"]
    x_max = info["x_max"]
    y_min = info["y_min"]
    y_max = info["y_max"]
    z_min = info["z_min"]
    z_max = info["z_max"]


    label_red = temp_path + "label_red.nii.gz"
    label_blue = temp_path + "label_blue.nii.gz"
    label_white = temp_path + "label_white.nii.gz"


    if label_file != None:
        run("fslmaths " + label_file + " -thr 3 -uthr 3 " + label_red, no_sim=True)       # GM
        run("fslmaths " + label_file + " -thr 1 -uthr 1 " + label_blue, no_sim=True)      # CSF
        run("fslmaths " + label_file + " -thr 2 -uthr 2 " + label_white, no_sim=True)     # WM



    jobs = []
    jobs.append("import gl")
    jobs.append("gl.resetdefaults()")
    jobs.append("gl.shaderquality1to10(10)")
    jobs.append("gl.colorbarposition(0)")
    jobs.append("gl.linewidth(0)")
    jobs.append("gl.smooth(0)")
    jobs.append("gl.bmpzoom(3)")


    jobs.append("gl.loadimage('<anatomical>')")
    if label_file != None:
        jobs.append("gl.overlayload('<label_red>')")
        jobs.append("gl.overlayload('<label_blue>')")
        jobs.append("gl.overlayload('<label_white>')")

        jobs.append("gl.minmax(1, 0.5, 3.2)")
        if skip_GM: # do not show GM, only lines.
            jobs.append("gl.opacity(1,0)")
        else:
            jobs.append("gl.opacity(1,100)")
        jobs.append("gl.colorname (1,'1red')")

        jobs.append("gl.minmax(2, 0.5, 1.2)")
        jobs.append("gl.opacity(2,100)")
        jobs.append("gl.colorname (2,'3blue')")

        jobs.append("gl.minmax(3, 0.5, 0.7)")
        jobs.append("gl.opacity(3,100)")
        jobs.append("gl.colorname (3,'bone')")



    jobs.append("x_min = <x_min>")
    jobs.append("x_max = <x_max>")
    jobs.append("y_min = <y_min>")
    jobs.append("y_max = <y_max>")
    jobs.append("z_min = <z_min>")
    jobs.append("z_max = <z_max>")
    jobs.append("cuts = <cuts>")

    jobs.append("for n, cut in enumerate(cuts):")
    jobs.append("\tx = ((1.0-cut)*x_min) + (cut*x_max)")
    jobs.append("\ty = ((1.0-cut)*y_min) + (cut*y_max)")
    jobs.append("\tz = ((1.0-cut)*z_min) + (cut*z_max)")
    jobs.append("\tgl.orthoviewmm(x,y,z)")

    # Axial (1), Coronal (2), Sagittal (4)
    if orientation in ['z', 'axial']:
        jobs.append("\tgl.view(1)")

    if orientation in ['y', 'coronal']:
        jobs.append("\tgl.view(2)")

    if orientation in ['x', 'sagittal']:
        jobs.append("\tgl.view(4)") # (sic!)




    jobs.append("\tgl.savebmp('<temp/>mricrogl_output_' + str(n+1) + '.png')")
    jobs.append("quit()")

    script = ""
    for job in jobs:
        job = job.replace("<x_min>", str(x_min))
        job = job.replace("<x_max>", str(x_max))
        job = job.replace("<y_min>", str(y_min))
        job = job.replace("<y_max>", str(y_max))
        job = job.replace("<z_min>", str(z_min))
        job = job.replace("<z_max>", str(z_max))
        job = job.replace("<cuts>", str(cuts))
        job = job.replace("<temp/>", temp_path)
        job = job.replace("<label_red>", label_red)     # GM
        job = job.replace("<label_blue>", label_blue)   # CSF
        job = job.replace("<label_white>", label_white) # WM
        job = job.replace("<anatomical>", anatomical)

        script += job + "\n"

    script_file = temp_path + "MRIcroGL_script.py"
    with open(script_file, "w") as file:
        file.write(script)
    run ("~/MRIcroGL/MRIcroGL " + script_file, no_sim=True)



    #########
    # join images
    img_list = []
    for n, cut in enumerate(cuts):
        input_file = (temp_path + "mricrogl_output_" + str(n+1) + ".png")
        input = Image.open(input_file)

        if n == 0:
            width, height = input.size
            width_new = width * len(cuts)
            height_new = height
            final_image = Image.new("RGB", (width_new, height_new), "black")

        final_image.paste(input, (n * width, 0))





    target_color = (236, 0, 0)
    new_color = (205, 62, 78)
    final_image = change_color(final_image, target_color, new_color)

    target_color = (0, 0, 182)
    new_color = (70, 130, 180)
    final_image = change_color(final_image, target_color, new_color)

    final_image.save(output_filename)



    return output_filename




def mricrogl_plot_layerfile (img_anatomical, img_layers, number_layers, output_filename=None, orientation='z', cuts=[0.20, 0.35, 0.5, 0.65, 0.80]):
    # Projects a layer file (intensities 1...N) onto an anatomical image and creates z plots.
    # Works only for layer files. Use a different function for labels or activities.

    #cuts = [0.20, 0.35, 0.5, 0.65, 0.80] # relative values between min and max on z axis

    if output_filename == None:
        output_filename = img_layers.replace(".nii.gz", "_MRIcroGL.png")


    info = fsl_info(img_anatomical)
    x_min = info["x_min"]
    x_max = info["x_max"]
    y_min = info["y_min"]
    y_max = info["y_max"]
    z_min = info["z_min"]
    z_max = info["z_max"]


    jobs = []
    jobs.append("import gl")
    jobs.append("gl.resetdefaults()")
    jobs.append("gl.shaderquality1to10(10)")
    jobs.append("gl.colorbarposition(0)")
    jobs.append("gl.linewidth(0)")
    jobs.append("gl.smooth(0)")
    jobs.append("gl.bmpzoom(3)")


    jobs.append("gl.loadimage('<anatomical>')")
    jobs.append("gl.overlayload('<layers>')")



    jobs.append("gl.minmax(1, 0, <number_layers>)") # must tstart with "0", otherwise the outermost layer will not be displayed.
    jobs.append("gl.opacity(1,100)")
    jobs.append("gl.colorname (1,'NIH')")





    jobs.append("x_min = <x_min>")
    jobs.append("x_max = <x_max>")
    jobs.append("y_min = <y_min>")
    jobs.append("y_max = <y_max>")
    jobs.append("z_min = <z_min>")
    jobs.append("z_max = <z_max>")
    jobs.append("cuts = <cuts>")

    jobs.append("for n, cut in enumerate(cuts):")
    jobs.append("\tx = ((1.0-cut)*x_min) + (cut*x_max)")
    jobs.append("\ty = ((1.0-cut)*y_min) + (cut*y_max)")
    jobs.append("\tz = ((1.0-cut)*z_min) + (cut*z_max)")
    jobs.append("\tgl.orthoviewmm(x,y,z)")



    # Axial (1), Coronal (2), Sagittal (4)
    if orientation in ['z', 'axial']:
        jobs.append("\tgl.view(1)")

    if orientation in ['y', 'coronal']:
        jobs.append("\tgl.view(2)")

    if orientation in ['x', 'sagittal']:
        jobs.append("\tgl.view(4)")  # (sic!)



    jobs.append("\tgl.savebmp('<temp/>mricrogl_output_' + str(n+1) + '.png')")
    jobs.append("quit()")

    script = ""
    for job in jobs:
        job = job.replace("<x_min>", str(x_min))
        job = job.replace("<x_max>", str(x_max))
        job = job.replace("<y_min>", str(y_min))
        job = job.replace("<y_max>", str(y_max))
        job = job.replace("<z_min>", str(z_min))
        job = job.replace("<z_max>", str(z_max))
        job = job.replace("<cuts>", str(cuts))
        job = job.replace("<temp/>", temp_path)
        job = job.replace("<layers>", img_layers)
        job = job.replace("<anatomical>", img_anatomical)
        job = job.replace("<number_layers>", str(number_layers))

        script += job + "\n"

    script_file = temp_path + "MRIcroGL_script.py"
    with open(script_file, "w") as file:
        file.write(script)
    run ("~/MRIcroGL/MRIcroGL " + script_file, no_sim=True)



    #########
    # join images
    img_list = []
    for n, cut in enumerate(cuts):
        input_file = (temp_path + "mricrogl_output_" + str(n+1) + ".png")
        input = Image.open(input_file)

        if n == 0:
            width, height = input.size
            width_new = width * len(cuts)
            height_new = height
            final_image = Image.new("RGB", (width_new, height_new), "black")

        final_image.paste(input, (n * width, 0))


    final_image.save(output_filename)
    return output_filename



def mricrogl_plot_activity (img_anatomical, img_activity, activity_thresholds, output_filename=None, orientation='z', cuts=[0.20, 0.35, 0.5, 0.65, 0.80]):
    if output_filename == None:
        output_filename = img_activity.replace(".nii.gz", "_MRIcroGL.png")

    [activity_min, activity_max] = activity_thresholds

    info = fsl_info(img_anatomical)
    x_min = info["x_min"]
    x_max = info["x_max"]
    y_min = info["y_min"]
    y_max = info["y_max"]
    z_min = info["z_min"]
    z_max = info["z_max"]


    jobs = []
    jobs.append("import gl")
    jobs.append("gl.resetdefaults()")
    jobs.append("gl.shaderquality1to10(10)")
    jobs.append("gl.colorbarposition(0)")
    jobs.append("gl.linewidth(0)")
    jobs.append("gl.smooth(0)")
    jobs.append("gl.bmpzoom(3)")


    jobs.append("gl.loadimage('<anatomical>')")
    jobs.append("gl.overlayload('<activity>')")

    jobs.append("gl.minmax(1, <activity_min>, <activity_max>)")
    jobs.append("gl.opacity(1,100)")
    jobs.append("gl.colorname (1,'NIH')")





    jobs.append("x_min = <x_min>")
    jobs.append("x_max = <x_max>")
    jobs.append("y_min = <y_min>")
    jobs.append("y_max = <y_max>")
    jobs.append("z_min = <z_min>")
    jobs.append("z_max = <z_max>")
    jobs.append("cuts = <cuts>")

    jobs.append("for n, cut in enumerate(cuts):")
    jobs.append("\tx = ((1.0-cut)*x_min) + (cut*x_max)")
    jobs.append("\ty = ((1.0-cut)*y_min) + (cut*y_max)")
    jobs.append("\tz = ((1.0-cut)*z_min) + (cut*z_max)")
    jobs.append("\tgl.orthoviewmm(x,y,z)")



    # Axial (1), Coronal (2), Sagittal (4)
    if orientation in ['z', 'axial']:
        jobs.append("\tgl.view(1)")

    if orientation in ['y', 'coronal']:
        jobs.append("\tgl.view(2)")

    if orientation in ['x', 'sagittal']:
        jobs.append("\tgl.view(4)")  # (sic!)

    if orientation in ['x', 'sagittal-noseleft']:
        jobs.append("\tgl.view(8)")



    jobs.append("\tgl.savebmp('<temp/>mricrogl_output_' + str(n+1) + '.png')")
    jobs.append("quit()")

    script = ""
    for job in jobs:
        job = job.replace("<x_min>", str(x_min))
        job = job.replace("<x_max>", str(x_max))
        job = job.replace("<y_min>", str(y_min))
        job = job.replace("<y_max>", str(y_max))
        job = job.replace("<z_min>", str(z_min))
        job = job.replace("<z_max>", str(z_max))
        job = job.replace("<cuts>", str(cuts))
        job = job.replace("<temp/>", temp_path)
        job = job.replace("<activity>", img_activity)
        job = job.replace("<anatomical>", img_anatomical)
        job = job.replace("<activity_min>", str(activity_min))
        job = job.replace("<activity_max>", str(activity_max))


        script += job + "\n"

    script_file = temp_path + "MRIcroGL_script.py"
    with open(script_file, "w") as file:
        file.write(script)
    run ("~/MRIcroGL/MRIcroGL " + script_file, no_sim=True)



    #########
    # join images
    img_list = []
    for n, cut in enumerate(cuts):
        input_file = (temp_path + "mricrogl_output_" + str(n+1) + ".png")
        input = Image.open(input_file)

        if n == 0:
            width, height = input.size
            width_new = width * len(cuts)
            height_new = height
            final_image = Image.new("RGB", (width_new, height_new), "black")

        final_image.paste(input, (n * width, 0))


    final_image.save(output_filename)
    return output_filename




def plot_activity_in_slices (anatomical, activity):

    spread          = 2.0 # spread >= 1: typically a number between 1 and 2 - spreads out the slices (1 = not at all, 2: by 200% of the original length in y direction)
    depth_range     = [0.05, 0.95]
    cuts = [0.05, 0.08, 0.11, 0.14, 0.30, 0.40, 0.42, 0.45, 0.47] # 0.0 = posterior, 1.0 anterior in original anatomical image
    number_slices = len(cuts)




    # determine thresholds for activity map:
    [min_intensity, max_intensity] = activity_thresholding_for_plotting(activity)

    info = fsl_info(anatomical)
    y_pixdim = info["y_pixdim"] # mm
    y_dim = info["y_dim"]       # slices
    y_dim_mm = y_dim * y_pixdim


    y_dim_new = int(spread * y_dim)
    y_dim_new_mm = y_dim_new * y_pixdim



    activity_anatomicalspace = transform_space(activity, anatomical, sampling="Li")     # required for "strip_image_in_y_direction"
    activity_masked = apply_mask(activity_anatomicalspace, anatomical, erode_voxels=True)        # no activity outside of the volume

    anatomical_stripped = temp_path + "image_stripped.nii.gz"
    processed_image = temp_path + "image_slices.nii.gz"
    processed_activity = temp_path + "activity_slices.nii.gz"

    [anatomical_stripped, activity_masked_stripped] = strip_image_in_y_direction (anatomical, activity_masked)


    nifti_edit_dimensions(anatomical_stripped, processed_image,    y_dim=y_dim_new) # with extended header such that there is more spce to move the volume towards anterior
    nifti_edit_dimensions(activity_masked,     processed_activity, y_dim=y_dim_new)




    
    result_list=[] # individual images of each slice need to be joined later
    for n, cut in enumerate(cuts):
        # cut (0..1) is the intended position in the anatomical image
        # cut_mm is the absolute position of the current cut in the original anatomical image
        # slicer_position_rel = relative slicer position between 0 (first cut) and 1 (last cut)
        # slicer_position_mm = absolute slicer position in teh current scenario
        # translation_mm = shift (from posterior to anterior) of the volume needed such the inteded cut occurs at the current slicer position

        cut_mm                      = cut * y_dim_mm
        slicer_position01           = pow(n/(number_slices-1.0), 1.3)                                               # over-represent smaller values: slices can be closer together at the right (farer) side of the image
        slicer_position_rel         = depth_range [0] + (depth_range[1]-depth_range[0]) * slicer_position01         # 0..1 relative to y_dim_new_mm
        depth                       = 1.0 - slicer_position_rel                                                     # slicer moves from 1 to 0 (posterior to anterior)
        slicer_position_mm          = y_dim_new_mm * slicer_position_rel
        translation_mm              = slicer_position_mm - cut_mm                                                   # translate volume such that target cut position and slicer position are aligned

        print ("\n\n n = " + str(n))
        print("cut: " + str(cut))
        print("cut_mm: " + str(cut_mm))
        print("slicer_position_rel: " + str(slicer_position_rel))
        print("depth: " + str(depth))
        print("slicer_position_mm: " + str(slicer_position_mm))
        print("translation_mm : " + str(translation_mm ))


        translated_image = (temp_path + "shifted_image_" + str(n) + ".nii.gz")
        translated_activity = (temp_path + "shifted_activity_" + str(n) + ".nii.gz")
        nifti_translate_image(processed_image,              translated_image,       y_mm=translation_mm) # translate to anterior
        nifti_translate_image(processed_activity,           translated_activity,    y_mm=translation_mm)


        if depth == 0.0:
            depth = 0.01 # do no set depth to 0.0, otherwise full volume will be shown.


        img_iteration_main              = temp_path + "shader_matte-output_" +  str(n) + ".png"         # main image with possible errors on the surface
        img_iteration_anatomical        = temp_path + "shader_matte-anatomical_" + str(n) + ".png"      # no overlay
        img_iteration_comparison        = temp_path + "shader_matte-comparison_" + str(n) + ".png"      # brighter surface
        img_iteration_surface           = temp_path + "shader_matte-surface_" + str(n) + ".png"         # darker surface



        jobs = []
        jobs.append("import gl")
        jobs.append("gl.resetdefaults()")
        jobs.append("gl.loadimage('<image>')")
        jobs.append("gl.overlayload('<overlay>')")
        jobs.append("gl.minmax(1, <activity_min>, <activity_max>)")
        jobs.append("gl.opacity(1,100)")
        jobs.append("gl.colorname (1,'NIH')")
        jobs.append("gl.shaderadjust('boundThresh', 0.35)")
        jobs.append("gl.shaderadjust('edgeThresh', 0.42)")
        jobs.append("gl.shaderadjust('edgeBoundMix',0.05)")
        jobs.append("gl.shaderadjust('colorTemp', 0.8)")
        jobs.append("gl.shadername('Matte')")
        jobs.append("gl.shaderquality1to10(10)")
        jobs.append("gl.shaderadjust('ambient', 1.0)")
        jobs.append("gl.shaderadjust('diffuse', 0.5)")
        jobs.append("gl.shaderadjust('specular', 0.5)")

        jobs.append("gl.shaderadjust('shininess', 1.0)") # likely now working
        jobs.append("gl.shaderadjust('Shininess', 1.0)")

        jobs.append("gl.shaderadjust('overlayFuzzy', 1.0)")
        jobs.append("gl.shaderadjust('overlayDepth', 0.0)")
        jobs.append("gl.shaderadjust('overlayClip', 1.0)")
        jobs.append("gl.linewidth(0)")
        jobs.append("gl.smooth(0)")
        jobs.append("gl.bmpzoom(6)")
        
        if n == 0:
            jobs.append("gl.clipthick(0.5)") # entire piece
        else:
            jobs.append("gl.clipthick(0.008)")
        jobs.append("gl.clipazimuthelevation(<depth>, 0, 180)")

        jobs.append("gl.savebmp('<img_iteration_main>')")

        jobs.append("gl.opacity(1,100)")
        jobs.append("gl.savebmp('<img_iteration_anatomical>')")



        jobs.append("gl.shaderadjust('ambient', 1.0)")
        jobs.append("gl.shaderadjust('diffuse', 0.0)")
        jobs.append("gl.shaderadjust('specular', 0.0)")
        jobs.append("gl.savebmp('<img_iteration_comparison>')")

        jobs.append("gl.shaderadjust('ambient', 0.0)")          # surface is darker now: can compare img_iteration_surface and img_iteration_anatomical to fix shader problem pixelwise
        jobs.append("gl.savebmp('<img_iteration_surface>')")





        jobs.append("quit()")

        script = ""
        for job in jobs:
            job = job.replace("<activity_min>", str(min_intensity))
            job = job.replace("<activity_max>", str(max_intensity))
            job = job.replace("<image>", translated_image)
            job = job.replace("<overlay>", translated_activity)
            job = job.replace("<depth>", str(round(depth, 2)))
            job = job.replace("<temp/>", temp_path)
            job = job.replace("<img_iteration_main>", img_iteration_main)
            job = job.replace("<img_iteration_anatomical>", img_iteration_anatomical)
            job = job.replace("<img_iteration_surface>", img_iteration_surface)
            job = job.replace("<img_iteration_comparison>", img_iteration_comparison)
            script += job + "\n"

        script_file = temp_path + "MRIcroGL_script.py"
        with open(script_file, "w") as file:
            file.write(script)
        run("~/MRIcroGL/MRIcroGL " + script_file, no_sim=True)

        img_iteration_main = fix_shader_problem(img_iteration_main, img_iteration_anatomical, img_iteration_comparison, img_iteration_surface)

        os.remove(translated_image)
        os.remove(translated_activity)
        #os.remove(img_iteration_anatomical)
        #os.remove(img_iteration_surface)

        result_list.append(img_iteration_main)

    for n, single_slice in enumerate(result_list):
        if n == 0:
            img_main = Image.open(single_slice)
            pixels_main = img_main.load()
            width, height = img_main.size
        else:
            img_top = Image.open(single_slice)
            pixels_top = img_top.load()

            for x in range(width):
                for y in range(height):
                    pixel_top = pixels_top[x, y]

                    if pixel_top != (0, 0, 0):
                        pixels_main[x, y] = pixel_top

    img_joined = temp_path + "slices_joined.png"
    img_main.save(img_joined)

    img_joined_finalized = temp_path + "slices_joined_finalized.png"
    finalize_sliced_plot (img_joined, anatomical_stripped, activity_masked_stripped, img_joined_finalized)




        