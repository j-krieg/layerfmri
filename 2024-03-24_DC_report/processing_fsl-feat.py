
############
# feat


def handler(signum, frame):
    try:
        with open (output_path + "/feat_results/..feat/report_log.html", "r") as file:
            html_content = file.read()
    except FileNotFoundError:
        html_content = ""

    soup = BeautifulSoup(html_content, 'html.parser')
    text = soup.get_text()
    text = '\n'.join(line for line in text.splitlines() if line.strip())

    try:
        with open(output_path + 'feat_results/start_time.txt', 'r') as file:
            start_time = int(file.read())
    except Exception as e:
        start_time = 0

    now_time = time.time()
    duration = int(now_time - start_time)
    time_format = sec2HHMMSS(duration)
    text += "\n\n running FSL-FEAT for " + str(time_format)

    print("\033[H\033[J")  # Clear the console output

    for keyword in ["ERROR MESSAGE", "abnormally", "Error", "error", "exited"]:
        if keyword in text:
            print(text + "\n\n\n### An error occured in FSL-FEAT.")
            quit()

    print(text)

    signal.alarm(1)  # Set up the next alarm


def run_feat(cmd): # output will be [[output_path + "activity_map.nii.gz"]]
    with open(output_path + 'feat_results/start_time.txt', 'w') as file:
        file.write(str(int(time.time())))

    signal.signal(signal.SIGALRM, handler)
    signal.alarm(1)
    output = os.popen(cmd).read()
    signal.alarm(0)
    print(output)

    jobs = []
    jobs.append("rm -rf -f" + output_path + "feat_results/feat_output")
    jobs.append("cp -r " + output_path + "feat_results/..feat/* " + output_path + "feat_results/feat_output")
    jobs.append("rm -rf -f" + output_path + "feat_results/..feat")

    if True: # delete large files that are not needed
        large_files = 	["filtered_func_data.nii.gz",
                         "prefiltered_func_data.nii.gz",
                         "prefiltered_func_data_thresh.nii.gz",
                         "stats/res4d.nii.gz",
                         "stats/threshac1.nii.gz"]

        for file in large_files:
            jobs.append("rm " + output_path + "feat_results/feat_output/" + file)


    jobs.append("cp " + output_path + "feat_results/feat_output/stats/cope1.nii.gz " + output_path + "activity_map.nii.gz")

    for job in jobs:
        run(job)

    print("\n\nFSL-FEAT completed successfully.")


def feat_processing(functional):
    #########################################
    # 1. Define inputs, settings, and outputs

    folder_current = str(os.getcwd() + "/")
    template_settings = "settings/template_design.fsf"

    # input files
    input_fmri = folder_current + functional # must be an absolute path for feat
    activity = ""
    if "vistim" in input_fmri:
        activity = "vistim"
        total_volumes = 150

    if "rft" in input_fmri:
        activity = "rft"
        total_volumes = 120

    templates = {"vistim":"own_template_vistim.tsv", "rft":"own_template_rft.tsv"}

    if not activity in templates:
        print("\n\n###Error: activity not found: " + activity)
        quit()
    else:
        input_task_events = "settings/" + templates[activity]

    repetition_time = 1.96 # number of seconds between two slices. => e.g., fslinfo pixdim4


    print("\n\ninput_task_events: " + str(input_task_events))

    # task_events has 3 columns: onset time, duration, and trial type. These columns are different from the input needd by feat
    # onset	duration	trial_type
    # 0     44.2        rest
    # 44.2  31.2        right_finger_tapping
    # 75.4  31.2        rest
    # 106.6 31.2        right_finger_tapping
    # 137.8 31.2        rest


    folder_output = folder_current + output_path + "feat_results/"  # must be an absolute path for feat
    os.popen("rm -rf -f " + folder_output).read()
    os.makedirs(folder_output)
    custom_design = folder_output + "design.fsf"

    #########################################
    # 2. Create the two .tsv files needed by feat. Note that "EV1.tsv" and "EV2.tsv" are different from, but generated based on "task_events.tsv":
    #   - separate files for active and resting phases
    #   - no header
    #   - different columns (1. onset time, 2. duration, 3. weighting factor)
    # "EV1.tsv": resting events, "EV2.tsv": activity events


    with open(input_task_events, "r") as file:
        lines = file.readlines()
        lines = lines[1:]  # remove header in the first line
        events = []
        for line in lines:
            columns = line.split()
            events.append(columns)

    events_resting = []  # for EV1.txt
    events_activity = []  # for EV2.txt

    for event in events:
        event_modified = [element if index != 2 else "1.0" for index, element in enumerate(
            event)]  # replace 3rd element by "1.0", because feat interprets this as a weighting factor
        if "rest" in event[2]:  # trial type
            events_resting.append(event_modified)  # e.g., "rest"
        else:
            events_activity.append(event_modified)  # e.g., "right_finger_tapping"

    EV1_file = folder_output + "EV1.txt"
    EV2_file = folder_output + "EV2.txt"

    for TSV_file, events in zip([EV1_file, EV2_file], [events_resting, events_activity]):
        with open(TSV_file, "w+") as file:
            for event in events:
                file.write('\t'.join(event) + '\n')

    #########################################
    # 3. create design.fsf from template

    with open("settings/template_design.fsf", "r") as file:
        settings_template = file.read()

        settings = settings_template.replace("<output_folder>", folder_output.rstrip("/"))  # without trailing "/"
        settings = settings.replace("<TSV_EV1>", EV1_file)
        settings = settings.replace("<TSV_EV2>", EV2_file)
        settings = settings.replace("<fmri_input>", input_fmri)
        settings = settings.replace("# ONLY TEMPLATE, MUST BE CHANGED.", "# ADAPTED FROM TEMPLATE.")
        settings = settings.replace("<total_volumes>", str(total_volumes)) # 120 for rft, 150 for vistim
        settings = settings.replace("<repetition_time>", "{:.6f}".format(repetition_time))

    with open(custom_design, "w+") as file:
        file.write(settings)

    feat_command = "feat " + custom_design
    # os.popen(feat_command).read()
    run_feat(feat_command)

