
import os, re
import phylo_tools as pt
import pandas as pd
#os.system("cat ~/Phylo_Evol_Timeseries/data/illumina_runs/161213/*/SampleSheet.csv > SampleSheets_161213.csv")
#os.system("cat ~/Phylo_Evol_Timeseries/data/illumina_runs/170303/*/SampleSheet.csv > SampleSheets_170303.csv")
#os.system("cat ~/Phylo_Evol_Timeseries/data/illumina_runs/170623/*/SampleSheet.csv > SampleSheets_170623.csv")
#os.system("cat ~/Phylo_Evol_Timeseries/data/illumina_runs/170721/*/SampleSheet.csv > SampleSheets_170721.csv")

# keep: run name, original sample name, new sample name, treatment, strain, replicate, day, index1, index2

# clean HCGS sheets
runs = ['SampleSheets_HCGB_161213', 'SampleSheets_HCGB_170303',
        'SampleSheets_HCGB_170623', 'SampleSheets_HCGB_170721',
        'GSF2056-run1-plates1-2-demultiplexing-summary',
        'GSF2056-run2-plates3-4-demultiplexing-summary',
        'SampleSheet-GSF2124-run3-plates1-2',
        'GSF2124 Lennon Run 3 Plates 3-4 Run Summary Sorted',
        'GSF2124-run5-plates5-6-demultiplexing-summay']

#library_folders = ['161213', '170303', '170623', '170721','GSF2056_run1',
#                    'GSF2056_run2', 'GSF2124_run3', 'GSF2124_run4',
#                    'GSF2124_run5', ]


def get_sample_names():
    run_dir_path = '/N/dc2/projects/muri2/Task2/Phylo_Evol_Timeseries'
    run_dirs = ['161213', '170303', '170623', '170721', 'GSF2056_run1',
                'GSF2056_run2', 'GSF2124_run3', 'GSF2124_run4', 'GSF2124_run5']
    sample_names_path = open(run_dir_path + '/data/illumina_runs/sample_names.txt', 'a')
    for run_dir in run_dirs:
        for path, subdirs, files in os.walk(run_dir_path + '/data/illumina_runs/' + run_dir):
            for name in files:
                new_path = os.path.join(path, name)
                if (run_dir == '161213') or (run_dir == '170303') or (run_dir == '170623') or (run_dir == '170721'):
                    new_path_split = new_path.rsplit('/', 2)
                    if new_path_split[-1] == 'SampleSheet.csv':
                        continue
                    sample_path = run_dir + '/' + new_path_split[1] + '/' + new_path_split[2]
                elif ('GSF' in run_dir):
                    if 'Undetermined' in new_path:
                        continue
                    new_path_split = new_path.rsplit('/', 2)
                    sample_path = new_path_split[1] + '/' + new_path_split[2]
                sample_names_path.write(sample_path + '\n')

    sample_names_path.close()


#def clean_metadata_runs():
#    for run in runs:


def merge_metadata():
    # first get dictionary for barcodes one and two for GSF2124, GSF2056
    GSF_files = ['GSF2056-run1-plates1-2-demultiplexing-summary',
            'GSF2056-run2-plates3-4-demultiplexing-summary',
            'SampleSheet-GSF2124-run3-plates1-2',
            'GSF2124 Lennon Run 3 Plates 3-4 Run Summary Sorted',
            'GSF2124-run5-plates5-6-demultiplexing-summay']
    ignore_lines = ['Undetermined', 'Sample', 'Lane', 'Sample_ID', ' Chemistry',
                    'Description', 'Assay', 'Application', 'Workflow', 'Date',
                    'Experiment Name', 'IEMFileVersion', 'Lane Summary',
                    '"GSF2124 Lennon Plates 5-6', 'GSF2124-plates5-6-run5 Summary',
                    '', 'Chemistry', 'GSF2056-run2-plates3-4 Lennon Summary',
                    'GSF2056-run1-plates1-2 Lennon/Shoemaker Summary']

    GSF_bc_dict = {}
    for GSF_file in GSF_files:
        GSF_file_ = open(pt.get_path() + '/data/library_metadata/' + GSF_file + '.csv', 'r')
        for GSF_line in GSF_file_:
            GSF_line = GSF_line.strip()#.split(',')
            if len(GSF_line) < 20:
                continue
            GSF_line = GSF_line.split(',')
            if GSF_line[0] in ignore_lines:
                continue
            if GSF_line[2] == 'Undetermined':
                continue
            if GSF_line[0] == '1':
                GSF_line = GSF_line[2:]

            #GSF_bc_dict[GSF_line[0]] = {}
            #GSF_bc_dict[GSF_line[0]]['BC1'] = {}

            #split the barcodes


    meta_dict = {}
    meta_path = open(pt.get_path() + '/data/library_metadata/' + 'sample_names.txt', 'r')
    lengthsss =[]
    ressssss = []
    for line in meta_path:
        line = line.strip()
        line_dash = line.split('/')
        run = line_dash[0]

        if 'GSF' not in run:
            run = 'HCGS' + run
        if '_' in run:
            run = run.replace('_', '-')
        # line
        file_name = line_dash[-1]
        file_name_spl = re.split('-|_',file_name)

        if file_name_spl[0] == 'GSF2124':
            if len(file_name_spl) == 9:
                line = file_name_spl[4]
                day = file_name_spl[5][1:]
                R = file_name_spl[-2]

            elif len(file_name_spl) == 10:
                line = file_name_spl[3] + file_name_spl[4] + file_name_spl[5]
                day = file_name_spl[6]
                R = file_name_spl[-2]

            elif len(file_name_spl) == 11:
                line = file_name_spl[4] + file_name_spl[5] + file_name_spl[6]
                day = file_name_spl[7]
                R = file_name_spl[-2]

        elif file_name_spl[0] == 'GSF2056':
            if len(file_name_spl) == 13:
                line = file_name_spl[3] + file_name_spl[7] + file_name_spl[8]
                day = file_name_spl[9]
                R = file_name_spl[-2]
                end = file_name_spl[-1]

            elif len(file_name_spl) == 14:
                line = file_name_spl[4] + file_name_spl[8] + file_name_spl[9]
                day = file_name_spl[10]
                R = file_name_spl[-2]
                end = file_name_spl[-1]

        elif 'HCGS' in run:
            if len(file_name_spl) == 8:
                line = file_name_spl[0]
                day = file_name_spl[1]
                BC1 = file_name_spl[3]
                BC2 = file_name_spl[4]
                R = file_name_spl[-2]
                end = file_name_spl[-1]

            elif len(file_name_spl) == 9:
                line = file_name_spl[1] + file_name_spl[2]
                day = file_name_spl[3][1:]
                BC1 = file_name_spl[4]
                BC2 = file_name_spl[5]
                R = file_name_spl[-2]
                end = file_name_spl[-1]

            elif len(file_name_spl) == 6:
                line = file_name_spl[0][1:]
                day = '100'
                BC1 = file_name_spl[1]
                BC2 = file_name_spl[2]
                R = file_name_spl[4]
                end = file_name_spl[-1]

            elif len(file_name_spl) == 7:
                line = file_name_spl[0]
                day = file_name_spl[1]
                BC1 = file_name_spl[2]
                BC2 = file_name_spl[3]
                R = file_name_spl[5]
                end = file_name_spl[-1]

        #ressssss.append(run)


        #meta_dict[line] = [run, line, day, BC1, BC2, R, rep]
    #print(set(lengthsss))
    #print(set(ressssss))

    #print(meta_dict)


    #df_paths['A'] = df_paths['file_path'].str.split('/', 1).str
    #df['First_Name'] = df.file_path.str.split('/', expand = True)[0]


    #run = runs[0]
    #if 'SampleSheets_HCGB_' in run:
    #    df = pd.read_csv(pt.get_path() + '/data/libraries/' + run + '.csv', sep = ',')
    #    df = df.iloc[::2, :]
    #    df_new = df["Index"].str.split("-", n = 1, expand = True)
    #    df["Index1"]= df_new[0]
    #    df["Index2"]= df_new[1]
    #    #df.drop(columns =["Index"], inplace = True)
    #    #out_name = pt.get_path() + '/data/libraries/' + run + '_clean.csv'
    #    #df.to_csv(out_name, sep = ',', index = True)
    #    if '161213' in run:
    #        df["New_sample_name"] = df['SampleID'] + '-' 'D100'
    #    print(df)

    #df = pd.read_csv(pt.get_path() + '/data/libraries/' + run + '.csv', sep = ',')

merge_metadata()

# create summary file of all data

#clean_hcgs_sheets()
