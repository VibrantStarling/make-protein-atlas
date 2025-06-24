import os
import sys



def extract_meta_data(inputfile,outprefix):
    cwd = os.getcwd()
    metadata_files = []
    raw_file_spectra_counts = {}
    raw_file_psm_counts = {}

    for f in os.listdir(cwd):
        if os.path.isfile(f) and "metadata.txt" in f:
            metadata_files.append(f)

    if len(metadata_files) == 0:
        exit(
            "No metadata files found, you need to run ThermoRawFileParser -i $file -m 1 on all raw files first to generate metadata text files")

    for file in metadata_files:
        f_m = open(file,"r")
        #print("file",file)

        raw_file_name = "error_not_found"
        for line in f_m:
            line = line.rstrip()

            if "RAW file path" in line:
                raw_file_name = line.split("=")[1]

                if raw_file_name[0:2] == "./":
                    raw_file_name = raw_file_name[2:] #ThermoFileReader has a small bug, in that if you use a whole directory mode it puts a prefix on file name

                #print("Raw",raw_file_name)
            elif "Number of MS2 spectra" in line:
                raw_file_spectra_counts[raw_file_name] = line.split("=")[1]
        f_m.close()

    #Now read the PSM tsv file

    if os.path.isfile(inputfile):
        f_psms = open(inputfile,"r")
        line_counter = 0
        for line in f_psms:
            line = line.rstrip()
            cells = line.split("\t")
            #header row
            if line_counter == 0:
                if cells[0] != "spectrum":
                    exit("Expecting first column to have header spectrum, exiting\n") #TODO recognise file headers
            else:
                spectrum = cells[0]
                raw_file_name = spectrum.split(".",1)[0]+".raw" #150923_KT_YiWa_10-1.00004.00004.2
                psm_count_per_raw_file = 0
                if raw_file_name not in raw_file_spectra_counts:
                    print("Didn't recognise raw file name in dictionary, exiting:",raw_file_name)
                    exit()
                if raw_file_name in raw_file_psm_counts:
                    psm_count_per_raw_file = raw_file_psm_counts[raw_file_name]
                psm_count_per_raw_file +=1
                raw_file_psm_counts[raw_file_name] = psm_count_per_raw_file
            line_counter+=1
    else:
        print("file not found:",inputfile)
        exit()

    f_out = open(outprefix + "_psm_stats.tsv","w")
    f_out.write("raw_file\tpsm_count\tspectra_count\trecovery\n")
    for raw_file in raw_file_spectra_counts.keys():

        if raw_file in raw_file_spectra_counts:
            spectra_count = raw_file_spectra_counts[raw_file]
            if raw_file in raw_file_psm_counts:
                psm_count = raw_file_psm_counts[raw_file]
                f_out.write(raw_file + "\t" + str(spectra_count) + "\t" + str(psm_count) + "\t" + str(float(psm_count)/float(spectra_count))+"\n")
    f_out.close()

    #for file in raw_file_psm_counts.keys():
    #    print(file,"\t",raw_file_spectra_counts[file])


if len(sys.argv)!= 3:
    print("Exit - expected usage\npython ",  sys.argv[0]," thresholded_psm_file.tsv file_prefix\n"
          "This script is to be run in a folder already containing one metadata.txt file from ThermoRawFile Parser per raw file\n"
                                                         ", and the thresholded PSM file in tsv format\n"
                                                         "file_prefix should be PXD code to give file name to write stats to ")
else:
    extract_meta_data(sys.argv[1],sys.argv[2])









