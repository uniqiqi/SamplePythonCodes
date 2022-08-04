import os
import pandas as pd
import csv
import argparse
import subprocess
import json
import numpy as np

def get_arguments():

    """
    input_file
    metagenome samples #to be analyzed # input different format
    antiSAMSH #antiSMASH sample or not T-> calc both whole and core F-> calc only whole
    cov_cutoff # coverage cutoff for whole BGC
    corecov_cutoff # core coverage cutoff for BGC
    Output_file
    """
    parser = argparse.ArgumentParser(description="",
                                    usage=''' 
    *** process the result.ALL file from BiG-MAP.map.py ***
    *** filter BGCs with low coverage score and re-arrange the data ***
    -B1 biom file for whole gene cluster mapping (required)
    -O output folder name (required): output folder name
    -B2 biom file for core gene cluster mapping (optional): include only when you want to include a filtration based on core coverage
    -c1 cov_cutoff (optional): default= 0.5
    -c2 corecov_cutoff (optional): default= 0.5
    
''')

    parser.add_argument("-B1", "--biom", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-B2", "--core_biom", help=argparse.SUPPRESS, required=False)
    parser.add_argument("-O", "--output_folder", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-c1", "--cov_cutoff", type=float, help=argparse.SUPPRESS, required=False, default=0.5)
    parser.add_argument("-c2", "--corecov_cutoff", type=float, help=argparse.SUPPRESS, required=False, default=0.5)
    return (parser.parse_args())



def convert_biom(biom, output_folder, core="whole"):
    # if core="whole", return json file named "whole_map_result", input file should be args.biom
    # if core="core", return json file named "core_map_result", input file should be args.core_biom
    json_file = os.path.join(output_folder, core + "_map_result.txt")
    cmd_convert_biom = f"biom convert -i {biom} -o {json_file} --table-type='Pathway table' --to-json"
    res_export = subprocess.check_output(cmd_convert_biom, shell=True)
    return(json_file)




def extract_biom(biom_dict):
    values = []
    gc_ids = []
    df_gc_ids = []
    sample_ids = []
    
    for record in biom_dict["rows"]:
        gc_ids.append(record["id"])
        regionno = record["id"].split('|')[1]
        product = record["id"].split('=')[1].split('--')[0]
        sourceo = record["id"].split('=')[2].split('--')[0]
        redundency = record["id"].rsplit('--')[-1]
        df_gc_ids.append(regionno + "|" + product + "|" + sourceo + "|" +redundency)
    for record in biom_dict["columns"]:
        sample_ids.append(record["id"])
    # generate a blank dictinary (#gc X #sample)    
    for id in gc_ids:
        dict = {'id': id,'values':{}}
        for sample in sample_ids:
            dict['values'][sample] = 0
        values.append(dict)
    
    # insert the nonzero data from biom file into the dictionary, otherwise numbers will be kept for zero
    data = biom_dict["data"]    #[rows,clolumns,values]
    for i in data:
        gc = gc_ids[i[0]]
        sample = sample_ids[i[1]]
        for record in values:
            if record["id"] == gc:
                record["values"][sample] = i[2]
    
    return(values, gc_ids, sample_ids, df_gc_ids)



def restructure_data(values, cov, corecov, sample_ids, gc_ids):
    restr_data = {}
    for sample in sample_ids:
        restr_data[sample] = {}
        for gc in gc_ids:
            restr_data[sample][gc] = [0,0,0]
    # insert coverage values into the restructured data
    for record in cov:
        gc_id = record["id"]
        for sample in record["metadata"].keys():
            #print(sample)
            #print(record["metadata"][sample])
            restr_data[sample][gc_id][0] = float(record["metadata"][sample])
    # by the same way, insert core coverage values into the restructured data
    for record in corecov:
        gc_id = record["id"]
        for sample in record["metadata"].keys():
            restr_data[sample][gc_id][1] = float(record["metadata"][sample])
    for record in values:
        gc_id = record["id"]
        for sample in record["values"].keys():
            restr_data[sample][gc_id][2] = float(record["values"][sample])
            
    return(restr_data)




def filter_coverage(restr_data, cutoff, core_cutoff):
    filtered_data = {}
    filter_notes = []
    filtered_data = restr_data
    for i in filtered_data.keys():
        for j in filtered_data[i].keys():
            if filtered_data[i][j][0] < cutoff:
                note = i + " - " + j + ": filter based on whole BGC coverage: " + str(filtered_data[i][j][0])
                filter_notes.append(note)
                filtered_data[i][j][0] = 0
                filtered_data[i][j][1] = 0
                filtered_data[i][j][2] = 0
            elif filtered_data[i][j][1] < core_cutoff:
                note = i + " - " + j + ": filter based on core BGC coverage: " + str(filtered_data[i][j][1])
                filter_notes.append(note)
                filtered_data[i][j][0] = 0
                filtered_data[i][j][1] = 0
                filtered_data[i][j][2] = 0
            else:
                note = i + " - " + j + ": RPKM value taken for heatmap based on whole coverage and core coverage: " \
                + str(filtered_data[i][j][0]) + " " + str(filtered_data[i][j][1])
                filter_notes.append(note)
    return(filtered_data, filter_notes)



def make_dataframe(dict, gc_ids):
    dataframe = pd.DataFrame(index=gc_ids)
    for i in dict.keys():
        colnames = [i + ".cov", i + ".corecov",i + ".RPKM"]
        data = pd.DataFrame.from_dict(dict[i], orient='index', columns= colnames)
        dataframe = pd.concat([dataframe, data], axis=1)  
    return(dataframe)




def group_samples(biom_dict, sample_ids):
    group_dict = {}
    sample_data = biom_dict["columns"]
    #print(sample_data)
    for record in sample_data:
        sample_id = record["id"]
        group = record["metadata"]["Isolation"]
        if group not in group_dict.keys():
            group_dict[group] = []
        group_dict[group].append(sample_id)
    return(group_dict)




def calculate_mean(data):
    mean_calc = data.replace(0, np.nan)
    mean = mean_calc.mean(axis=1)
    return(mean)
    
def calculate_perc(data):
    bin_calc = np.sign(data)
    perc = bin_calc.mean(axis=1)
    return(perc)




def norm_log2_data(df):
    """
    This function is cirectly cited from BiG-MAP.analyse.py 
    to conduct RPKM value normalization based on the filtered data
    ----------
    Normalization and log2 convertion as performed by metagenomeseq
    ----------
    df
        dataframe, rows with samples, columns with GC and RPKM values
    returns
    ----------
    norm_df = dataframe with normalized and log2 converted RPKM values
    """
    norm_df = pd.DataFrame()

    df_adj = df.replace(0, np.nan)
    quantile_dict = (df_adj.quantile(axis=0)).to_dict() # calculate the quantile
    # determine the numeric values in the dataframe
    numeric_cols = [col for col in df_adj if df_adj[col].dtype.kind != 'O']
    df_adj[numeric_cols] -= np.finfo(float).eps # substract the machine epsilon
    # calculate the normalization factor by taking the sum of the values below the quantile \
    # normalize the data by dividing the counts by the normalization factor
    for sample in quantile_dict.keys():
        normfac = df_adj[sample][df_adj[sample] <= quantile_dict[sample]].sum()
        norm_df = norm_df.append(df_adj[sample]/(normfac/1000))
    norm_df = norm_df.T
    norm_df[numeric_cols] += 1 # add 1 as correction for log2 calc
    norm_df = (np.log2(norm_df)).replace(np.nan, 0.0) # calculate log2
    return norm_df




def main():
    # get arguments
    args = get_arguments()
    print(" ")
    # clarify arguments
    if args.core_biom:
        print("processing data from " + args.biom + " and " + args.core_biom)
        print("whole BGC coverage cutoff=" + str(args.cov_cutoff) + ", core BGC cutoff= " + str(args.corecov_cutoff))
        print("whole BGC coverage cutoff=" + str(args.cov_cutoff) + ", core BGC cutoff= " + str(args.corecov_cutoff))
        print("All files output to " + args.output_folder)
        cov_cutoff = args.cov_cutoff
        corecov_cutoff = args.corecov_cutoff
    else:
        print("processing data from " + args.biom)
        print("will not compute coverage of core biosynthetic genes")
        print("whole BGC coverage cutoff=" + str(args.cov_cutoff))
        print("All files output to " + str(args.output_folder))
        cov_cutoff = args.cov_cutoff
        corecov_cutoff = 0
    
    print(" ")
    print("==============================================================")
    print(" ")
    # test output
    print("whole BGC cutoff: " + str(cov_cutoff))
    print("core BGC cutoff: " + str(corecov_cutoff))

    # convert biom file to dictionary
    # parse data from biom file
    whole_json = convert_biom(args.biom, args.output_folder, "whole")
    with open (whole_json, "r") as file:
        whole_dict = json.load(file)
        
    values, gc_ids, sample_ids, df_gc_ids = extract_biom(whole_dict)
    cov = whole_dict["rows"]
    
    if args.core_biom:
        core_json = convert_biom(args.core_biom, args.output_folder, "core")
        with open (core_json, "r") as file:
            core_dict = json.load(file)
        corecov = core_dict["rows"]
    else:
        # if core_biom not provided in the arguments, use whole cov for parsing to maintain the structure
        corecov = whole_dict["rows"]
        
    # restructure RPKM, cov and corecov to a dict with defined index and column names
    restr_data = restructure_data(values, cov, corecov, sample_ids, gc_ids)
    filtered_data, filtered_notes = filter_coverage(restr_data, cov_cutoff, corecov_cutoff)
    # data output
    with open(args.output_folder + "/filter_notes.txt", "w") as file:
        for line in filtered_notes:
            file.write("%s\n" % line)

    # covert restructured dict data to dataframe
    original_df = make_dataframe(restr_data, gc_ids)
    original_df.index = df_gc_ids
    filtered_df = make_dataframe(filtered_data, gc_ids)
    filtered_df.index = df_gc_ids
    # data output
    original_df.to_csv(args.output_folder + "/orginal_data.csv")
    filtered_df.to_csv(args.output_folder + "/filtered_data.csv")
    print(" ")
    print("=================================================")

    # analyzing meta sample groups by reading the meta data from biom file
    sample_groups = group_samples(whole_dict, sample_ids)
    for group in sample_groups.keys():
        print(group + ": ")
        print(sample_groups[group])
        # output sample groups in the present data

    # output processed data by metasample groups which can be directly use for downstreamm heatmap drawing
    # data contains narmolized RPKM value, 
    # and average whole/core coverage values across all non-zero metasamples in this group
    NORMED_DATA = pd.DataFrame(index = filtered_df.index.values)
    for group_name in sample_groups.keys():
        g_r_mean = group_name + ".RPKM_mean"
        g_perc =  group_name + ".perc"
        g_cov_mean = group_name + ".cov_mean"
        g_corecov_mean = group_name+ ".corecov_mean"
        unnm_outfile = group_name + "_RPKM_data.csv"
        nm_outfile = group_name + "_heatmap_data.csv"
        print("Now processing data from [ " + group_name + " ]")
        print("The summerized data will be output to: " + unnm_outfile)
        print("The normalized data will be output to: " + nm_outfile)
        print(" ")
        value_cols = []
        cov_cols = []
        corecov_cols = []
        for sample in sample_groups[group_name]:
            value_cols.append(sample + ".RPKM")
            cov_cols.append(sample + ".cov")
            corecov_cols.append(sample + ".corecov")
        values_data = pd.DataFrame(filtered_df[[col for col in filtered_df if col in value_cols]])
        normed_data = norm_log2_data(values_data)
        cov_data = pd.DataFrame(filtered_df[[col for col in filtered_df if col in cov_cols]])
        corecov_data = pd.DataFrame(filtered_df[[col for col in filtered_df if col in corecov_cols]])
        
        ### ADD RPKM_mean to the dataframe for sorting
        normed_mean = calculate_mean(normed_data)
        normed_data[g_r_mean] = normed_mean
        values_data[g_r_mean] = calculate_mean(values_data)
        perc = calculate_perc(values_data)
        normed_data[g_perc] = perc
        values_data[g_perc] = perc
        cov_mean = calculate_mean(cov_data)
        normed_data[g_cov_mean] = cov_mean
        values_data[g_cov_mean] = cov_mean
        corecov_mean = calculate_mean(corecov_data)
        normed_data[g_corecov_mean] = corecov_mean
        values_data[g_corecov_mean] = corecov_mean
        normed_data = normed_data.sort_values(by=[g_perc, g_r_mean, g_cov_mean, g_corecov_mean], ascending = False)
        values_data = values_data.sort_values(by=[g_perc, g_r_mean, g_cov_mean, g_corecov_mean], ascending = False)
        normed_data.to_csv(args.output_folder + "/" + nm_outfile)
        values_data.to_csv(args.output_folder + "/" + unnm_outfile)
        NORMED_DATA = NORMED_DATA.merge(normed_data, left_index=True, right_index=True, copy = False)
        NORMED_DATA.to_csv(args.output_folder + "/compiled_data.csv")
        print(" ")
        

if __name__ == "__main__":
    main()

