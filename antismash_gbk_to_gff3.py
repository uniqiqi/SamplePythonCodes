### the script is used to convert antiSMASH generated gene cluster genbank file to GFF3 format ###
### for downstream roary and scoary analysis ###
### only CDS feature in the genbank file is considered given the ultimate goal of this project ###


"""
###############
gff file format
###############

##gff-version 3
##sequence-region [locus_tag_prefix1] [cluster_start1] [cluster_end1]
##sequence-region [locus_tag_prefix1] [cluster_start1] [cluster_end1]




##column info
[locus_tag_prefix1]	[source]	[feature_type]	[feature_start]	[feature_end]	[score]	[strand]	[phase]	[attributes]
[locus_tag_prefix1]	[source]	[feature_type]	[feature_start]	[feature_end]	[score]	[strand]	[phase]	[attributes]
[locus_tag_prefix1]	[source]	[feature_type]	[feature_start]	[feature_end]	[score]	[strand]	[phase]	[attributes]
...
...
...
[locus_tag_prefix2]	[source]	[feature_type]	[feature_start]	[feature_end]	[score]	[strand]	[phase]	[attributes]
[locus_tag_prefix2]	[source]	[feature_type]	[feature_start]	[feature_end]	[score]	[strand]	[phase]	[attributes]
[locus_tag_prefix2]	[source]	[feature_type]	[feature_start]	[feature_end]	[score]	[strand]	[phase]	[attributes]

...
##FASTA
>[locus_tag_prefix]
whole fasta sequence record of 1
>[locus_tag_prefix]
whole fasta sequence record of 2
"""

import os
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import fnmatch
import argparse


def get_arguments():
	parser = argparse.ArgumentParser(description="",
									usage=''' 
    *** convert antiSMASH generated gbk file to gff format ***
    -F folder that contains the antiSMASH result from a bacetria genome
    -O output foldername
''')

	parser.add_argument("-F", "--folder", help=argparse.SUPPRESS, required=True)
	parser.add_argument("-O", "--outdir", help=argparse.SUPPRESS, required=True)

	return (parser.parse_args())

def wrapped_fasta(string, every=100):
	"""
	truncate fastan sequence for every 100 nucleotdes
	"""
	wrapped_fasta = "\n".join(string[i:i+every] for i in range(0, len(string), every)) + "\n"
	return(wrapped_fasta)

def parse_feature(feature):
	"""
	extact information for every CDS feature in the antiSAMSH generated genbank file
	---------
	input: CDS feature from an antiSMASH BGC genbank file
	---------
	output:
		features_info: a dictionary containing information of each CDS feature that gff reuiqred 
	"""
	feature_info = {}
	### exclude truncated protein/genes on the contig edge
	start = str(int(feature.location.start)+1)
	feature_info["start"] = start

	end = str(int(feature.location.end))
	feature_info["end"] = end

	strand = "+" if feature.strand == 1 else "-"
	feature_info["strand"] = strand

	if "locus_tag" in feature.qualifiers:
		item2 = "locus_tag=" + ''.join(feature.qualifiers['locus_tag'])
	else:
		item2 = "locus_tag=" + item1
	feature_info["locus_tag"] = item2

	if "product" in feature.qualifiers:
		item3 = "product=" + ''.join(feature.qualifiers['product'])
	else:
		item3 = "product=" + "no_prediction"
	feature_info["product"] = item3

	if "protein_id" in feature.qualifiers:
		item4 = "protein_id=" + ''.join(feature.qualifiers['protein_id'])
	else:
		item4 = "protein_id=" + "no_protein_id"
	feature_info["protein_id"] = item4

	if "gene_kind" in feature.qualifiers:
		gene_kind = ''.join(feature.qualifiers['gene_kind'])
		if gene_kind == "biosynthetic" or gene_kind == "biosynthetic-additional":
			item5 = "gene_kind=" + gene_kind
		else:
			item5 = "gene_kind=" + "other_genes"
	else:
		item5 = "gene_kind=" + "other_genes"
	feature_info["gene_kind"] = item5

	return(feature_info)


def get_genome_name(string):
	if "/" in string:
		if string.endswith('/') == True:
			genome_name = string.rsplit('/')[-2]
		else:
			genome_name = string.rsplit('/')[-1]
	else:
		genome_name = string

	return(genome_name)


def main():

	args = get_arguments()
	genome_name = get_genome_name(args.folder)

	dir_in_str = args.folder
	outputfile = args.outdir + "/" + genome_name + ".gff"
	print("output .gff file to " + outputfile)
	### write first line
	### divide each gff file into three part
	### append info from each BGC file to each part after reading in the genbank file
	### when writing output, write the three part in order to a final gff file for each genome

	title = []
	title.append("##gff-version 3\n")
	info = []
	fasta_sequence = []
	region_record = 0

	directory = os.fsencode(dir_in_str)
	for file in os.listdir(directory):
		filename = os.fsdecode(file)
		filepath = os.path.join(dir_in_str, filename)
		if fnmatch.fnmatch(filename, '*.region*.gbk'):
			region_record += 1
			### generate accession prefix first ###
			prefix = filename.rsplit('.', 1)[0]
			if "..." in prefix:
				prefix = prefix.split('...',1)[0] + "." + prefix.split('...',1)[1]
			### generate FASTA header based on prefix
			fasta_header = ">" + prefix + "\n"
			### read in the record ####
			record = SeqIO.read(filepath,"genbank")
			### get whole length of the record ###
			cluster_length = str(len(record.seq))
			### get the whole sequence of the record ###
			sequence = str(record.seq)
			title.append("##sequence-region " + prefix + " 1 " + cluster_length + "\n")
			fasta_sequence.append(fasta_header)
			fasta_sequence.append(wrapped_fasta(sequence))
			fea_num = 0
			for feature in record.features:
				if feature.type == "CDS":
					fea_num = fea_num + 1
					feat = parse_feature(feature)
					id = "ID=" + prefix + "_" + str(fea_num).zfill(5)
					info_line = prefix  + "\tgenbank\tCDS\t" + feat["start"] + "\t" + feat["end"] + "\t.\t" + feat["strand"] + "\t0\t"
					attributes = id + ";" + feat["locus_tag"] + ";" + feat["product"] + ";" + feat["protein_id"] + ";" + feat["gene_kind"]
					info_line = info_line + attributes + "\n"
					info.append(info_line)
	if region_record != 0:
		###output to file ###
		with open(outputfile, "w") as output:
			for cluster in title:
				output.write(cluster)
			for feature in info:
				output.write(feature)
			for seq in fasta_sequence:
				output.write(seq)

		output.close()


if __name__ == "__main__":
	main()