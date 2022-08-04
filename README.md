# SamplePythonCodes
Two example python codes from the contributor. More systematic analysis codes will be available once the research paper is published. 


antismash_gbk_to_gff3.py -- 
The script is used to convert antiSMASH generated gene cluster genbank file to GFF3 format for downstream roary and scoary analysis. 
Only CDS feature in the genbank file is considered given the ultimate goal of this project.

data_processing_for_heatmap.py -- 
Conver Biom file to three dataframes containing whole BGC coverage score, core BGC coverage score and a normalized RPKM value for heatmap plotting. 
