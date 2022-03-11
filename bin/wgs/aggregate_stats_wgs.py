#!/usr/bin/env python

# aggregate_stats.py OUT INP_QC INP_DUP INP_ALIGNMENT_STAT
# A script to parse the quality and statistics files
# and aggregate relevant metrics into an output file.

# Parameters:
# out = output file
# inp_qc = *stat file output by qualtool
# inp_dup = *.dat output picard mark duplicates
# inp_hs = *Metricsfile.txt outoput by Picard Alignment

import sys

if len(sys.argv) < 5:
    print >>sys.stderr, "Commandline arguments missing:\nFormat: aggregate_stats.py OUT INP_QC INP_DUP INP_ALIGNMENT_STAT\nout = output file\ninp_qc = *stat file output by qualtool\ninp_dup = *.dat output by Picard MarkDuplicates\ninp_hs = *AlignmentMetrics.txt outoput by Picard CollectAlignmentSummaryMetrics"
    sys.exit()
    
out = open(sys.argv[1],"w")
inp_qc  = open(sys.argv[2],"r")
inp_dup = open(sys.argv[3],"r")
#inp_hs  = open(sys.argv[4],"r")

qc_out = [None, None]
read_data = False
for line in inp_qc:
    line = line.strip()
    elems = line.split("\t")
    
    if line.startswith("QC statistics"):
        read_data = True
    
    if line.startswith("Detailed QC statistics"):
        break
    
    if read_data:
        if None not in qc_out:
            break
        if line.startswith("Total number of reads"):
            try:
                elems = line.split() 
                qc_out[0] = str(int(elems[-1]) + int(elems[-2]))
            except Exception:
                qc_out[0] = "NA"
        if line.startswith("Total number of HQ filtered reads"):
            try:
                elems = line.split() 
                qc_out[1] = str(int(elems[-1]) + int(elems[-2]))
            except Exception:
                qc_out[1] = "NA"
print >>out, "Total number of reads\t%s\nTotal number of HQ filtered reads\t%s" %(qc_out[0],qc_out[1])

####

data_lines_dup = []
for line in inp_dup:
    line = line.strip()
    if line and not(line.startswith("#")):
        data_lines_dup.append(line)

col_names = data_lines_dup[0].split("\t")
col_values = data_lines_dup[1].split("\t")
for i,n in enumerate(col_names):
        if n in ["PERCENT_DUPLICATION"]:
            print >>out, "%s\t%s" %(n,col_values[i])

####

data_lines = []
with open(sys.argv[4], "r") as inp_hs:
    for line in inp_hs:
        if line.startswith("## METRICS CLASS"):
            data_lines.append(next(inp_hs, '').strip())
            data_lines.append(next(inp_hs, '').strip())
            data_lines.append(next(inp_hs, '').strip())
            data_lines.append(next(inp_hs, '').strip())

if len(data_lines) != 4:
    print >>sys.stderr, "AlignmentMetrics.txt is invalid"
else:
    col_names = data_lines[0].split("\t")
    col_values = data_lines[3].split("\t")    
    for i,n in enumerate(col_names):
        if n in  ["PF_READS_ALIGNED", "PCT_PF_READS_ALIGNED"]:
            print >>out, "%s\t%s" %(n,col_values[i])

    #data_lines[3] is for paired data. Could be an issue here with SE data. This requires testing. 

####

data_lines = []
with open(sys.argv[5], "r") as inp_cov:
    for line in inp_cov:
        if line.startswith("## METRICS CLASS"):
            data_lines.append(next(inp_cov, '').strip())
            data_lines.append(next(inp_cov, '').strip())


#####
if len(data_lines) != 2:
    print >>sys.stderr, "WGS_Metrics.txt is invalid"
else:
    col_names = data_lines[0].split("\t")
    col_values = data_lines[1].split("\t")    
    for i,n in enumerate(col_names):
        if n in  ["MEAN_COVERAGE", "SD_COVERAGE", "MEDIAN_COVERAGE"] or n.endswith("X"):
            print >>out, "%s\t%s" %(n,col_values[i])

