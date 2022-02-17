#!/usr/bin/env python

# aggregate_CGA_stats.py OUT INP_QC INP_HS
# A script to parse the quality and hybrid-selection statistics files
# and aggregate relevant metrics into an output file.

# Parameters:
# out = output file
# inp_qc = *stat file output by NGSQCToolkit
# inp_hs = *Metricsfile.txt outoput by Picard CalculateHsMetrics 

import sys

if len(sys.argv) < 4:
    print >>sys.stderr, "Commandline arguments missing:\nFormat: aggregate_CGA_stats.py OUT INP_QC INP_HS\nout = output file\ninp_qc = *stat file output by NGSQCToolkit\ninp_hs = *Metricsfile.txt outoput by Picard CalculateHsMetrics"
    sys.exit()
    
out = open(sys.argv[1],"w")
inp_qc = open(sys.argv[2],"r")
inp_dup = open(sys.argv[3],"r")
inp_hs = open(sys.argv[4],"r")

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

data_lines = []
# for line in inp_hs:
#    print(line)
#    line = line.strip()
#    print(line)
#    if line and not(line.startswith("#")):
#        data_lines.append(line)

with open(sys.argv[4], "r") as inp_hs:
    for line in inp_hs:
        #print(line.strip())
        if line.startswith("## METRICS CLASS"):
            #print(line)
            data_lines.append(next(inp_hs, '').strip())
            data_lines.append(next(inp_hs, '').strip())
            #line.strip())
            #data_lines.append()
            #

print(data_lines)

#print(len(data_lines))

if len(data_lines) != 2:
    print >>sys.stderr, "CoverageMetrics.txt is invalid"
else:
    col_names = data_lines[0].split("\t")
    col_values = data_lines[1].split("\t")
    for i,n in enumerate(col_names):
        if n in ["PF_UNIQUE_READS", "PCT_PF_UQ_READS_ALIGNED", "PCT_SELECTED_BASES", "MEAN_TARGET_COVERAGE"] or n.startswith("PCT_TARGET_BASES"):
            print >>out, "%s\t%s" %(n,col_values[i])


