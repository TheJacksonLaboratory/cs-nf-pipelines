#!/usr/bin/env python

################################################################################
##################### COPYRIGHT ################################################
# New York Genome Center

# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2018) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.

# Version: 0.1 (2018-12-06)
# Author: Kanika Arora (karora@nygenome.org)
##################### COPYRIGHT ################################################
################################################################################


import sys
import argparse
import os
import re

class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, '\nERROR: %s\n\n' % (message))


def check_if_exists(file_or_dir_path, type="file"):
    if type=="file":
        if not os.path.isfile(file_or_dir_path):
            print("ERROR: Required file {0} does not exist. Cannot run.".format(file_or_dir_path))
            sys.exit(1)
    elif type=="directory":
        if not os.path.isdir(file_or_dir_path):
            print("ERROR: Required file {0} does not exist. Cannot run.".format(file_or_dir_path))
            sys.exit(1)
    else:
        if not os.path.exists(file_or_dir_path):
            print("ERROR: {0} does not exist. Cannot run.".format(file_or_dir_path))
            sys.exit(1)

def compute_vaf(alt_count,dp):
    '''
        Compute VAF from dp and alt_count.
    '''
    if not isinstance(alt_count, int):
        if str(int(alt_count)) == alt_count:
            alt_count=int(alt_count)
        else:
            raise ValueError("alt_count should be an integer. "+alt_count+" provided. Cannot run")
    if not isinstance(dp, int):
        if str(int(dp)) == dp:
            dp=int(dp)
        else:
            raise ValueError("dp should be an integer. "+dp+" provided. Cannot run")
    if alt_count > dp:
        raise ValueError("alt_count {0} is greater than depth {1}.".format(alt_count,dp))
    return (0 if dp==0 else round(float(alt_count)/dp,4))

def parse_format_return_allele_counts(ref,alt,format_dict,caller):
    '''
        Return allele counts, depth and allele fraction as reported by the given caller.
        If these fields are not present as is, they are computed from other fields from that caller.
    '''
    AD = format_dict[caller+"_AD"] if caller+"_AD" in format_dict else "."
    DP = format_dict[caller+"_DP"] if caller+"_DP" in format_dict else "."
    AF = format_dict[caller+"_AF"] if caller+"_AF" in format_dict else "."
    if caller == "lancet":
        AF=str(compute_vaf(AD.split(",")[1],DP))
    if caller == "strelka2":
        ## compute as suggested in https://github.com/Illumina/strelka/blob/master/docs/userGuide/README.md#somatic
        if len(ref)==1 and len(alt)==1:
            ## SNV
            AD_ref=format_dict["strelka2_"+ref+"U"].split(",")[0]
            AD_alt=format_dict["strelka2_"+alt+"U"].split(",")[0]
            DP=str(int(format_dict["strelka2_AU"].split(",")[0])+int(format_dict["strelka2_CU"].split(",")[0])+int(format_dict["strelka2_GU"].split(",")[0])+int(format_dict["strelka2_TU"].split(",")[0]))
        else:
            ## INDEL
            AD_ref=format_dict["strelka2_TAR"].split(",")[0]
            AD_alt=format_dict["strelka2_TIR"].split(",")[0]
            DP=str(int(AD_ref)+int(AD_alt))
        AD=AD_ref+","+AD_alt
        AF=str(compute_vaf(AD_alt,DP))
    return (AD,DP,AF)



def __main__():
    parser = ArgumentParser(prog='add_final_allele_counts',
                            description='Picks final values for AD and DP based on set caller priority.', epilog='',
                            formatter_class=lambda prog: argparse.ArgumentDefaultsHelpFormatter(prog, max_help_position=100, width=150))
    parser.add_argument('-v', '--vcf', help = 'SNV VCF file.', required=True)
    parser.add_argument('-o', '--output', help = 'Output VCF file.', required=True)
    parser.add_argument('-p', '--priority', help = 'Comma-separated prioritized list of sources (callers) for picking final allele counts for the variants.', default='nygc,strelka2,mutect2,lancet')
    args=parser.parse_args()
    VCF=args.vcf
    OUT=args.output
    check_if_exists(VCF)
    CALLER_PRIORITY=args.priority.split(",")
    f=open(VCF)
    o=open(OUT,"w")
    seen=0
    for line in f:
        if line.startswith("#"):
            o.write(line)
            if line.startswith("##FORMAT") and seen==0:
                o.write('##INFO=<ID=AlleleCountSource,Number=1,Type=String,Description="Name of caller or method that was used for the AD,DP and AF fields">\n')
                o.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles.">\n')
                o.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth of coverage: Number of (filtered) reads covering site.">\n')
                o.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction of alternate allele in the sample.">\n')
                seen=1
        else:
            toks=line.strip().split("\t")
            info_pattern=r"(\w+)=([^;]+);*"
            matches = re.findall(info_pattern,toks[7])
            info_dict = {a:b for a,b in matches}
            normal_format=dict(zip(toks[8].split(":"),toks[9].split(":")))
            tumor_format=dict(zip(toks[8].split(":"),toks[10].split(":")))
            if 'called_by' in info_dict:
                called_by=info_dict["called_by"].split(",")
            else:
                called_by = []
            if "nygc_AD" in normal_format:
                called_by.append("nygc")
            chosen_caller=""
            for caller in CALLER_PRIORITY:
                if caller in called_by:
                    chosen_caller=caller
                    break
            if chosen_caller=="":
                o.write(line)
            else:
                toks[7]=toks[7]+";AlleleCountSource="+chosen_caller
                (normal_AD,normal_DP,normal_AF)=parse_format_return_allele_counts(toks[3],toks[4],normal_format,chosen_caller)
                (tumor_AD,tumor_DP,tumor_AF)=parse_format_return_allele_counts(toks[3],toks[4],tumor_format,chosen_caller)
                toks[8]=toks[8]+":AD:DP:AF"
                toks[9]=toks[9]+':{0}:{1}:{2}'.format(normal_AD,normal_DP,normal_AF)
                toks[10]=toks[10]+':{0}:{1}:{2}'.format(tumor_AD,tumor_DP,tumor_AF)
                o.write("\t".join(toks)+"\n")
    f.close()
    o.close()

if __name__ == "__main__":
    __main__()