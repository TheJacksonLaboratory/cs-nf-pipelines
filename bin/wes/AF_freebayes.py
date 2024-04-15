#!/usr/bin/env python3
"""

GOALS of script:
    recompute the locus depth from the allele-depths

    Add Estimated Allele Frequency (ALT_AF) to the info cell

The FORMAT column indicates the order of the fields in the following,
sample, column.

There are several columns that need to be found:

    AD: Allele Depth. Comma delim with values for REF and each ALT. Required to adjust multi-allelic printout to one per line.
    AD: Allele Depth. Comma delim with values for REF and each ALT. Required to adjust multi-allelic printout to one per line.
    AC: Total number of alternate alleles in called genotypes. Required to adjust multi-allelic printout to one per line.
    RO: Reference allele observation count
    QA: Sum of quality of the alternate observations. Comma delim with values each ALT.
    GL: Genotype Likelihood. Comma delim with 3 values per allele. Required to adjust multi-allelic printout to one per line.

Note that we use AO and RO here as they "only include reads deemed good enough for an allele call."

For each ALT allele: AF = AO_allele / ( sum(AO_alleles) + RO ) 

Total depth is sum(AO_alleles) + RO

"""
import sys

# expecting sys.argv to look something like
if len(sys.argv) != 4:
    print("Incorrect number of args!")
    print("expected usage:")
    print("allele_depth_min_and_AF_from_ADs.py inputFile outputFile sample_column_number")
    raise Exception("Incorrect number of args!")

ALT_ALLELE_INDEX = 4
FILTER_CELL_INDEX = 6
INFO_CELL_INDEX = 7
FORMAT_CELL_INDEX = 8
SAMPLE_DATA_INDEX = int(sys.argv[3])

inp = open(sys.argv[1], 'r')
out = open(sys.argv[2], 'w')

NEW_INFO_HEADERS = ['##INFO=<ID=ALT_AF,Number=A,Type=Float,'
                    'Description="Estimated Allele Frequency, '
                    'for each ALT allele, in the same order as '
                    'listed">',
                    '##INFO=<ID=DP_HQ,Number=.,Type=Integer,'
                    'Description="HQ Read depth; sum of allelic '
                    'depths">']
vcf_headers = []
for line in inp:
    line = line.strip()
    if line.startswith("#"):
        vcf_headers.append(line)
        # adding new info headers just before the #CHROM line
        if line.startswith("#CHROM"):
            # remove any existing INFO headers for DP_HQ and ALT_AF
            vcf_headers = [e for e in vcf_headers if
                           "##INFO=<ID=DP_HQ" not in e]
            vcf_headers = [e for e in vcf_headers if
                           "##INFO=<ID=ALT_AF" not in e]

            # add the new DP_HQ and ALT_AF INFO headers just before
            # the #CHROM line
            for e in NEW_INFO_HEADERS:
                vcf_headers.insert(-1, e)
            # write out all the headers
            for e in vcf_headers:
                print(e, file=out)
        continue
    elems = line.split("\t")
    info = elems[INFO_CELL_INDEX]
    formatTokens = elems[FORMAT_CELL_INDEX].split(":")

    AO_Index = formatTokens.index('AO')
    RO_Index = formatTokens.index('RO')

    aos = elems[SAMPLE_DATA_INDEX].split(":")[AO_Index].\
        split(',')

    ro = elems[SAMPLE_DATA_INDEX].split(":")[RO_Index]

    for a, ao in enumerate(aos):
        aos[a] = int(ao)
    totalAlleleDepth = sum(aos) + int(ro)

    alternativeAlleles = elems[ALT_ALLELE_INDEX].split(",")
    alternateAlleleFrequencies = []
    for a, ao in enumerate(aos):
        # if the Allele depths are all zero, just call it zero,
        # we can't divide by 0
        if totalAlleleDepth == 0:
            alternateAlleleFrequencies.append("0")
        else:
            alternateAlleleFrequencies.append(str(
                round(ao / totalAlleleDepth, 3)))

    # update the info cell to have DP_HQ or ALT_AF entries, taking
    # care to remove any existing ones
    info_cell_items = elems[INFO_CELL_INDEX].split(";")
    # remove any existing DP_HQ or ALT_AF entries
    info_cell_items = [e for e in info_cell_items if "DP_HQ=" not in e]
    info_cell_items = [e for e in info_cell_items if "ALT_AF=" not in e]

    # add new entry for DP_HQ 
    info_cell_items.append("DP_HQ=" + str(totalAlleleDepth))
    
    # add new entries for ALT_AF separately for each alternative allele. 
    # This splitting ensures that the each alternative allele and its allele frequency are on a new row, and will thus get annotated (downstream) independently.
    for a, alt in enumerate(alternativeAlleles):
        tmp_elems = list(elems)
        tmp_elems[ALT_ALLELE_INDEX] = alt
        tmp_info_cell_items = list(info_cell_items)
        tmp_info_cell_items.append("ALT_AF=" + alternateAlleleFrequencies[a])
        tmp_elems[INFO_CELL_INDEX] = ";".join(tmp_info_cell_items)
        tmp_format_cell_items = list(tmp_elems[SAMPLE_DATA_INDEX].split(":"))
        tmp_elems[SAMPLE_DATA_INDEX] = ":".join(tmp_format_cell_items)
        print("\t".join(tmp_elems), file=out)
