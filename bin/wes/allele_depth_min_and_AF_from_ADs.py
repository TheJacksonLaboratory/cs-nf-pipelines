#!/usr/bin/env python3
"""

GOALS of script:
    recompute the locus depth from the allele-depths, and filter
    based on a minimumTotalAlleleDepth

    Add Estimated Allele Frequency (ALT_AF) to the info cell

The FORMAT column indicates the order of the fields in the following,
sample, column.
We need to find the AD (allele depth) field in the FORMAT, and
then process the corresponding field in the sample column.  It will
be two comma separated integers.
Add them together and compare to the required minimum depth.
If greater or equal, output the line as-is.
If less, change the FILTER column value to MinDP and output the line.

"""
import sys

ALT_ALLELE_INDEX = 4
FILTER_CELL_INDEX = 6
INFO_CELL_INDEX = 7
FORMAT_CELL_INDEX = 8
SAMPLE_DATA_INDEX = int(sys.argv[4])

# expecting sys.argv to look something like
if len(sys.argv) != 5:
    print("Incorrect number of args!")
    print("expected usage:")
    print("allele_depth_min_and_AF_from_ADs.py inputFile outputFile minDP sample_column_number")
    raise Exception("Incorrect number of args!")

inp = open(sys.argv[1], 'r')
out = open(sys.argv[2], 'w')
minimumTotalAlleleDepth = int(sys.argv[3])

NEW_INFO_HEADERS = ['##INFO=<ID=ALT_AF,Number=A,Type=Float,'
                    'Description="Estimated Allele Frequency, '
                    'for each ALT allele, in the same order as '
                    'listed">',
                    '##INFO=<ID=DP_HQ,Number=1,Type=Integer,'
                    'Description="HQ Read depth; sum of allelic '
                    'depths">',
                    '##FILTER=<ID=minDP,Description="DP_HQ < %s">' %minimumTotalAlleleDepth]
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
            vcf_headers = [e for e in vcf_headers if
                           "##FILTER=<ID=minDP" not in e]
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
    alleleDepthIndex = formatTokens.index('AD')
    ads = elems[SAMPLE_DATA_INDEX].split(":")[alleleDepthIndex].\
        split(',')
    for a, ad in enumerate(ads):
        ads[a] = int(ad)
    totalAlleleDepth = sum(ads)
    # update filter cell with minDP to fail if the totalAlleleDepth
    # is less than the specified value
    if minimumTotalAlleleDepth is not None and \
            totalAlleleDepth < minimumTotalAlleleDepth:
        # ignore case when looking for an existing minDP in the
        # filter cell
        if "minDP".lower() in (elems[FILTER_CELL_INDEX]).lower():
            # already contains minDP so don't add another
            pass
        elif elems[FILTER_CELL_INDEX] == ".":
            elems[FILTER_CELL_INDEX] = "minDP"
        elif elems[FILTER_CELL_INDEX] == "PASS":
            elems[FILTER_CELL_INDEX] = "minDP"
        else:
            elems[FILTER_CELL_INDEX] += ";minDP"
    alternativeAlleles = elems[ALT_ALLELE_INDEX].split(",")
    alternateAlleleFrequencies = []
    for a, ad in enumerate(ads):
        # if the Allele depths are all zero, just call it zero,
        # we can't divide by 0
        if totalAlleleDepth == 0:
            alternateAlleleFrequencies.append("0")
        else:
            alternateAlleleFrequencies.append(str(
                round(ad * (100.0 / totalAlleleDepth))))
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
        tmp_info_cell_items.append("ALT_AF=" + alternateAlleleFrequencies[a+1])
        tmp_elems[INFO_CELL_INDEX] = ";".join(tmp_info_cell_items)
        tmp_format_cell_items = list(tmp_elems[SAMPLE_DATA_INDEX].split(":"))
        tmp_format_cell_items[alleleDepthIndex] = "%s,%s" %(ads[0],ads[a+1])
        tmp_elems[SAMPLE_DATA_INDEX] = ":".join(tmp_format_cell_items)
        tmp_elems[SAMPLE_DATA_INDEX].split(":")[alleleDepthIndex] = "%s,%s" %(ads[0],ads[a+1])
        print("\t".join(tmp_elems), file=out)