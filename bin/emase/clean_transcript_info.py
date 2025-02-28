import sys
import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser(prog = 'clean_transcript_info.py',
                                description = 'This script takes the "pooled" transcript list, parses by haplotype, and ensures that all haplotypes have all transcripts.\nFurther, it then takes the complete transcript list and parses it to remove haplotype IDs. This second file is for use in bam2emase')

parser.add_argument('--input-transcript-list', action='store', dest='input_transcript_list', type=str, required=True, help='Input transcript list from prepare-emase run with multiple genomes')
parser.add_argument('--haplotype-transcript-output', action='store', dest='haplotype_outfile', type=str, required=True, default=None, help='Output file to replace "emase.pooled.transcripts.info"')
parser.add_argument('--full-transcript-output', action='store', dest='transcript_outfile', type=str, required=True, default=None, help='Output file to replace "emase.transcripts.info"')
parser.add_argument('--haplotype-list', action='store', dest='haplotypes', type=str, required=True, default='', help='Haplotype list (e.g., A,B,C,D,E,F,G,H)')
# parser.add_argument('--exclude_partial-trascripts', dest='exclude_partials', default=False, action='store_true', help='Remove transcripts not shared among all haplotypes')

args = parser.parse_args()

# if (exclude_partials): 
#     pass
#     # If we don't want partial transcripts (i.e., transcripts not found in all), should we remove them? 

haplotypes = args.haplotypes.split(",")
    # split haplotype lists

transcript_table = pd.read_csv(args.input_transcript_list, sep='\t', header = None, names = ['transcript_id', 'length'])
    # import transcript table

transcript_table_working = transcript_table.copy()
    # avoid pandas slice error by making a copy, also original is used in final cat statement. 

transcript_table_working[['Transcript', 'Haplotype']] = transcript_table_working['transcript_id'].str.rsplit(pat='_', n=1, expand=True)
    # split transcript and haplotype ID into new columns. 

pivot_transcripts = transcript_table_working.pivot(index='Transcript', columns='Haplotype', values='length')
    # pivot table generated on Transcript
    # Transcript x Haplotype. NAs in table are missing transcripts within that haplotype. 

missing_transcripts = pivot_transcripts[pivot_transcripts.isna().any(axis=1)]
    # subset pivot to only columns with NA data (i.e., missing transcripts within that haplotype). 

if (not missing_transcripts.empty):
    # if the transcript list is incomplete, add back those that are missing. 

    stacked = pd.melt(missing_transcripts.reset_index(), id_vars='Transcript',value_vars=haplotypes)
        # convert wide to long on transcript. 

    transcripts_to_add = stacked[stacked['value'].isna()].copy()
    transcripts_to_add['transcript_id'] = transcripts_to_add[['Transcript', 'Haplotype']].agg('_'.join, axis=1)
    transcripts_to_add['value'] = 0.0
    transcripts_to_add['length'] = transcripts_to_add['value']
    transcripts_to_merge = transcripts_to_add[["transcript_id", "length"]]
        # Filter the long list to only NA. Reconstruct a 'transcript_haplotype' string. Add value = 0.0. change 'value' to 'length'. Keep just 2 columns. 
        # long list with NA are transcripts missing from the original 'emase.pooled.transcripts.info'

    updated_transcript_table = pd.concat([transcript_table, transcripts_to_merge])
    updated_transcript_table.to_csv(args.haplotype_outfile, sep='\t', index = False, header = False)
        # concat the missing transcripts to to the original table, and output to file. 
else:
    # else, if all transcripts are present, just write the file out as it is. 
    transcript_table.to_csv(args.haplotype_outfile, sep='\t', index = False, header = False)

full_transcript_list = pd.DataFrame(pivot_transcripts.sort_index().index)
full_transcript_list['length'] = 0.0
full_transcript_list.to_csv(args.transcript_outfile, sep='\t', index = False, header = False)
    # the 'index' list from the pivot table is a list of all transcripts across all haplotypes. 
    # this list can be used to generate the file to replace "emase.transcripts.info"
    # transcript length is meaningless here. It is not used in bam2emase or bam-to-emase, only the first column 'ID' is used. 