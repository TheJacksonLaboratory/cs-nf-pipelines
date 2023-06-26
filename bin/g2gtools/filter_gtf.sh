#!/bin/sh

gtf_in=$1

# Define string patterns for GTF tags
# NOTES:
# - Since GENCODE release 31/M22 (Ensembl 97), the "lincRNA" and "antisense"
#   biotypes are part of a more generic "lncRNA" biotype.
# - These filters are relevant only to GTF files from GENCODE. The GTFs from
#   Ensembl release 98 have the following differences:
#   - The names "gene_biotype" and "transcript_biotype" are used instead of
#     "gene_type" and "transcript_type".
#   - Readthrough transcripts are present but are not marked with the
#     "readthrough_transcript" tag.
BIOTYPE_PATTERN="($2)"

FILE_NAME=`basename $1 .gtf`

GENE_PATTERN="gene_biotype \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_biotype \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""

# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_biotype (biotype)
#   - allowable transcript_biotype (biotype)
#   - no "readthrough_transcript" tag
# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.
cat "$gtf_in" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | sed -E 's/.*gene_id "([^"]+)".*/\1/' \
    | sort \
    | uniq \
    > gene_allowlist

# Filter the GTF file based on the gene allowlist
gtf_filtered=$FILE_NAME".filtered.gtf"

# Copy header lines beginning with "#"
grep -E "^#" "$gtf_in" > "$gtf_filtered"

# ripgrep is much faster
rg -Ff gene_allowlist "$gtf_in" >> "$gtf_filtered"
