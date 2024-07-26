def help() {
    println '''
Parameter        | Default | Description
-----------------|---------|---------------------------------------------------------------------------
--bpm_file       | /<PATH> | The path to the BPM file.
--egt_file       | /<PATH> | The path to the EGT file.
-w               | /<PATH> | The directory for intermediary files and Nextflow processes. Ensure ample storage.
--help           | false   | Print this help message and exit.

--bpm            | /<PATH> | The path to the BPM file.
--csv            | /<PATH> | The path to the CSV file.
--egt            | /<PATH> | The path to the EGT file.
--gtcs           | /<PATH> | The path to GTC output.
--fasta-ref      | /<PATH> | The path to the reference FASTA file.
--extra          | /<PATH> | The path to the output directory.

// Example usage of BCFTOOLS:
// ---------------------------
// bcftools sort -Ou -T ./bcftools. | \\
// bcftools norm --no-version -Ob -c x -f <fasta> | \\
// tee bcftools_convert.bcf | \\
// bcftools index --force --output bcftools_convert.bcf.csi
// bcftools convert -O v -o bcftools_convert.vcf bcftools_convert.bcf
'''
}
