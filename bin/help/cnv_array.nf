def help() {
    println '''
Parameter | Default | Description

--bpm_file | /<PATH> | The path to the BPM file.
--egt_file | /<PATH> | The path to the EGT file.
-w | /<PATH> | The directory for intermediary files and Nextflow processes. This directory can become quite large. Ensure ample storage.
--help | false | Print this help message and exit.

+gtc2vcf --no-version -Ou \
--bpm | /<PATH> | The path to the BPM file 
--csv | /<PATH> | The path to csv file 
--egt | /<PATH> | The patht to egt file 
--gtcs | /<PATH> | The path to gtgc output 
--fasta-ref | /<PATH> | The path to reference 
--extra | /<PATH> | The path to output directory    

bcftools sort -Ou -T ./bcftools. | \
bcftools norm --no-version -Ob -c x -f ${fasta} | \
tee bcftools_convert.bcf | \
    
bcftools index --force --output bcftools_convert.bcf.csi
bcftools convert -O v -o bcftools_convert.vcf bcftools_convert.bcf
    """

'''
}

