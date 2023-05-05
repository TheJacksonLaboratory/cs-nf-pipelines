process ANNOTATE_BOOLEAN_PEAKS {
    tag "${antibody}"
    
    cpus 1
    memory 5.GB
    time '04:00:00'

    container 'ubuntu:20.04'

    input:
    tuple val(antibody), path(boolean_txt), path(homer_peaks)

    output:
    path '*.boolean.annotatePeaks.txt', emit: annotate_peaks_txt

    script:
    prefix="\$(echo ${boolean_txt} | sed 's/.boolean.txt//g')" 
    """
    cut -f2- ${homer_peaks} | awk 'NR==1; NR > 1 {print \$0 | "sort -T '.' -k1,1 -k2,2n"}' | cut -f6- > tmp.txt
    paste ${boolean_txt} tmp.txt > ${prefix}.boolean.annotatePeaks.txt
    """
}
