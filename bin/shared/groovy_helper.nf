// Helper Function: Return A Match Of a Pattern
def find(files, match) {
    for (i in files) {
       if (i ==~ match) {
           println(i)
           return(i)
       }
    }
}

// my_genome_bam = find(files, /.*genome.bam/)
