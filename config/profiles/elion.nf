// profile credit to https://github.com/TheJacksonLaboratory/nf-tenx

singularity {
   enabled = true
   autoMounts = true
   cacheDir = '/projects/omics_share/meta/containers'
 }

process {
    executor = 'slurm'
    queue = 'batch'
    clusterOptions = '-q normal'
    module = 'slurm'
}

executor {
    $slurm {
        queueSize = 250
        // The number of tasks the executor will handle in a parallel manner
        submitRateLimit = '1 / 2 s'
        // Determines the max rate of job submission per time unit, for example '10sec' eg. max 10 jobs per second or '1/2 s' i.e. 1 job submissions every 2 seconds.
    }
}