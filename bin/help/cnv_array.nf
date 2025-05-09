def help() {
    println '''
Parameter | Default | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /fastscratch or other directory with ample storage.
--gtc_csv | '/projects/compsci/omics_share/human/GRCh38/supporting_files/cnv_array/HumanCytoSNP-12v2-1_L1.csv' | Genotype Call (GTC) manifest for IDAT conversion. Provided by Illumina.
--bpm_file | '/projects/compsci/omics_share/human/GRCh38/supporting_files/cnv_array/HumanCytoSNP-12v2-1_L1.bpm' | Manifest file describing the SNP or probe content on a BeadChip. Provided by Illumina.
--egt_file | '/projects/compsci/omics_share/human/GRCh38/supporting_files/cnv_array/HumanCytoSNP-12v2-1_L.egt' | Cluster file describing the cluster positions for the Illumina genotyping array. Provided by Illumina.
--ref_fa | '/projects/compsci/omics_share/human/GRCh38/genome/sequence/gatk/Homo_sapiens_assembly38.fasta' | The reference fasta file. Reference FASTA build should match Illumina provided files.
--snp_platform | 'IlluminaCytoSNP' | SNP platform supported by ASCAT. See full supported list here: https://github.com/VanLoo-lab/ascat?tab=readme-ov-file#supported-arrays-without-matched-germline
--gc_file | '/projects/compsci/omics_share/human/GRCh38/supporting_files/cnv_array/HumanCytoSNP-12v2-1_GCcontent_validSNPloci.txt' | ASCAT’s GC correction file, generated from scripts at https://github.com/VanLoo-lab/ascat/tree/master/LogRcorrection
--rt_file | '/projects/compsci/omics_share/human/GRCh38/supporting_files/cnv_array/HumanCytoSNP-12v2-1_ReplicationTiming_SNPloci_hg38.txt' | ASCAT’s replication timing file, generated from scripts at https://github.com/VanLoo-lab/ascat/tree/master/LogRcorrection
--chrArm | '/projects/compsci/omics_share/human/GRCh38/supporting_files/cnv_array/GRCh38_chromosome_arm.txt' | Chromosome arm locations, used in CNV segment annotation.
--cnvGeneFile | '/projects/compsci/omics_share/human/GRCh38/supporting_files/cnv_array/biomaRt_GRCh38_ensemblv102_CNVgeneAnnotations_primaryChroms.txt'
'''
}
