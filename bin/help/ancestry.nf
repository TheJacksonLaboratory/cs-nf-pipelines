def help(){
  println '''
Parameter | Default | Description

--pubdir | /<PATH> | The directory that the saved outputs will be stored.
--organize_by | sample | How to organize the output folder structure. Options: sample or analysis.
--cacheDir | /projects/omics_share/meta/containers | This is directory that contains cached Singularity containers. JAX users should not change this parameter.
-w | /<PATH> | The directory that all intermediary files and nextflow processes utilize. This directory can become quite large. This should be a location on /flashscratch or other directory with ample storage.

--csv_input | null | Provide a CSV manifest file with the header: "sampleID,lane,fastq_1,fastq_2". See the repository wiki for an example file. Fastq_2 is optional and used only in PE data. Fastq files can either be absolute paths to local files, or URLs to remote files. If remote URLs are provided, `--download_data` must be specified.
--download_data | null | Requires `--csv_input`. When specified, read data in the CSV manifest will be downloaded from provided URLs. 
--genotype_targets | '/projects/compsci/omics_share/human/GRCh38/supporting_files/ancestry_panel/snp_panel_v2_targets_annotations.snpwt.bed.gz' | Target SNP bed file for the ancestry panel. Can contain annotation information. 
--snpID_list | '/projects/compsci/omics_share/human/GRCh38/supporting_files/ancestry_panel/snp_panel_v2.list' | Target SNPs in list used in BCFtools filtering step
--snp_annotations | '/projects/compsci/omics_share/human/GRCh38/supporting_files/ancestry_panel/snp_panel_v2_targets_annotations.snpwt.bed.gz' | Target SNP bed file with annotations for the ancestry panel.
--snpweights_panel | '/projects/compsci/omics_share/human/GRCh38/supporting_files/ancestry_panel/ancestry_panel_v2.snpwt' | SNP weights panel in the appropriate format
'''
}
