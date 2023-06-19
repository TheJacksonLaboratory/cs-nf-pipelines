#!/usr/bin/env Rscript

#zout <- file("messages.Rout", open = "wt")
#sink(zout, type = "message")

library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

if(length(args) == 1) {
  sample_name <- args[1]
} else {
  stop("Incorrect number of arguments. Please supply only the sample name as an argument.")
}

sv_summary_fh <- paste(sample_name, "survivor_summary.csv", sep = ".")
sv_olap_fh <- paste(sample_name, "merged.overlap.annotated.txt", sep = ".")
surv_out_fh <- paste(sample_name, "survivor_joined_results.csv", sep = ".")


sv_summary <- read_csv(sv_summary_fh,
                       col_types = cols("chr" = "c"))


ins <- read_delim(paste(sample_name, "ins.bed", sep="."), delim = "\t",
                  col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

del <- read_delim(paste(sample_name, "del.bed", sep="."), delim = "\t",
                  col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

tra <- read_delim(paste(sample_name, "tra.bed", sep="."), delim = "\t",
                  col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

dup <- read_delim(paste(sample_name, "dup.bed", sep="."), delim = "\t",
                  col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

inv <- read_delim(paste(sample_name, "inv.bed", sep="."), delim = "\t",
                  col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

ins_c <- read_delim(paste(sample_name, "ins.c.bed", sep="."), delim = "\t",
                    col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

del_c <- read_delim(paste(sample_name, "del.c.bed", sep="."), delim = "\t",
                    col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

inv_c <- read_delim(paste(sample_name, "inv.c.bed", sep="."), delim = "\t",
                    col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

ins_genes <- read_delim(paste(sample_name, "ins.genes.bed", sep="."), delim = "\t",
                        col_names = c("chr", "start", "stop", "sv_name",
                                      "ref_chr", "ref_start",
                                      "ref_stop", "gene"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c",
                                   "ref_chr" = "c",
                                   "ref_start" = "d",
                                   "ref_stop" = "d",
                                   "gene" = "c"))

del_genes <- read_delim(paste(sample_name, "del.genes.bed", sep="."), delim = "\t",
                        col_names = c("chr", "start", "stop", "sv_name",
                                      "ref_chr", "ref_start",
                                      "ref_stop", "gene"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c",
                                   "ref_chr" = "c",
                                   "ref_start" = "d",
                                   "ref_stop" = "d",
                                   "gene" = "c"))

tra_genes <- read_delim(paste(sample_name, "tra.genes.bed", sep="."), delim = "\t",
                        col_names = c("chr", "start", "stop", "sv_name",
                                      "ref_chr", "ref_start",
                                      "ref_stop", "gene"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c",
                                   "ref_chr" = "c",
                                   "ref_start" = "d",
                                   "ref_stop" = "d",
                                   "gene" = "c"))

dup_genes <- read_delim(paste(sample_name, "dup.genes.bed", sep="."), delim = "\t",
                        col_names = c("chr", "start", "stop", "sv_name",
                                      "ref_chr", "ref_start",
                                      "ref_stop", "gene"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c",
                                   "ref_chr" = "c",
                                   "ref_start" = "d",
                                   "ref_stop" = "d",
                                   "gene" = "c"))

inv_genes <- read_delim(paste(sample_name, "inv.genes.bed", sep="."), delim = "\t",
                        col_names = c("chr", "start", "stop", "sv_name",
                                      "ref_chr", "ref_start",
                                      "ref_stop", "gene"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c",
                                   "ref_chr" = "c",
                                   "ref_start" = "d",
                                   "ref_stop" = "d",
                                   "gene" = "c"))

ins_exons <- read_delim(paste(sample_name, "ins.exons.bed", sep="."), delim = "\t",
                        col_names = c("chr", "start", "stop", "sv_name",
                                      "ref_chr", "ref_start",
                                      "ref_stop", "direction",
                                      "gene_id", "transcript_id",
                                      "exon_num", "gene"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c",
                                   "ref_chr" = "c",
                                   "ref_start" = "d",
                                   "ref_stop" = "d",
                                   "direction" = "c",
                                   "gene_id" = "c",
                                   "transcript_id" = "c",
                                   "exon_num" = "d",
                                   "gene" = "c"))

del_exons <- read_delim(paste(sample_name, "del.exons.bed", sep="."), delim = "\t",
                        col_names = c("chr", "start", "stop", "sv_name",
                                      "ref_chr", "ref_start",
                                      "ref_stop", "direction",
                                      "gene_id", "transcript_id",
                                      "exon_num", "gene"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c",
                                   "ref_chr" = "c",
                                   "ref_start" = "d",
                                   "ref_stop" = "d",
                                   "direction" = "c",
                                   "gene_id" = "c",
                                   "transcript_id" = "c",
                                   "exon_num" = "d",
                                   "gene" = "c"))

tra_exons <- read_delim(paste(sample_name, "tra.exons.bed", sep="."), delim = "\t",
                        col_names = c("chr", "start", "stop", "sv_name",
                                      "ref_chr", "ref_start",
                                      "ref_stop", "direction",
                                      "gene_id", "transcript_id",
                                      "exon_num", "gene"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c",
                                   "ref_chr" = "c",
                                   "ref_start" = "d",
                                   "ref_stop" = "d",
                                   "direction" = "c",
                                   "gene_id" = "c",
                                   "transcript_id" = "c",
                                   "exon_num" = "d",
                                   "gene" = "c"))

dup_exons <- read_delim(paste(sample_name, "dup.exons.bed", sep="."), delim = "\t",
                        col_names = c("chr", "start", "stop", "sv_name",
                                      "ref_chr", "ref_start",
                                      "ref_stop", "direction",
                                      "gene_id", "transcript_id",
                                      "exon_num", "gene"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c",
                                   "ref_chr" = "c",
                                   "ref_start" = "d",
                                   "ref_stop" = "d",
                                   "direction" = "c",
                                   "gene_id" = "c",
                                   "transcript_id" = "c",
                                   "exon_num" = "d",
                                   "gene" = "c"))

inv_exons <- read_delim(paste(sample_name, "inv.exons.bed", sep="."), delim = "\t",
                        col_names = c("chr", "start", "stop", "sv_name",
                                      "ref_chr", "ref_start",
                                      "ref_stop", "direction",
                                      "gene_id", "transcript_id",
                                      "exon_num", "gene"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c",
                                   "ref_chr" = "c",
                                   "ref_start" = "d",
                                   "ref_stop" = "d",
                                   "direction" = "c",
                                   "gene_id" = "c",
                                   "transcript_id" = "c",
                                   "exon_num" = "d",
                                   "gene" = "c"))

ins_regulatory <- read_delim(paste(sample_name, "ins.regulatory.bed", sep="."), delim = "\t",
                        col_names = c("chr", "start", "stop", "sv_name",
                                      "ref_chr", "ref_build",
                                      "ref_type", "ref_start",
                                      "ref_stop", "blank_1",
                                      "blank_2", "blank_3",
                                      "meta_string"),
                  col_types = cols("chr" = "c",
                                   "start" = "c",
                                   "stop" = "c",
                                   "sv_name" = "c",
                                   "ref_chr" = "c",
                                   "ref_build" = "c",
                                   "ref_type" = "c",
                                   "ref_start" = "c",
                                   "ref_stop" = "c",
                                   "blank_1" = "c",
                                   "blank_2" = "c",
                                   "blank_3" = "c",
                                   "meta_string"  = "c")) %>%
                  {if(nrow(.) != 0)
                  separate(., meta_string, 
                  into = c("ID", "bound_end", "bound_start", "description", "feature_type"), 
                            sep = ";") %>% 
                  mutate(., ID = str_replace(ID, "ID=regulatory_region:", ""),
                            description = str_replace(description, "description=", ""),
                            reg_out = str_c(ID, description, sep = ":")) 
                   else(.)}

del_regulatory <- read_delim(paste(sample_name, "del.regulatory.bed", sep="."), delim = "\t",
                        col_names = c("chr", "start", "stop", "sv_name",
                                      "ref_chr", "ref_build",
                                      "ref_type", "ref_start",
                                      "ref_stop", "blank_1",
                                      "blank_2", "blank_3",
                                      "meta_string"),
                  col_types = cols("chr" = "c",
                                   "start" = "c",
                                   "stop" = "c",
                                   "sv_name" = "c",
                                   "ref_chr" = "c",
                                   "ref_build" = "c",
                                   "ref_type" = "c",
                                   "ref_start" = "c",
                                   "ref_stop" = "c",
                                   "blank_1" = "c",
                                   "blank_2" = "c",
                                   "blank_3" = "c",
                                   "meta_string"  = "c")) %>%
                  {if(nrow(.) != 0)
                  separate(., meta_string, 
                  into = c("ID", "bound_end", "bound_start", "description", "feature_type"), 
                            sep = ";") %>% 
                  mutate(., ID = str_replace(ID, "ID=regulatory_region:", ""),
                            description = str_replace(description, "description=", ""),
                            reg_out = str_c(ID, description, sep = ":")) 
                   else(.)}

inv_regulatory <- read_delim(paste(sample_name, "inv.regulatory.bed", sep="."), delim = "\t",
                        col_names = c("chr", "start", "stop", "sv_name",
                                      "ref_chr", "ref_build",
                                      "ref_type", "ref_start",
                                      "ref_stop", "blank_1",
                                      "blank_2", "blank_3",
                                      "meta_string"),
                  col_types = cols("chr" = "c",
                                   "start" = "c",
                                   "stop" = "c",
                                   "sv_name" = "c",
                                   "ref_chr" = "c",
                                   "ref_build" = "c",
                                   "ref_type" = "c",
                                   "ref_start" = "c",
                                   "ref_stop" = "c",
                                   "blank_1" = "c",
                                   "blank_2" = "c",
                                   "blank_3" = "c",
                                   "meta_string"  = "c")) %>%
                  {if(nrow(.) != 0)
                  separate(., meta_string, 
                  into = c("ID", "bound_end", "bound_start", "description", "feature_type"), 
                            sep = ";") %>% 
                  mutate(., ID = str_replace(ID, "ID=regulatory_region:", ""),
                            description = str_replace(description, "description=", ""),
                            reg_out = str_c(ID, description, sep = ":")) 
                   else(.)}             

dup_regulatory <- read_delim(paste(sample_name, "dup.regulatory.bed", sep="."), delim = "\t",
                        col_names = c("chr", "start", "stop", "sv_name",
                                      "ref_chr", "ref_build",
                                      "ref_type", "ref_start",
                                      "ref_stop", "blank_1",
                                      "blank_2", "blank_3",
                                      "meta_string"),
                  col_types = cols("chr" = "c",
                                   "start" = "c",
                                   "stop" = "c",
                                   "sv_name" = "c",
                                   "ref_chr" = "c",
                                   "ref_build" = "c",
                                   "ref_type" = "c",
                                   "ref_start" = "c",
                                   "ref_stop" = "c",
                                   "blank_1" = "c",
                                   "blank_2" = "c",
                                   "blank_3" = "c",
                                   "meta_string"  = "c")) %>%
                  {if(nrow(.) != 0)
                  separate(., meta_string, 
                  into = c("ID", "bound_end", "bound_start", "description", "feature_type"), 
                            sep = ";") %>% 
                  mutate(., ID = str_replace(ID, "ID=regulatory_region:", ""),
                            description = str_replace(description, "description=", ""),
                            reg_out = str_c(ID, description, sep = ":")) 
                   else(.)}

tra_regulatory <- read_delim(paste(sample_name, "tra.regulatory.bed", sep="."), delim = "\t",
                        col_names = c("chr", "start", "stop", "sv_name",
                                      "ref_chr", "ref_build",
                                      "ref_type", "ref_start",
                                      "ref_stop", "blank_1",
                                      "blank_2", "blank_3",
                                      "meta_string"),
                  col_types = cols("chr" = "c",
                                   "start" = "c",
                                   "stop" = "c",
                                   "sv_name" = "c",
                                   "ref_chr" = "c",
                                   "ref_build" = "c",
                                   "ref_type" = "c",
                                   "ref_start" = "c",
                                   "ref_stop" = "c",
                                   "blank_1" = "c",
                                   "blank_2" = "c",
                                   "blank_3" = "c",
                                   "meta_string"  = "c")) %>%
                  {if(nrow(.) != 0)
                  separate(., meta_string, 
                  into = c("ID", "bound_end", "bound_start", "description", "feature_type"), 
                            sep = ";") %>% 
                  mutate(., ID = str_replace(ID, "ID=regulatory_region:", ""),
                            description = str_replace(description, "description=", ""),
                            reg_out = str_c(ID, description, sep = ":")) 
                   else(.)}
                                                       
ins_full <- ins %>%
    mutate(., sv_type = "insertion") %>%
    {if(nrow(ins_c) != 0) 
      mutate(., beck = if_else(sv_name %in% ins_c$sv_name,  "Y", "N"))
     else(mutate(., beck = "N")) } %>%  
    {if(nrow(ins_genes) != 0)
      left_join(., select(ins_genes, sv_name, gene))
     else(.) } %>%
    {if(nrow(ins_exons) != 0)
      left_join(., select(ins_exons, sv_name, gene, exon_num) %>% 
              distinct() %>%
                arrange(., sv_name, gene, exon_num) %>%
                group_by(., sv_name, gene) %>% 
                summarize(., exons = paste(exon_num, collapse = ";")))
     else(.) } %>%
    {if(nrow(ins_regulatory) != 0)
       mutate(., regulatory_regions = if_else(sv_name %in% ins_regulatory$sv_name,  "Y", "N"))
        else(mutate(., regulatory_regions = "N")) }

del_full <- del %>%
    mutate(., sv_type = "deletion") %>%
    {if(nrow(del_c) != 0) 
      mutate(., beck = if_else(sv_name %in% del_c$sv_name,  "Y", "N"))
     else(mutate(., beck = "N")) } %>%
    {if(nrow(del_genes) != 0)
      left_join(., select(del_genes, sv_name, gene))
     else(.) } %>%
    {if(nrow(del_exons) != 0)
      left_join(., select(del_exons, sv_name, gene, exon_num) %>% 
              distinct() %>%
                arrange(., sv_name, gene, exon_num) %>%
                group_by(., sv_name, gene) %>% 
                summarize(., exons = paste(exon_num, collapse = ";")))
    else(.) } %>%
    {if(nrow(del_regulatory) != 0)
       mutate(., regulatory_regions = if_else(sv_name %in% del_regulatory$sv_name,  "Y", "N"))
        else(mutate(., regulatory_regions = "N")) }                

tra_full <- tra %>%
    mutate(., sv_type = "translocation") %>%
    mutate(., beck = NA) %>%
    {if(nrow(tra_genes) != 0)
      left_join(., select(tra_genes, sv_name, gene))
     else(.) } %>%
    {if(nrow(tra_exons) != 0)
      left_join(., select(tra_exons, sv_name, gene, exon_num) %>% 
              distinct() %>%
                arrange(., sv_name, gene, exon_num) %>%
                group_by(., sv_name, gene) %>% 
                summarize(., exons = paste(exon_num, collapse = ";")))
     else(.)} %>%
    {if(nrow(tra_regulatory) != 0)
       mutate(., regulatory_regions = if_else(sv_name %in% tra_regulatory$sv_name,  "Y", "N"))
        else(mutate(., regulatory_regions = "N")) }

inv_full <- inv %>%
    mutate(., sv_type = "inversion") %>%
    {if(nrow(inv_c) != 0)
      mutate(., beck = if_else(sv_name %in% inv_c$sv_name, 
                                  "Y", "N"))
      else(mutate(., beck = "N")) } %>%
    {if(nrow(inv_genes) != 0)
      left_join(., select(inv_genes, sv_name, gene))
     else(.) } %>%
    {if(nrow(inv_exons) != 0)
      left_join(., select(inv_exons, sv_name, gene, exon_num) %>% 
              distinct() %>%
                arrange(., sv_name, gene, exon_num) %>%
                group_by(., sv_name, gene) %>% 
                summarize(., exons = paste(exon_num, collapse = ";")))
     else(.)} %>%
    {if(nrow(inv_regulatory) != 0)
       mutate(., regulatory_regions = if_else(sv_name %in% inv_regulatory$sv_name,  "Y", "N"))
        else(mutate(., regulatory_regions = "N")) }

dup_full <- dup %>%
    mutate(., sv_type = "duplication") %>%
    mutate(., beck = NA) %>%
    {if(nrow(dup_genes) != 0)
      left_join(., select(dup_genes, sv_name, gene))
     else(.) } %>%
    {if(nrow(dup_exons) != 0)
      left_join(., select(dup_exons, sv_name, gene, exon_num) %>% 
              distinct() %>%
                arrange(., sv_name, gene, exon_num) %>%
                group_by(., sv_name, gene) %>% 
                summarize(., exons = paste(exon_num, collapse = ";")))
     else(.)} %>%
    {if(nrow(dup_regulatory) != 0)
       mutate(., regulatory_regions = if_else(sv_name %in% dup_regulatory$sv_name,  "Y", "N"))
        else(mutate(., regulatory_regions = "N")) }

bind_rows(ins_full, del_full, tra_full, inv_full, dup_full) %>%
    mutate(., pos = start + 1) %>%
    mutate(., sv_size = abs(stop - start) + 1) %>%
    relocate(chr, pos, sv_name, sv_size, sv_type, 
             beck, regulatory_regions, gene, exons) %>%
    select(., -start, -stop) %>%
    arrange(., chr, pos) %>%
    mutate(., gene = case_when(sv_size > 5000000 ~ "sv_too_large",
                               TRUE ~ gene),
              exons = case_when(sv_size > 5000000 ~ "sv_too_large",
                               TRUE ~ exons),
              regulatory_regions = case_when(sv_size > 5000000 ~ "sv_too_large",
                               TRUE ~ regulatory_regions)) %>%
    distinct() %>%
    write_csv(., surv_out_fh)