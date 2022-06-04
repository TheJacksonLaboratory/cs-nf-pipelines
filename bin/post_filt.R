#!/usr/bin/env Rscript

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

ins_s <- read_delim(paste(sample_name, "ins.s.bed", sep="."), delim = "\t",
                    col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

ins_e <- read_delim(paste(sample_name, "ins.e.bed", sep="."), delim = "\t",
                    col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

del_s <- read_delim(paste(sample_name, "del.s.bed", sep="."), delim = "\t",
                    col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

del_e <- read_delim(paste(sample_name, "del.e.bed", sep="."), delim = "\t",
                    col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

tra_e <- read_delim(paste(sample_name, "tra.e.bed", sep="."), delim = "\t",
                    col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

inv_e <- read_delim(paste(sample_name, "inv.e.bed", sep="."), delim = "\t",
                    col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

dup_e <- read_delim(paste(sample_name, "dup.e.bed", sep="."), delim = "\t",
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

ins_full <- ins %>%
    mutate(., sv_type = "insertion") %>%
    {if(nrow(ins_s) != 0) 
      mutate(., sanger = if_else(sv_name %in% ins_s$sv_name,  "N", "Y"))
     else(mutate(., sanger = "N")) } %>%
    {if(nrow(ins_e) != 0)
      mutate(., ensembl = if_else(sv_name %in% ins_e$sv_name, 
                                  "N", "Y"))
      else(mutate(., ensembl = "N")) } %>%
    {if(nrow(ins_genes) != 0)
      left_join(., select(ins_genes, sv_name, gene))
     else(.) } %>%
    {if(nrow(ins_exons) != 0)
      left_join(., select(ins_exons, sv_name, gene, exon_num) %>% 
              distinct() %>%
                arrange(., sv_name, gene, exon_num) %>%
                group_by(., sv_name, gene) %>% 
                summarize(., exons = paste(exon_num, collapse = ";")))
    else(.) }

del_full <- del %>%
    mutate(., sv_type = "deletion") %>%
    {if(nrow(del_s) != 0) 
      mutate(., sanger = if_else(sv_name %in% del_s$sv_name,  "N", "Y"))
     else(mutate(., sanger = "N")) } %>%
    {if(nrow(del_e) != 0)
      mutate(., ensembl = if_else(sv_name %in% del_e$sv_name, 
                                  "N", "Y"))
      else(mutate(., ensembl = "N")) } %>%
    {if(nrow(del_genes) != 0)
      left_join(., select(del_genes, sv_name, gene))
     else(.) } %>%
    {if(nrow(del_exons) != 0)
      left_join(., select(del_exons, sv_name, gene, exon_num) %>% 
              distinct() %>%
                arrange(., sv_name, gene, exon_num) %>%
                group_by(., sv_name, gene) %>% 
                summarize(., exons = paste(exon_num, collapse = ";")))
    else(.) }              

tra_full <- tra %>%
    mutate(., sv_type = "translocation") %>%
    mutate(., sanger = NA) %>%
    {if(nrow(tra_e) != 0)
      mutate(., ensembl = if_else(sv_name %in% tra_e$sv_name, 
                                  "N", "Y"))
      else(mutate(., ensembl = "N")) } %>%
    {if(nrow(tra_genes) != 0)
      left_join(., select(tra_genes, sv_name, gene))
     else(.) } %>%
    {if(nrow(tra_exons) != 0)
      left_join(., select(tra_exons, sv_name, gene, exon_num) %>% 
              distinct() %>%
                arrange(., sv_name, gene, exon_num) %>%
                group_by(., sv_name, gene) %>% 
                summarize(., exons = paste(exon_num, collapse = ";")))
     else(.)}

inv_full <- inv %>%
    mutate(., sv_type = "inversion") %>%
    mutate(., sanger = NA) %>%
    {if(nrow(inv_e) != 0)
      mutate(., ensembl = if_else(sv_name %in% inv_e$sv_name, 
                                  "N", "Y"))
      else(mutate(., ensembl = "N")) } %>%
    {if(nrow(inv_genes) != 0)
      left_join(., select(inv_genes, sv_name, gene))
     else(.) } %>%
    {if(nrow(inv_exons) != 0)
      left_join(., select(inv_exons, sv_name, gene, exon_num) %>% 
              distinct() %>%
                arrange(., sv_name, gene, exon_num) %>%
                group_by(., sv_name, gene) %>% 
                summarize(., exons = paste(exon_num, collapse = ";")))
     else(.)}

dup_full <- dup %>%
    mutate(., sv_type = "duplication") %>%
    mutate(., sanger = NA) %>%
    {if(nrow(dup_e) != 0)
      mutate(., ensembl = if_else(sv_name %in% dup_e$sv_name, 
                                  "N", "Y"))
      else(mutate(., ensembl = "N")) } %>%
    {if(nrow(dup_genes) != 0)
      left_join(., select(dup_genes, sv_name, gene))
     else(.) } %>%
    {if(nrow(dup_exons) != 0)
      left_join(., select(dup_exons, sv_name, gene, exon_num) %>% 
              distinct() %>%
                arrange(., sv_name, gene, exon_num) %>%
                group_by(., sv_name, gene) %>% 
                summarize(., exons = paste(exon_num, collapse = ";")))
     else(.)}

bind_rows(ins_full, del_full, tra_full, inv_full, dup_full) %>%
    mutate(., pos = start + 1) %>%
    mutate(., sv_size = abs(stop - start) + 1) %>%
    relocate(chr, pos, sv_name, sv_size, sv_type, 
             sanger, ensembl, gene, exons) %>%
    select(., -start, -stop) %>%
    left_join(., read_delim(sv_olap_fh,
                        delim = "\t",
                        col_types = cols("chr" = "c")) %>%
             select(., -chr, -pos), 
             by = c("sv_name" = "SV")) %>%
    arrange(., chr, pos) %>%
    write_csv(., surv_out_fh)