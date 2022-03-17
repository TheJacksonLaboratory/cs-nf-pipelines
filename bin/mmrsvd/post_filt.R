#!/usr/bin/env Rscript

library(tidyverse)

sv_summary <- read_csv("summary.csv")

ins <- read_delim("ins.bed", delim = "\t",
                  col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

del <- read_delim("del.bed", delim = "\t",
                  col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

tra <- read_delim("tra.bed", delim = "\t",
                  col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

dup <- read_delim("dup.bed", delim = "\t",
                  col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

inv <- read_delim("inv.bed", delim = "\t",
                  col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

ins_s <- read_delim("ins.s.bed", delim = "\t",
                    col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

ins_e <- read_delim("ins.e.bed", delim = "\t",
                    col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

del_s <- read_delim("del.s.bed", delim = "\t",
                    col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

del_e <- read_delim("del.e.bed", delim = "\t",
                    col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

tra_e <- read_delim("tra.e.bed", delim = "\t",
                    col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

inv_e <- read_delim("inv.e.bed", delim = "\t",
                    col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))

dup_e <- read_delim("dup.e.bed", delim = "\t",
                    col_names = c("chr", "start", "stop", "sv_name"),
                  col_types = cols("chr" = "c",
                                   "start" = "d",
                                   "stop" = "d",
                                   "sv_name" = "c"))


ins_genes <- read_delim("ins.genes.bed", delim = "\t",
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

del_genes <- read_delim("del.genes.bed", delim = "\t",
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

tra_genes <- read_delim("tra.genes.bed", delim = "\t",
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

dup_genes <- read_delim("dup.genes.bed", delim = "\t",
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

inv_genes <- read_delim("inv.genes.bed", delim = "\t",
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

ins_exons <- read_delim("ins.exons.bed", delim = "\t",
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

del_exons <- read_delim("del.exons.bed", delim = "\t",
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

tra_exons <- read_delim("tra.exons.bed", delim = "\t",
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

dup_exons <- read_delim("dup.exons.bed", delim = "\t",
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

inv_exons <- read_delim("inv.exons.bed", delim = "\t",
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
    mutate(., sanger = if_else(sv_name %in% ins_s$sv_name, 
                               "N", "Y")) %>%
    mutate(., ensembl = if_else(sv_name %in% ins_e$sv_name, 
                                "N", "Y")) %>%
    left_join(., select(ins_genes, sv_name, gene)) %>%
    left_join(., select(ins_exons, sv_name, gene, exon_num) %>% 
              distinct() %>%
              arrange(., sv_name, gene, exon_num) %>%
              group_by(., sv_name, gene) %>% 
              summarize(., exons = paste(exon_num, collapse = ";")))

del_full <- del %>%
    mutate(., sv_type = "deletion") %>%
    mutate(., sanger = if_else(sv_name %in% del_s$sv_name, 
                               "N", "Y")) %>%
    mutate(., ensembl = if_else(sv_name %in% del_e$sv_name, 
                                "N", "Y")) %>%
    left_join(., select(del_genes, sv_name, gene)) %>%
    left_join(., select(del_exons, sv_name, gene, exon_num) %>% 
              distinct() %>%
              arrange(., sv_name, gene, exon_num) %>%
              group_by(., sv_name, gene) %>% 
              summarize(., exons = paste(exon_num, collapse = ";")))

tra_full <- tra %>%
    mutate(., sv_type = "translocation") %>%
    mutate(., sanger = NA) %>%
    mutate(., ensembl = if_else(sv_name %in% tra_e$sv_name, 
                                "N", "Y")) %>%
    left_join(., select(tra_genes, sv_name, gene)) %>%
    left_join(., select(tra_exons, sv_name, gene, exon_num) %>% 
              distinct() %>%
              arrange(., sv_name, gene, exon_num) %>%
              group_by(., sv_name, gene) %>% 
              summarize(., exons = paste(exon_num, collapse = ";")))

inv_full <- inv %>%
    mutate(., sv_type = "inversion") %>%
    mutate(., sanger = NA) %>%
    mutate(., ensembl = if_else(sv_name %in% tra_e$sv_name, 
                                "N", "Y")) %>%
    left_join(., select(inv_genes, sv_name, gene)) %>%
    left_join(., select(inv_exons, sv_name, gene, exon_num) %>% 
                  distinct() %>%
                  arrange(., sv_name, gene, exon_num) %>%
                  group_by(., sv_name, gene) %>% 
                  summarize(., exons = paste(exon_num, collapse = ";")))

dup_full <- dup %>%
    mutate(., sv_type = "duplication") %>%
    mutate(., sanger = NA) %>%
    mutate(., ensembl = if_else(sv_name %in% tra_e$sv_name, 
                                "N", "Y")) %>%
    left_join(., select(dup_genes, sv_name, gene)) %>%
    left_join(., select(dup_exons, sv_name, gene, exon_num) %>% 
                  distinct() %>%
                  arrange(., sv_name, gene, exon_num) %>%
                  group_by(., sv_name, gene) %>% 
                  summarize(., exons = paste(exon_num, collapse = ";")))

bind_rows(ins_full, del_full, tra_full, inv_full, dup_full) %>%
    mutate(., pos = start + 1) %>%
    mutate(., sv_size = abs(stop - start) + 1) %>%
    relocate(chr, pos, sv_name, sv_size, sv_type, 
             sanger, ensembl, gene, exons) %>%
    select(., -start, -stop) %>%
    write_csv(., "survivor_results.csv")