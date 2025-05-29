library(data.table)
library(tidyverse)
library(signature.tools.lib)
library(BSgenome.Hsapiens.UCSC.hg38)
# library(biomaRt)
library(GenomicRanges)

setwd("~/manuscript/Biopsy analysis of trial S1616/tables")
patients <- fread("table_s1.tsv")
samples <- fread("table_s2.tsv")
load("colors.Rda")

# Set factor levels
patients$arm <- factor(patients$arm, levels = c("Combination", "Ipilimumab"))
patients$resp <- factor(patients$resp, levels = c("CR","PR","SD","PD"))
patients$resp2 <- factor(ifelse(patients$resp %in% c('CR','PR'), 'CR/PR', 'SD/PD'), levels = c('CR/PR','SD/PD'))
samples$timepoint <- factor(samples$timepoint, levels = c('baseline','prior to pd1','ontx c2'))


## Get WES data
wes <- samples %>% left_join(patients) %>% filter(!is.na(wes.id))
nrow(wes)
wes %>% group_by(subject) %>% summarise %>% nrow
wes %>% group_by(subject, arm) %>% summarise %>% group_by(arm) %>% count

# Choose representative samples for wes
# Note that these manual filters were set after review
WESrep <- wes %>% mutate(percpos.s100 = ifelse(is.na(percpos.s100), -1, percpos.s100)) %>%
  arrange(subject, timepoint, desc(percpos.s100), desc(wes.meanX_coverage)) %>%
  filter(subject != "PT0689" | timepoint != "ontx c2") %>%
  filter(subject != "PT0741" | timepoint != "baseline") %>%
  filter(subject != "PT0621" | timepoint != "ontx c2") %>% 
  filter(subject != "PT0616" | ! timepoint %in% c("ontx c2","prior to pd1")) %>%
  filter(subject != "PT0692" | timepoint != "ontx c2") %>%
  group_by(subject) %>% slice_head(n = 1) %>%
  mutate(percpos.s100 = ifelse(percpos.s100 == -1, NA, percpos.s100))

load("~/manuscript/s1616_tm/data/20240801_variants.Rda")
load("~/manuscript/s1616_tm/data/20240801_cn.Rda")

WESvar <- WESvar %>% inner_join(wes)
WESseg <- WESseg %>% inner_join(wes)

# > gene_loc %>% filter(hgnc_symbol == "B2M")
# ensembl_gene_id hgnc_symbol chromosome_name start_position end_position
# 1 ENSG00000166710         B2M           chr15       44711358     44718851
SEG %>% filter(sample == "SC004-Pool2-PT0765-030520-BX-B4-DNA", chromosome == "chr15", start.pos<44711358, end.pos>44718851) %>%
  mutate(len = (end.pos-start.pos+1)/1e3)

## Tumor mutation burden
tmb <- VAR %>% filter(tumorVAF>0.05) %>% group_by(sample) %>% count(name = "nvar_0.05") %>%
  inner_join(WESvar) %>%
  rowwise %>% mutate(TMB = `nvar_0.05`/(wes.perc_50X*35694053/1e6)) %>% ungroup# The size of the targeted Exome
# check
table_s2 <- fread("table_s2.tsv")
table_s2 %>% select(specimen, TMB) %>% left_join(select(tmb, specimen, TMB), by = c("specimen")) %>% mutate(diff = TMB.x-TMB.y)

## Get canonical drivers
ckb <- fread("ckb_gene_variant_table.tsv")
ckb %>% group_by(Gene) %>% summarise %>% nrow
ckb <- ckb %>% mutate(Variant_ = ifelse(grepl("^(\\w+\\d+)(\\w+|\\*|\\?)\\b", Variant), gsub("^(\\w+\\d+)(\\w+|\\*|\\?)\\b", "\\1", Variant), NA))
ckb <- tibble(Gene = "NRAS", "Variant" = "GQ60-61GK", "Impact" = "missense", "Protein Effect" = "loss of function - predicted") %>% full_join(ckb)
ckb_kras <- ckb %>% filter(Gene == "NRAS") %>% mutate(Gene = "KRAS")
ckb_other <- c('B2M','JAK1')

# Perfect Gene/Mutation match for known function (KRAS is not included in the publicly available database, but known to be similar to HRAS/NRAS)
VARperfect <- VAR %>% filter(tumorVAF>0.05) %>% inner_join(ckb, by = c("SYMBOL" = "Gene", "Mutation" = "Variant")) %>%
  mutate(in_ckb = TRUE)
VARperfect_kras <- VAR %>% filter(tumorVAF>0.05) %>% inner_join(ckb_kras, by = c("SYMBOL" = "Gene", "Mutation" = "Variant")) %>%
  mutate(in_ckb = NA)
VARnf1 <- VAR %>% filter(SYMBOL == "NF1" & tumorVAF>0.05) %>% 
  mutate(`Protein Effect` = ifelse(IMPACT == "HIGH", "loss of function - predicted", "unknown"),
         in_ckb = NA)
VAR_ <- VARperfect %>% full_join(VARperfect_kras) %>% full_join(VARnf1)
# Add in anything outside of ckb but in corresponding genes
VAR_ <- VAR %>% filter(tumorVAF>0.05 & SYMBOL %in% union(ckb$Gene, c(ckb_kras$Gene, ckb_other))) %>% anti_join(VAR_) %>% 
  mutate(`Protein Effect` = ifelse(IMPACT == "HIGH", "loss of function - predicted", "unknown"), in_ckb = NA) %>%
  full_join(VAR_)

## Mutation signatures
hg38 <- BSgenome.Hsapiens.UCSC.hg38
VAR5VAF <- VAR %>% filter(tumorVAF>0.05) %>% inner_join(WESvar)
VAR5VAF <- split(VAR5VAF, VAR5VAF$sample)

# Annotate SNVs
snv_cats <- lapply(VAR5VAF, function(prep){
  tab_snv <- prep %>% filter(nchar(REF) == 1 & nchar(ALT) == 1) %>% dplyr::select(CHROM, POS, REF, ALT) %>% distinct()
  colnames(tab_snv) <- c('chr','position','REF','ALT')
  if(nrow(tab_snv)==0){
    return(NULL)
  } else {
    cat_snv <- tabToSNVcatalogue(tab_snv, genome.v = "hg38")
    cat_return <- cat_snv$catalogue %>% rownames_to_column("nt_change")
    return(cat_return)
  }
})
snv_cats <- snv_cats %>% enframe(name = "sample") %>% unnest(value)
dnv_cats <- lapply(VAR5VAF, function(prep){
  tab_dnv <- prep %>% filter(grepl("MNV", FILTER) & nchar(REF) == 2 & nchar(ALT) == 2) %>% dplyr::select(CHROM, POS, REF, ALT) %>% distinct()
  colnames(tab_dnv) <- c('chr','position','REF','ALT')
  if(nrow(tab_dnv)==0){
    return(NULL)
  } else {
    cat_dnv <- tabToDNVcatalogue(tab_dnv)
    cat_return <- cat_dnv$DNV_catalogue %>% rownames_to_column("nt_change") %>% rename("sample" = "catalogue")
    return(cat_return)
  }
})
dnv_cats <- dnv_cats %>% enframe(name = "sample") %>% unnest(value)
cats <- snv_cats %>% full_join(dnv_cats)
cats <- cats %>% select(nt_change) %>% distinct %>% rowwise %>% mutate(sample = list(names(VAR5VAF))) %>% unnest(sample) %>%
  full_join(cats) %>% mutate(catalogue = ifelse(is.na(catalogue), 0, catalogue))
signatures <- cats %>%
  mutate(status = case_when(nt_change == "CC>TT" ~ "muttype_CC>TT",
                            grepl("\\[C>T\\]", nt_change) ~ "muttype_C>T",
                            nt_change == "T[T>A]T" ~ "muttype_T[T>A]T (UV, SBS7c)",
                            TRUE ~ "muttype_Other")) %>%
  group_by(sample, status) %>% summarise(n = sum(catalogue)) %>% 
  group_by(sample) %>% mutate(n = n/sum(n)) %>% spread(status, n)
signatures_ <- WESvar %>% select(sample, wes.id) %>% left_join(signatures) %>% mutate(sample = NULL)

## Pull retrospective dataset
retro_ann <- data.table::fread("~/repos/MORRISON-1-public/WES/WES-CancerCell-MORRISON1-metadata.tsv")
retro_incl <- retro_ann %>% filter(grepl("-pre", timepoint.id),
                                   sample.tumor.type %in% c('cutaneous','mucosal','unknown'), 
                                   treatment.regimen.name %in% c('Combo','PD1'))
retro_var <- data.table::fread("~/repos/MORRISON-1-public/WES/data/WES-CancerCell-MORRISON1-variants.tsv")
# retro_var <- data.table::fread("~/Downloads/fat-variants-final.tsv")
retro_segs <- data.table::fread("~/Google Drive/Shared drives/Galvez-Campbell/Morrison- Processed Data/morrison_segments.tsv")

# TMB
nonsilent_consequences <- c('transcript_ablation','splice_acceptor_variant','splice_donor_variant','stop_gained',
                            'frameshift_variant','stop_lost','start_lost','transcript_amplification',
                            'inframe_insertion','inframe_deletion','missense_variant','protein_altering_variant')
morrison_wesqc <- data.table::fread("~/dat/shared/morrison/int_data/wes-qc-table.tsv") %>%
  dplyr::select(sample, depth)
#
morrison_ <- retro_var %>% left_join(retro_ann) %>% filter(collapsed.consequences %in% nonsilent_consequences)
morrison_ <- morrison_ %>% filter(vaf>0.05)
morrison_tmb <- morrison_ %>% 
  group_by(sample.id, subject.id, cohort, previous.treatment, treatment.regimen.name, timepoint.id, bor, tmb) %>%
  count %>% filter(grepl("-pre", timepoint.id))
# Coverage
coverage_dir <- "~/dat/shared/morrison/int_data/cumulative-depths/"
coverage_file <- function(sample) {
  paste(coverage_dir, sample, ".sample_cumulative_coverage_counts", sep="")
}
coverage_data <- function(sample) {
  d = read.delim(coverage_file(sample))
  d$sample.id <- sample
  d
}
# Morrison TMB
morrison_tmb <- morrison_tmb %>%
  rowwise %>%
  mutate(bases_50x = as.numeric(coverage_data(sample.id)['gte_50']),
         TMB = n*1e6/bases_50x)
# tmb_retro %>% full_join(old_script, by = c("sample.id"), suffix = c(".first", ".old_script")) %>% 
#   full_join(morrison_tmb, by = c("sample.id")) %>% 
#   select(sample.id, matches("TMB"), matches("n\\."), n)

priordatasets_ann <- morrison_tmb %>% left_join(morrison_wesqc, by = c("sample.id" = "sample")) %>%
  dplyr::select(sample.id, subject.id, cohort, previous.treatment, treatment.regimen.name, bases_50x, TMB, depth) %>%
  mutate(included_in_s1616comparison = ifelse(sample.id %in% retro_incl$sample.id, TRUE, FALSE))

RETROperfect <- retro_var %>% filter(vaf>0.05) %>% inner_join(ckb, by = c("gene.hgnc.symbol" = "Gene", "aa.change" = "Variant")) %>%
  mutate(in_ckb = TRUE)
RETROperfect_kras <- retro_var %>% filter(vaf>0.05) %>% inner_join(ckb_kras, by = c("gene.hgnc.symbol" = "Gene", "aa.change" = "Variant")) %>%
  mutate(in_ckb = NA)
RETROnf1 <- retro_var %>% filter(gene.hgnc.symbol == "NF1" & vaf>0.05) %>%
  mutate(`Protein Effect` = case_when(collapsed.consequences %in% c('splice_acceptor_variant','splice_donor_variant','stop_gained','frameshift_variant',
                                                                    'stop_lost','start_lost') ~ "loss of function - predicted",
                                      collapsed.consequences %in% c('inframe_insertion','inframe_deletion','missense_variant','protein_altering_variant') ~ "unknown",
                                      TRUE ~ as.character(NA)),
         in_ckb = NA) %>%
  filter(!is.na(`Protein Effect`))
RETRO_ <- RETROperfect %>% full_join(RETROperfect_kras) %>% full_join(RETROnf1)
#
RETRO_ <- retro_var %>% filter(vaf>0.05 & gene.hgnc.symbol %in% union(ckb$Gene, c(ckb_kras$Gene, ckb_other))) %>% anti_join(RETRO_) %>%
  mutate(`Protein Effect` = case_when(collapsed.consequences %in% c('splice_acceptor_variant','splice_donor_variant','stop_gained','frameshift_variant',
                                                                    'stop_lost','start_lost') ~ "loss of function - predicted",
                                      collapsed.consequences %in% c('inframe_insertion','inframe_deletion','missense_variant','protein_altering_variant') ~ "unknown",
                                      TRUE ~ as.character(NA)),
         in_ckb = NA) %>%
  filter(!is.na(`Protein Effect`)) %>%
  full_join(RETRO_)
# 
prep_segs <- retro_segs %>% filter(! (A==1 & B==1)) %>%
  mutate(status = case_when(CNt == 0 ~ "homozygous loss",
                            B == 0 ~ "loh",
                            CNt > 4 ~ "total cn>4",
                            TRUE ~ "other alteration")) %>%
  filter(status %in% c('homozygous loss','total cn>4'))
segs_gr <- GRanges(seqnames = prep_segs$chromosome, IRanges(start = prep_segs$start.pos, end = prep_segs$end.pos))
mcols(segs_gr) <- prep_segs[,c('CNt','A','B','file','status')]
#
ensembl <- biomaRt::useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="https://www.ensembl.org")
ensembl <- biomaRt::useDataset("hsapiens_gene_ensembl", mart=ensembl)
gene_loc <- biomaRt::getBM(attributes = c('ensembl_gene_id','hgnc_symbol','chromosome_name','start_position','end_position'),
                  filter = 'transcript_is_canonical', values = TRUE,
                  mart = ensembl)
gene_loc <- gene_loc %>% filter(chromosome_name %in% c(1:22,'X','Y')) %>% mutate(chromosome_name = paste0("chr", chromosome_name))
gene_gr <- GRanges(seqnames = gene_loc$chromosome_name, IRanges(start = gene_loc$start_position, end = gene_loc$end_position))
mcols(gene_gr) <- gene_loc[,c('ensembl_gene_id','hgnc_symbol')]
#
find_genes<- findOverlaps(gene_gr, segs_gr)
gene_segs <- gene_gr[queryHits(find_genes)]
mcols(gene_segs) <- c(mcols(gene_segs), as.data.frame(segs_gr[subjectHits(find_genes)]))
genes_cn <- data.frame(gene_segs) %>%
  filter(seqnames %in% paste0("chr", 1:22) & hgnc_symbol != "")
RETROCN_ <- genes_cn %>% filter(hgnc_symbol %in% union(ckb$Gene, c(ckb_kras$Gene, ckb_other)))


# Update Table S1
selected_wes <- WESrep %>% dplyr::select(subject, wes.id) %>% rename("wes.id" = "representative_wes")
patients <- patients %>% left_join(selected_wes)
# write.table(patients, file = "~/manuscript/Biopsy analysis of trial S1616/tables/table_s1.tsv", sep = "\t", row.names = F)

# Update Table S2
samples_ <- samples %>% left_join(select(tmb, specimen, nvar_0.05, TMB))
samples_ <- samples_ %>% left_join(signatures_)
# write.table(samples_, file = "~/manuscript/Biopsy analysis of trial S1616/tables/table_s2.tsv", sep = "\t", row.names = F)


# Export tables
# All variants with driver annotation
prep <- VAR %>% full_join(VAR_)
prep <- WESvar %>% dplyr::select(sample, wes.id) %>% inner_join(prep) %>% mutate(sample = NULL, name = NULL)
# write.table(prep, file = "~/manuscript/Biopsy analysis of trial S1616/tables/s1616_full_variants.tsv", sep = '\t', quote = F, row.names = F)


# All copy number statuses
prep <- WESseg %>% dplyr::select(sample, wes.id) %>% inner_join(GENECN) %>% mutate(sample = NULL)
# write.table(prep, file = "~/manuscript/Biopsy analysis of trial S1616/tables/s1616_genecn.tsv", sep = '\t', quote = F, row.names = F)

# Morrison updated TMB
# write.table(priordatasets_ann, file = "~/manuscript/Biopsy analysis of trial S1616/tables/priordatasets_annotation.tsv", sep = '\t', quote = F, row.names = F)

# Morrison variant annotation for genomic comparison
prep <- retro_var %>% full_join(RETRO_)
prep <- priordatasets_ann %>% ungroup %>% select(sample.id) %>% inner_join(prep)
# write.table(prep, file = "~/manuscript/Biopsy analysis of trial S1616/tables/priordatasets_drivervariants.tsv", sep = '\t', quote = F, row.names = F)

# Morrison copy number statuses
prep <- RETROCN_ %>% mutate(sample.id = gsub("_segments.txt", "", file), file = NULL)
prep <- priordatasets_ann %>% ungroup %>% select(sample.id) %>% inner_join(prep)
# write.table(prep, file = "~/manuscript/Biopsy analysis of trial S1616/tables/priordatasets_genecn.tsv", sep = '\t', quote = F, row.names = F)



##### END #####
