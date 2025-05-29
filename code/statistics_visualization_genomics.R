library(tidyverse)
library(ggpubr)
library(broom)
library(patchwork)

options(future.globals.maxSize = 24000 * 1024^2)

setwd("~/manuscript/Biopsy analysis of trial S1616/tables")
load("colors.Rda")
load("~/dat/shared/palettes.Rda")

# Load in patient and sample annotation
patients <- data.table::fread("table_s1.tsv")
samples <- data.table::fread("table_s2.tsv")
var <- data.table::fread("s1616_full_variants.tsv")
var <- var %>% mutate(tumorVAF = ifelse(is.na(tumorVAF), as.numeric(tumorAD)/as.numeric(tumorDP), as.numeric(tumorVAF)))
cn <- data.table::fread("s1616_genecn.tsv")
# Load in variant annotation information
ckb <- data.table::fread("ckb_gene_variant_table.tsv")
ckb_genes <- unique(c(unique(ckb$Gene), 'KRAS','B2M','JAK1','JAK2'))

# Set factor levels
patients$arm <- factor(patients$arm, levels = c("Combination", "Ipilimumab"))
patients$subtype <- factor(patients$subtype, levels = c("Cutaneous","Mucosal","Unknown"))
patients$prior_mapki <- factor(patients$prior_mapki, levels = c("Dabrafenib + Trametinib"))
patients$resp <- factor(patients$resp, levels = c("CR","PR","SD","PD"))
patients$resp2 <- factor(ifelse(patients$resp %in% c('CR','PR'), 'CR/PR', 'SD/PD'), levels = c('CR/PR','SD/PD'))
samples$timepoint <- factor(samples$timepoint, levels = c('baseline','prior to pd1','ontx c2'))

# Confirm sample numbers
patients %>% filter(!is.na(representative_wes)) %>%
  group_by(arm) %>% count
wes <- samples %>% filter(!is.na(wes.id)) %>% left_join(patients) # All wes samples
wes_rep <- wes %>% filter(wes.id == representative_wes) # Only representative wes samples (deduplicate at the patient level)
#
var_rep <- var %>% inner_join(wes_rep)
# Manual fixes: # SC005-Pool2-PT0563-081417-BXA2-DNA MNV ENST00000369535.5:c.180_181delinsTA results in Q61K Mutation
var_rep <- var_rep %>%
  mutate(`Protein Effect` = ifelse(HGVSc == "ENST00000369535.5:c.180_181delinsTA", "loss of function - predicted", `Protein Effect`))
var_ann <- tibble(`Protein Effect` = c('gain of function','gain of function - predicted','loss of function','loss of function - predicted',
                                           'no effect','no effect - predicted','unknown'),
                      description = c('GOF','GOF','LOF','LOF',NA,NA,NA))
var_rep <- var_rep %>% left_join(var_ann)
#
cn_rep <- cn %>% inner_join(wes_rep)
cn_ann <- data.frame(status = c('homozygous loss','loh','total cn>4'),
                     description = c('Deep del','LOH','CN>4'))
cn_rep <- cn_rep %>% left_join(cn_ann)

#### Statistics ####
glm_res <- data.table::fread("~/manuscript/Biopsy analysis of trial S1616/tables/table_s3.tsv")

# Compare TMB, mutational signatures, and drivers between arms
compare_means(TMB ~ arm, wes_rep, method = 'wilcox')
wes_ <- wes_rep %>% gather(muttype, value, matches('muttype_')) %>%
  mutate(group = ifelse(arm == "Ipilimumab", "Ipi", ifelse(resp2 == "CR/PR", "Comb. CR/PR", "Comb. SD/PD")))
wes_$group <- factor(wes_$group, levels = c("Comb. CR/PR", "Comb. SD/PD", "Ipi"))
arm_comp <- compare_means(value ~ arm, wes_, method = 'wilcox', group.by = 'muttype') %>%
  mutate(y = case_when(muttype == "muttype_C>T" ~ 0.9,
                       muttype == "muttype_CC>TT" ~ 0.125,
                       muttype == "muttype_Other" ~ 1.10,
                       muttype == "muttype_T[T>A]T (UV, SBS7c)" ~ 0.03),
         ysub = y-0.025*y,
         ysubdir = ysub-0.025*y)
resp_comp <- compare_means(value ~ group, wes_, method = 'wilcox', group.by = 'muttype') %>%
  filter(grepl("Comb", group1) & grepl("Comb", group2)) %>%
  mutate(y = case_when(muttype == "muttype_C>T" ~ 0.9,
                       muttype == "muttype_CC>TT" ~ 0.125,
                       muttype == "muttype_Other" ~ 1.10,
                       muttype == "muttype_T[T>A]T (UV, SBS7c)" ~ 0.03),
         ysub = y-0.025*y,
         ysubdir = ysub-0.025*y,
         diff_y = 0.05*y,
         y = y+diff_y, ysub = ysub+diff_y)
ggplot(wes_, aes(x = group, y = value, fill = group)) +
  geom_boxplot(fill = NA, outlier.color = NA, size = 0.25) + 
  ggbeeswarm::geom_beeswarm(size = 1, shape = 21, colour = 'black', stroke = 0.1) +
  #
  geom_text(data = resp_comp, aes(x = 1.5, y = y+0.05*y, label = p.format), inherit.aes = FALSE, size = 5*0.36) +
  geom_segment(data = resp_comp, aes(y = y, yend = y, x = 1, xend = 2), inherit.aes = FALSE, size = 0.25) +
  geom_segment(data = resp_comp, aes(y = y, yend = ysub, x = 1, xend = 1), inherit.aes = FALSE, size = 0.25) +
  geom_segment(data = resp_comp, aes(y = y, yend = ysub, x = 2, xend = 2), inherit.aes = FALSE, size = 0.25) +
  #
  geom_text(data = arm_comp, aes(x = 2.25, y = y+0.05*y, label = p.format), inherit.aes = FALSE, size = 5*0.36) +
  geom_segment(data = arm_comp, aes(y = y, yend = y, x = 1.5, xend = 3), inherit.aes = FALSE, size = 0.25) +
  geom_segment(data = arm_comp, aes(y = y, yend = ysub, x = 1.5, xend = 1.5), inherit.aes = FALSE, size = 0.25) +
  geom_segment(data = arm_comp, aes(y = y, yend = ysub, x = 3, xend = 3), inherit.aes = FALSE, size = 0.25) +
  geom_segment(data = arm_comp, aes(y = ysub, yend = ysub, x = 1, xend = 2), inherit.aes = FALSE, size = 0.25) +
  geom_segment(data = arm_comp, aes(y = ysubdir, yend = ysub, x = 1, xend = 1), inherit.aes = FALSE, size = 0.25) +
  geom_segment(data = arm_comp, aes(y = ysubdir, yend = ysub, x = 2, xend = 2), inherit.aes = FALSE, size = 0.25) +
  #
  scale_fill_manual(values = as.character(c(resp2_colors, arm_colors$Ipilimumab))) +
  facet_wrap(~ gsub("muttype_", "", muttype), scales = 'free_y', nrow = 2) +
  labs(y = "Fraction of SNVs/DNVs", x = "Patient subset") +
  theme_bw() +
  theme(text = element_text(size = 5), title = element_text(size = 7), legend.position = 'none',
        axis.ticks = element_line(size = 0.1),
        panel.border = element_rect(linewidth = 0.1, colour = 'black'),
        strip.background = element_blank(), strip.clip = 'off',
        panel.spacing = unit(2, 'pt'), plot.margin = unit(c(0,0,0,0), 'pt'))
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/extfigure2c.pdf", height = 2.5, width = 3.75, units = 'in', family = 'ArialMT')

xlabs <- wes_rep %>% group_by(arm) %>% count %>%
  mutate(label = paste0(str_to_title(arm), "\n(N=",n,")"))
ggplot(wes_rep, aes(x = arm, y = TMB, fill = arm, colour = arm)) +
  geom_boxplot(fill = NA, outlier.color = NA, size = 0.25) + 
  ggbeeswarm::geom_beeswarm(size = 1, shape = 21, colour = 'black', stroke = 0.1) +
  stat_compare_means(comparisons = list(c('Combination','Ipilimumab')), size = 5*0.36) +
  scale_fill_manual(values = as.character(arm_colors)) +
  scale_colour_manual(values = as.character(arm_colors)) +
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  scale_x_discrete(labels = xlabs$label) +
  labs(x = "Patient subset", y = "TMB (Mut/Mb)") +
  theme_bw() +
  theme(text = element_text(size = 5), title = element_text(size = 7), legend.position = 'none',
        axis.ticks = element_line(size = 0.1),
        panel.border = element_rect(linewidth = 0.1, colour = 'black'))
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/extfigure2b.pdf", height = 2.5, width = 1.25, units = 'in', family = 'ArialMT')

# # GLM: gene mutation rates across all genes
# prep_glm <- wes_rep %>% group_by(arm, subject) %>% summarise
# prep_glm <- var_rep %>% filter(tumorVAF>0.05) %>% group_by(SYMBOL) %>% summarise %>% cross_join(prep_glm)
# prep_glm %>% group_by(SYMBOL) %>% summarise %>% nrow
# # 17,116 genes mutated in at least one sample
# prep_glm <- var_rep %>% filter(tumorVAF>0.05) %>% group_by(arm, subject, SYMBOL) %>% summarise(mut = 1) %>%
#   full_join(prep_glm) %>% mutate(mut = ifelse(is.na(mut), 0, 1))
# # Gene must be mutated in at least 5% of patients
# prep_glm <- prep_glm %>% group_by(SYMBOL) %>% summarise(mut = sum(mut)) %>% filter(mut>nrow(wes_rep)*0.05) %>%
#   select(SYMBOL) %>% inner_join(prep_glm)
# prep_glm %>% group_by(SYMBOL) %>% summarise %>% nrow
# # 9,025 genes mutated in at least 4 samples (>5% of patients)
# prep_glm <- prep_glm %>% mutate(arm_ = ifelse(arm == "Combination", 1, 0))
# run_glm <- prep_glm %>% ungroup %>% nest_by(SYMBOL) %>%
#   mutate(fitMut = list(glm(arm_ ~ mut, data = data, family = binomial()))) %>%
#   summarize(tidy(fitMut))
# glm <- run_glm %>% filter(term == 'mut') %>%
#   ungroup %>%
#   mutate(q.gene.mut = p.adjust(p.value, method='BH')) %>%
#   arrange(q.gene.mut)
# # Save all glm objects
# glm_res <- glm %>% mutate(comparison = "combination_v_ipilimumab")

# # GLM: known genetic alteration rates across all genes
# prep_glm <- wes_rep %>% group_by(arm, subject) %>% summarise
# prep_glm <- var_rep %>% filter(tumorVAF>0.05) %>% 
#   full_join(cn_rep, by = c(intersect(colnames(var_rep), colnames(cn_rep)), "SYMBOL" = "hgnc_symbol")) %>%
#   filter(!is.na(description)) %>% unite(altered_gene, SYMBOL, description, sep = " ", remove = FALSE) %>%
#   group_by(altered_gene) %>% summarise %>% cross_join(prep_glm)
# prep_glm %>% group_by(altered_gene) %>% summarise %>% nrow
# # 87,795 alterations mutated in at least one sample
# prep_glm <- var_rep %>% filter(tumorVAF>0.05) %>% 
#   full_join(cn_rep, by = c(intersect(colnames(var_rep), colnames(cn_rep)), "SYMBOL" = "hgnc_symbol")) %>%
#   filter(!is.na(description)) %>% unite(altered_gene, SYMBOL, description, sep = " ", remove = FALSE) %>%
#   group_by(arm, subject, altered_gene) %>% summarise(gene_alteration = 1) %>%
#   full_join(prep_glm) %>% mutate(gene_alteration = ifelse(is.na(gene_alteration), 0, 1))
# # Gene must be mutated in at least 5% of patients
# prep_glm <- prep_glm %>% group_by(altered_gene) %>% summarise(gene_alteration = sum(gene_alteration)) %>% filter(gene_alteration>nrow(wes_rep)*0.05) %>%
#   select(altered_gene) %>% inner_join(prep_glm)
# prep_glm %>% group_by(altered_gene) %>% summarise %>% nrow
# # 52,407 genes mutated in at least 4 samples (>5% of patients)
# prep_glm <- prep_glm %>% mutate(arm_ = ifelse(arm == "Combination", 1, 0))
# run_glm <- prep_glm %>% ungroup %>% nest_by(altered_gene) %>%
#   mutate(fitMut = list(glm(arm_ ~ gene_alteration, data = data, family = binomial()))) %>%
#   summarize(tidy(fitMut))
# glm <- run_glm %>% filter(term == 'gene_alteration') %>%
#   ungroup %>%
#   mutate(q.gene.mut = p.adjust(p.value, method='BH')) %>%
#   arrange(q.gene.mut)
# # Save to full glm results
# glm_res <- glm %>% mutate(comparison = "combination_v_ipilimumab") %>%
#   full_join(glm_res)

## Load prior dataset
prior <- data.table::fread("priordatasets_annotation.tsv")
prior <- prior %>% filter(included_in_s1616comparison) %>% mutate(dataset = "prior") #%>%
  # filter(previous.treatment == "naive", treatment.regimen.name == "PD1")
comp_prior <- wes_rep %>% mutate(dataset = arm) %>% 
  mutate(resp2 = ifelse(resp %in% c('CR','PR'), 'CR/PR', 'SD/PD')) %>% 
  mutate(dataset_ = ifelse(dataset == "prior", dataset, ifelse(arm == "Combination", resp2, "Ipi")))
comp_prior$dataset_ <- factor(comp_prior$dataset_, levels = c('prior','CR/PR','SD/PD','Ipi'))
comp_prior %>% group_by(dataset) %>% summarise(mean = mean(TMB), median = median(TMB))
comp_prior %>% group_by(dataset_) %>% summarise(mean = mean(TMB), median = median(TMB))
stat_tmbprior <- compare_means(TMB ~ dataset, comp_prior, method = 'wilcox')
stat_tmbbor <- compare_means(TMB ~ resp2, comp_prior, method = 'wilcox', group.by = c('dataset', 'arm'))
xlabs <- comp_prior %>% 
  group_by(dataset_) %>% count %>% ungroup %>% 
  mutate(label = ifelse(dataset_ == "prior", paste0(str_to_title(dataset_),"\n(N=",n,")"), as.character(dataset_))) %>%
  mutate(label = ifelse(dataset_ == "prior", label, 
                        paste0("S1616 ", 
                               ifelse(dataset_ == "Ipi", "Ipi\n", "Comb.\n"), 
                               ifelse(dataset_ == "Ipi", "", as.character(dataset_)), 
                               " (N=",n,")")))
ggplot(comp_prior, aes(x = dataset_, y = TMB, fill = dataset_)) +
  geom_boxplot(fill = NA, outlier.color = NA, size = 0.25) + 
  ggbeeswarm::geom_beeswarm(shape = 21, size = 1, stroke = 0.1) +
  #
  geom_segment(size = 0.10, data = stat_tmbprior, aes(x = 1, xend = 2, y = 800, yend = 800), inherit.aes = FALSE) +
  geom_segment(size = 0.10, data = tibble(x = c(1,2), y = c(700,700), yend = c(800, 800)), lineend = "square",
               aes(x = x, xend = x, y = y, yend = yend), inherit.aes = FALSE) +
  geom_text(data = filter(stat_tmbbor, dataset == 'Combination', arm == "Combination"), 
            aes(x = 1.5, y = 1050, label = p.format), inherit.aes = FALSE, size = 5*0.36) +
  #
  geom_segment(size = 0.10, data = stat_tmbprior, aes(x = 1.5, xend = 3, y = 1500, yend = 1500), inherit.aes = FALSE) +
  geom_segment(size = 0.10, data = stat_tmbprior, aes(x = 1, xend = 2, y = 1400, yend = 1400), inherit.aes = FALSE) +
  geom_segment(size = 0.10, data = tibble(x = c(1.5, 3), y = c(1400,1400), yend = c(1500, 1500)), lineend = "square",
               aes(x = x, xend = x, y = y, yend = yend), inherit.aes = FALSE) +
  geom_segment(size = 0.10, data = tibble(x = c(1, 2), y = c(1400,1400), yend = c(1300, 1300)), lineend = "square",
               aes(x = x, xend = x, y = y, yend = yend), inherit.aes = FALSE) +
  geom_text(data = stat_tmbprior, aes(x = 2.25, y = 1850, label = p.format), inherit.aes = FALSE, size = 5*0.36) +
  #
  scale_fill_manual(values = c(as.character(resp2_colors), as.character(arm_colors$Ipilimumab))) +
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  scale_x_discrete(labels = xlabs$label) +
  labs(x = "Patient subset", y = "TMB (Mut/Mb)") +
  theme_bw() +
  theme(text = element_text(size = 5), title = element_text(size = 7), legend.position = 'none',
        axis.ticks = element_line(size = 0.1),
        panel.border = element_rect(linewidth = 0.1, colour = 'black'))
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figure1B.pdf", height = 2.5, width = 2.5, units = 'in', family = 'ArialMT')

comp_prior <- wes_rep %>% mutate(dataset = "s1616") %>% 
  full_join(prior, by = c("wes.id" = "sample.id", "dataset", "TMB", "wes.meanX_coverage" = "depth", "resp" = "bor")) %>%
  mutate(resp2 = ifelse(resp %in% c('CR','PR'), 'CR/PR', 'SD/PD')) %>% 
  mutate(dataset_ = ifelse(dataset == "prior", dataset, "S1616"))
comp_prior$dataset_ <- factor(comp_prior$dataset_, levels = c('prior','S1616'))
comp_prior %>% group_by(dataset_) %>% summarise(mean = mean(TMB), median = median(TMB))
stat_tmbprior <- compare_means(TMB ~ dataset, comp_prior, method = 'wilcox')
xlabs <- comp_prior %>% 
  group_by(dataset_) %>% count %>% ungroup %>% 
  mutate(label = ifelse(dataset_ == "prior", paste0(str_to_title(dataset_),"\n(N=",n,")"), as.character(dataset_))) %>%
  mutate(label = ifelse(dataset_ == "prior", label, 
                        paste0("S1616 ", 
                               ifelse(dataset_ == "Ipi", "Ipi\n", "Comb.\n"), 
                               ifelse(dataset_ == "Ipi", "", as.character(dataset_)), 
                               " (N=",n,")")))
ggplot(comp_prior, aes(x = dataset_, y = TMB, fill = dataset_)) +
  geom_boxplot(fill = NA, outlier.color = NA, size = 0.25) + 
  ggbeeswarm::geom_beeswarm(shape = 21, size = 1, stroke = 0.1) +
  #
  geom_segment(size = 0.10, data = stat_tmbprior, aes(x = 1, xend = 2, y = 800, yend = 800), inherit.aes = FALSE) +
  geom_segment(size = 0.10, data = tibble(x = c(1, 2), y = c(700,700), yend = c(800, 800)), lineend = "square",
               aes(x = x, xend = x, y = y, yend = yend), inherit.aes = FALSE) +
  geom_text(data = stat_tmbprior,
            aes(x = 1.5, y = 1050, label = p.format), inherit.aes = FALSE, size = 5*0.36) +
  #
  scale_fill_manual(values = c('white', 'grey50')) +
  scale_y_log10(labels = scales::label_number(drop0trailing = TRUE)) +
  scale_x_discrete(labels = xlabs$label) +
  labs(x = "Cohort", y = "TMB (Mut/Mb)") +
  theme_bw() +
  theme(text = element_text(size = 5), title = element_text(size = 7), legend.position = 'none',
        axis.ticks = element_line(size = 0.1),
        panel.border = element_rect(linewidth = 0.1, colour = 'black'))
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/extfigure2b.pdf", height = 2.5, width = 1.25, units = 'in', family = 'ArialMT')


## Compare alteration rates to prior dataset
prior_var <- data.table::fread("priordatasets_drivervariants.tsv")
prior_var <- prior_var %>% left_join(var_ann)
prior_cn <- data.table::fread("priordatasets_genecn.tsv")
prior_cn <- prior_cn %>% left_join(cn_ann)

# s1616 alteration rates
s1616_altrate <- wes_rep %>% group_by(arm, subject) %>% summarise
s1616_altrate <- var_rep %>% filter(tumorVAF>0.05) %>% 
  full_join(cn_rep, by = c(intersect(colnames(var_rep), colnames(cn_rep)), "SYMBOL" = "hgnc_symbol")) %>%
  filter(!is.na(description) & description != "LOH") %>% unite(altered_gene, SYMBOL, description, sep = " ", remove = FALSE) %>%
  group_by(SYMBOL, altered_gene) %>% summarise %>% cross_join(s1616_altrate)
s1616_altrate <- var_rep %>% filter(tumorVAF>0.05) %>% 
  full_join(cn_rep, by = c(intersect(colnames(var_rep), colnames(cn_rep)), "SYMBOL" = "hgnc_symbol")) %>%
  filter(!is.na(description) & description != "LOH") %>% unite(altered_gene, SYMBOL, description, sep = " ", remove = FALSE) %>%
  group_by(arm, subject, SYMBOL, altered_gene) %>% summarise(gene_alteration = 1) %>%
  full_join(s1616_altrate) %>% mutate(gene_alteration = ifelse(is.na(gene_alteration), 0, 1))
s1616_altrate <- s1616_altrate %>% group_by(SYMBOL, altered_gene) %>% summarise(gene_alteration = sum(gene_alteration)) %>% 
  mutate(alteration_freq = gene_alteration/nrow(wes_rep))

# prior alteration rates
prior_altrate <- prior %>% group_by(sample.id, subject.id) %>% summarise
prior_altrate <- prior_var %>% filter(vaf>0.05) %>% 
  full_join(prior_cn, by = c(intersect(colnames(prior_var), colnames(prior_cn)), "gene.hgnc.symbol" = "hgnc_symbol")) %>%
  filter(!is.na(description) & description != "LOH") %>% unite(altered_gene, gene.hgnc.symbol, description, sep = " ", remove = FALSE) %>%
  group_by(gene.hgnc.symbol, altered_gene) %>% summarise %>% cross_join(prior_altrate)
prior_altrate <- prior_var %>% filter(vaf>0.05) %>% 
  full_join(prior_cn, by = c(intersect(colnames(prior_var), colnames(prior_cn)), "gene.hgnc.symbol" = "hgnc_symbol")) %>%
  filter(!is.na(description) & description != "LOH") %>% unite(altered_gene, gene.hgnc.symbol, description, sep = " ", remove = FALSE) %>%
  group_by(sample.id, subject.id, gene.hgnc.symbol, altered_gene) %>% summarise(gene_alteration = 1) %>%
  full_join(prior_altrate) %>% mutate(gene_alteration = ifelse(is.na(gene_alteration), 0, 1))
prior_altrate <- prior_altrate %>% group_by(gene.hgnc.symbol, altered_gene) %>% summarise(gene_alteration = sum(gene_alteration)) %>% 
  mutate(alteration_freq = gene_alteration/nrow(prior))

# compare alteration rates
comp_altrate <- s1616_altrate %>% full_join(prior_altrate, by = c("SYMBOL" = "gene.hgnc.symbol", "altered_gene"), suffix = c(".s1616", ".prior"))
comp_altrate[is.na(comp_altrate)] <- 0
comp_altrate %>% filter(SYMBOL %in% ckb_genes) %>% 
  mutate(diff = abs(alteration_freq.s1616-alteration_freq.prior)) %>% 
  arrange(desc(diff))

comp_altrate %>% filter(SYMBOL %in% c(ckb_genes, 'NF1')) %>%
  mutate(color = ifelse(SYMBOL %in% c('BRAF','NRAS','NF1'), "A", "B")) %>%
  arrange(desc(color), -alteration_freq.s1616, -alteration_freq.prior) %>%
  ggplot(., aes(x = alteration_freq.prior, y = alteration_freq.s1616, colour = color)) +
  geom_abline(linetype = 'dashed', colour = 'grey80') +
  geom_point(size = 1) +
  ggrepel::geom_text_repel(aes(label = altered_gene), size = 5*0.36, max.overlaps = 25, force = 5) +
  coord_fixed(xlim = c(0, 0.4), ylim = c(0, 0.4)) +
  scale_colour_manual(values = c('red','grey20')) +
  labs(x = paste0("Frac. altered in prior dataset (N=", nrow(prior), ")"),
       y = paste0("Frac. altered in S1616 dataset (N=", nrow(wes_rep), ")")) +
  theme_bw() +
  theme(text = element_text(size = 5), title = element_text(size = 7),
        legend.position = 'none',
        axis.ticks = element_line(size = 0.1),
        panel.border = element_rect(linewidth = 0.1, colour = 'black'))
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figure1C.pdf", height = 2.5, width = 2.5, units = 'in', family = 'ArialMT')

# Compare TMB, mutational signatures, and drivers between response to combination
compare_means(TMB ~ resp2, wes_rep, method = 'wilcox', group.by = 'arm')
ggplot(wes_rep, aes(x = resp2, y = TMB)) + geom_boxplot() + scale_y_log10() + facet_wrap(~ arm)
wes_rep %>% group_by(arm, resp2) %>% summarise(median = median(TMB), n = n())
wes_ <- wes_rep %>% gather(muttype, value, matches('muttype_'))
compare_means(value ~ resp2, wes_, method = 'wilcox', group.by = c('muttype', 'arm')) %>%
  filter(arm == "Combination")

# GLM: gene mutation rates across all genes by response
wes_rep_ <- wes_rep %>% filter(arm == "Combination")
wes_rep_ %>% group_by(resp2) %>% count
# prep_glm <- wes_rep_ %>% filter(arm == "Combination") %>% group_by(arm, resp2, subject) %>% summarise
# prep_glm <- var_rep %>% filter(arm == "Combination") %>% filter(tumorVAF>0.05) %>% group_by(SYMBOL) %>% summarise %>% cross_join(prep_glm)
# prep_glm %>% group_by(SYMBOL) %>% summarise %>% nrow
# # 16,214 genes mutated in at least one sample
# prep_glm <- var_rep %>% filter(arm == "Combination") %>% filter(tumorVAF>0.05) %>% group_by(arm, resp2, subject, SYMBOL) %>% summarise(mut = 1) %>%
#   full_join(prep_glm) %>% mutate(mut = ifelse(is.na(mut), 0, 1))
# # Gene must be mutated in at least 5% of patients (at least 3)
# prep_glm <- prep_glm %>% group_by(SYMBOL) %>% summarise(mut = sum(mut)) %>% filter(mut>nrow(wes_rep_)*0.05) %>%
#   select(SYMBOL) %>% inner_join(prep_glm)
# prep_glm %>% group_by(SYMBOL) %>% summarise %>% nrow
# # 9,094 genes mutated in at least 3 samples (>5% of patients)
# prep_glm <- prep_glm %>% mutate(resp_ = ifelse(resp2 == "CR/PR", 1, 0))
# run_glm <- prep_glm %>% ungroup %>% nest_by(SYMBOL) %>%
#   mutate(fitMut = list(glm(resp_ ~ mut, data = data, family = binomial()))) %>%
#   summarize(tidy(fitMut))
# glm <- run_glm %>% filter(term == 'mut') %>%
#   ungroup %>%
#   mutate(q.gene.mut = p.adjust(p.value, method='BH')) %>%
#   arrange(q.gene.mut)
# # Save all glm objects
# glm_res <- glm %>% mutate(comparison = "combinationCRPR_v_combinationSDPD") %>%
#   full_join(glm_res)

# # GLM: known genetic alteration rates across all genes
# prep_glm <- wes_rep_ %>% filter(arm == "Combination") %>% group_by(arm, resp2, subject) %>% summarise
# prep_glm <- var_rep %>% filter(arm == "Combination") %>% filter(tumorVAF>0.05) %>% 
#   full_join(cn_rep, by = c(intersect(colnames(var_rep), colnames(cn_rep)), "SYMBOL" = "hgnc_symbol")) %>%
#   filter(!is.na(description)) %>% unite(altered_gene, SYMBOL, description, sep = " ", remove = FALSE) %>%
#   group_by(altered_gene, SYMBOL, description) %>% summarise %>% cross_join(prep_glm)
# prep_glm %>% group_by(altered_gene, SYMBOL, description) %>% summarise %>% nrow
# # 87,789 alterations mutated in at least one sample
# prep_glm <- var_rep %>% filter(arm == "Combination") %>% filter(tumorVAF>0.05) %>% 
#   full_join(cn_rep, by = c(intersect(colnames(var_rep), colnames(cn_rep)), "SYMBOL" = "hgnc_symbol")) %>%
#   filter(!is.na(description)) %>% unite(altered_gene, SYMBOL, description, sep = " ", remove = FALSE) %>%
#   group_by(arm, resp2, subject, altered_gene, SYMBOL, description) %>% summarise(gene_alteration = 1) %>%
#   full_join(prep_glm) %>% mutate(gene_alteration = ifelse(is.na(gene_alteration), 0, 1))
# # Gene must be mutated in at least 5% of patients
# prep_glm <- prep_glm %>% group_by(altered_gene) %>% summarise(gene_alteration = sum(gene_alteration)) %>% 
#   filter(gene_alteration>nrow(wes_rep_)*0.05) %>%
#   select(altered_gene) %>% inner_join(prep_glm)
# prep_glm %>% group_by(altered_gene) %>% summarise %>% nrow
# # 59,775 genes mutated in at least 4 samples (>5% of patients)
# prep_glm <- prep_glm %>% mutate(resp_ = ifelse(resp2 == "CR/PR", 1, 0))
# run_glm <- prep_glm %>% ungroup %>% nest_by(altered_gene) %>%
#   mutate(fitMut = list(glm(resp_ ~ gene_alteration, data = data, family = binomial()))) %>%
#   summarize(tidy(fitMut))
# glm <- run_glm %>% filter(term == 'gene_alteration') %>%
#   ungroup %>%
#   mutate(q.gene.mut = p.adjust(p.value, method='BH')) %>%
#   arrange(q.gene.mut)
# glm %>% filter(altered_gene %in% c('NRAS LOF','BRAF GOF','NF1 LOF','KIT GOF'))
# prep_glm %>% filter(altered_gene == "BRAF GOF", gene_alteration == 1) %>%
#   left_join(patients) %>% group_by(resp) %>% count
# prep_glm %>% filter(altered_gene == "BRAF GOF", gene_alteration == 1) %>%
#   left_join(patients) %>% group_by(prior_mapki) %>% count
# prep_glm %>% filter(altered_gene %in% c('NRAS LOF','BRAF GOF','NF1 LOF','KIT GOF')) %>% 
#   group_by(altered_gene, resp2, gene_alteration) %>% count %>% 
#   spread(gene_alteration, n) %>% mutate(frac = `1`/(`1`+`0`))

var_rep %>% filter(SYMBOL == "BRAF", description == "GOF") %>% group_by(subject, arm, resp, prior_mapki) %>% count %>%
  group_by(arm, prior_mapki) %>% count
# Save to full glm results
# glm_res <- glm %>% mutate(comparison = "combinationCRPR_v_combinationSDPD") %>%
#   full_join(glm_res)

glm_res %>% group_by(term, comparison) %>% count
# glm_res %>% select(comparison, term, SYMBOL, altered_gene, estimate, std.error, statistic, p.value, q.gene.mut) %>%
#   mutate(SYMBOL = ifelse(is.na(SYMBOL), gsub("^(\\w+.+) [G|L|D|C].+\\b", "\\1", altered_gene), SYMBOL)) %>%
#   write.table("~/manuscript/Biopsy analysis of trial S1616/tables/table_s3.tsv", sep = '\t', quote = F, row.names = F)

plot_glm <- glm_res %>% filter(term == "gene_alteration", comparison == "combinationCRPR_v_combinationSDPD") %>% 
  filter(! grepl("LOH", altered_gene)) %>%
  mutate(SYMBOL = ifelse(is.na(SYMBOL), gsub("^(\\w+.+) [G|L|D|C].+\\b", "\\1", altered_gene), SYMBOL)) %>% 
  mutate(marker = ifelse(SYMBOL %in% c(ckb_genes, 'NF1'), TRUE, FALSE), 
         colour = ifelse(marker & estimate>0, "A", ifelse(marker & estimate<0, "B", "C"))) %>%
  filter(abs(estimate)<3) %>%
  group_by(estimate, p.value, marker) %>% slice_head(n = 1) %>% arrange(marker, p.value)
ggplot(plot_glm, aes(x = estimate, y = -log10(p.value), colour = colour, fill = colour, alpha = colour)) +
  geom_hline(linetype = 'dashed', colour = 'grey50', size = 0.1, yintercept = -log10(0.05)) +
  geom_point(size = 1, shape = 21, stroke = 0.1) +
  ggrepel::geom_text_repel(data = filter(plot_glm, marker),
                           aes(label = altered_gene), size = 5*0.36, max.overlaps = 10, colour = 'black', show.legend = FALSE,
                           min.segment.length = 0.1) +
  scale_fill_manual(values = c(as.character(resp2_colors), 'white'),
                    labels = c('Driver associated with CR/PR','Driver associated with SD/PD','Other alterations')) +
  scale_colour_manual(values = c('black','black','black'),
                      labels = c('Driver associated with CR/PR','Driver associated with SD/PD','Other alterations')) +
  scale_alpha_manual(values = c(1, 1, 0.8),
                     labels = c('Driver associated with CR/PR','Driver associated with SD/PD','Other alterations')) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 1.5)) +
  labs(y = "-log10(p value)", x = "Estimate\n(Alteration v. response to Combination)",
       colour = "Direction", fill = "Direction", alpha = "Direction") +
  guides(colour = guide_legend(nrow = 3), alpha = guide_legend(nrow = 3),
         fill = guide_legend(nrow = 3)) +
  theme_bw() +
  theme(text = element_text(size = 5), title = element_text(size = 7),
        axis.ticks = element_line(size = 0.1),
        panel.border = element_rect(linewidth = 0.1, colour = 'black'),
        legend.position = 'bottom', legend.key.height = unit(0.25, 'line'),
        legend.key.spacing = unit(1, 'pt'), legend.margin = margin(0, 0, 0, 0))
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figure1D.pdf", height = 2.5, width = 2.375, unit = 'in', family="ArialMT")

# plot_glm <- glm_res %>% filter(term == "mut", comparison == "combination_v_ipilimumab") %>% 
#   mutate(SYMBOL = ifelse(is.na(SYMBOL), gsub("^(\\w+.+) [G|L|D|C].+\\b", "\\1", altered_gene), SYMBOL)) %>%
#   mutate(marker = ifelse(SYMBOL %in% c(ckb_genes, 'NF1'), TRUE, FALSE),
#          colour = ifelse(marker & estimate>0, "A", ifelse(marker & estimate<0, "B", "C"))) %>%
#   # filter(abs(estimate)<3) %>%
#   group_by(estimate, p.value, marker) %>% slice_head(n = 1) %>% arrange(marker, p.value)
# ggplot(plot_glm, aes(x = estimate, y = -log10(p.value), colour = colour, fill = colour, alpha = colour)) +
#   geom_hline(linetype = 'dashed', colour = 'grey50', size = 0.1, yintercept = -log10(0.05)) +
#   geom_point(shape = 21, stroke = 0.1) +
#   ggrepel::geom_text_repel(data = filter(plot_glm, marker),
#                            aes(label = SYMBOL), size = 5*0.36, max.overlaps = 10, colour = 'black', show.legend = FALSE,
#                            min.segment.length = 0.1) +
#   scale_fill_manual(values = c(as.character(arm_colors), 'white'),
#                     labels = c('Mutant gene associated with Combination arm','Mutant gene associated with Ipilimumab arm','Other genes')) +
#   scale_colour_manual(values = c('black','black','black'),
#                       labels = c('Mutant gene associated with Combination arm','Mutant gene associated with Ipilimumab arm','Other genes')) +
#   scale_alpha_manual(values = c(1, 1, 0.8),
#                      labels = c('Mutant gene associated with Combination arm','Mutant gene associated with Ipilimumab arm','Other genes')) +
#   coord_cartesian(xlim = c(-3, 3), ylim = c(0, 1.5)) +
#   labs(y = "-log10(p value)", x = "Estimate\n(Alteration v. S1616 Arm)",
#        colour = "Direction", fill = "Direction", alpha = "Direction") +
#   guides(colour = guide_legend(nrow = 3), alpha = guide_legend(nrow = 3),
#          fill = guide_legend(nrow = 3)) +
#   theme_bw() +
#   theme(text = element_text(size = 5), title = element_text(size = 7),
#         axis.ticks = element_line(size = 0.1),
#         panel.border = element_rect(linewidth = 0.1, colour = 'black'),
#         legend.position = 'bottom', legend.key.height = unit(0.25, 'line'),
#         legend.key.spacing = unit(1, 'pt'), legend.margin = margin(0, 0, 0, 0))
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/extfigure2d.pdf", height = 2.5, width = 2.375, unit = 'in', family="ArialMT")
plot_glm <- glm_res %>% filter(term == "gene_alteration", comparison == "combination_v_ipilimumab") %>% 
  filter(! grepl("LOH", altered_gene)) %>%
  mutate(SYMBOL = ifelse(is.na(SYMBOL), gsub("^(\\w+.+) [G|L|D|C].+\\b", "\\1", altered_gene), SYMBOL)) %>% 
  mutate(marker = ifelse(SYMBOL %in% c(ckb_genes, 'NF1'), TRUE, FALSE), 
         colour = ifelse(marker & estimate>0, "A", ifelse(marker & estimate<0, "B", "C"))) %>%
  filter(abs(estimate)<3) %>%
  group_by(estimate, p.value, marker) %>% slice_head(n = 1) %>% arrange(marker, p.value)
ggplot(plot_glm, aes(x = estimate, y = -log10(p.value), colour = colour, fill = colour, alpha = colour)) +
  geom_hline(linetype = 'dashed', colour = 'grey50', size = 0.1, yintercept = -log10(0.05)) +
  geom_point(size = 1, shape = 21, stroke = 0.1) +
  ggrepel::geom_text_repel(data = filter(plot_glm, marker),
                           aes(label = altered_gene), size = 5*0.36, max.overlaps = 10, colour = 'black', show.legend = FALSE,
                           min.segment.length = 0.1) +
  scale_fill_manual(values = c(as.character(arm_colors), 'white'),
                    labels = c('Alteration associated with Combination arm','Alteration associated with Ipilimumab arm','Other alterations')) +
  scale_colour_manual(values = c('black','black','black'),
                      labels = c('Alteration associated with Combination arm','Alteration associated with Ipilimumab arm','Other alterations')) +
  scale_alpha_manual(values = c(1, 1, 0.8),
                     labels = c('Alteration associated with Combination arm','Alteration associated with Ipilimumab arm','Other alterations')) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(0, 1.5)) +
  labs(y = "-log10(p value)", x = "Estimate\n(Alteration v. S1616 Arm)",
       colour = "Direction", fill = "Direction", alpha = "Direction") +
  guides(colour = guide_legend(nrow = 3), alpha = guide_legend(nrow = 3),
         fill = guide_legend(nrow = 3)) +
  theme_bw() +
  theme(text = element_text(size = 5), title = element_text(size = 7),
        axis.ticks = element_line(size = 0.1),
        panel.border = element_rect(linewidth = 0.1, colour = 'black'),
        legend.position = 'bottom', legend.key.height = unit(0.25, 'line'),
        legend.key.spacing = unit(1, 'pt'), legend.margin = margin(0, 0, 0, 0))
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/extfigure2d.pdf", height = 2.5, width = 2.375, unit = 'in', family="ArialMT")


##### Absence of genetic drivers in responding tumors #####
wes_paired <- wes %>% group_by(subject, arm, resp, resp2, timepoint) %>% count %>% spread(timepoint, n, fill = 0) %>%
  filter(baseline>0 & `ontx c2`>0)
wes_paired %>% group_by(arm, resp2) %>% count

var_paired <- var %>% 
  left_join(var_ann) %>% left_join(wes) %>% inner_join(wes_paired) %>% 
  select(subject, wes.id, CHROM, POS, REF, ALT, SYMBOL, Mutation, in_ckb, `Protein Effect`, description, tumorVAF)
driver_paired <- var_paired %>% filter((SYMBOL %in% ckb_genes & !is.na(description)) | (SYMBOL == "NF1" & description == "LOF"))

order_wes <- wes %>% inner_join(wes_paired) %>% 
  select(subject, arm, resp, resp2, representative_wes, wes.id, timepoint, wes.sequenza_purity)
order_wes$timepoint_ <- factor(order_wes$timepoint, levels = c('prior to pd1','baseline','ontx c2'))
order_wes <- order_wes %>% arrange(arm, resp2, subject, timepoint_, desc(wes.sequenza_purity))
order_wes %>% filter(! subject %in% driver_paired$subject) # PT0524 doesn't have any driver mutations
#
driver_paired <- driver_paired %>% filter((tumorVAF>0.05 & in_ckb) | (SYMBOL == "NF1" & description == "LOF")) %>%
  select(-wes.id, -tumorVAF) %>% distinct %>% full_join(order_wes) %>%
  left_join(driver_paired) %>% filter(!is.na(SYMBOL))
# Narrow down so there's up to 2 results shown per gene
check_numbers <- driver_paired %>% group_by(SYMBOL, subject, wes.id) %>% mutate(n_hits = n()) %>% filter(n_hits>1)
# Remove from driver_paired
remove <- check_numbers %>% ungroup %>% 
  select(SYMBOL, subject, CHROM, POS, REF, ALT, SYMBOL, Mutation) %>% distinct
driver_paired <- driver_paired %>% anti_join(remove)
# Only add back top 2
add_back <- check_numbers %>%
  arrange(subject, wes.id, in_ckb, desc(tumorVAF)) %>% group_by(subject, SYMBOL) %>% slice_head(n = 2) %>% 
  select(subject, CHROM, POS, REF, ALT, SYMBOL, Mutation) %>% distinct %>% 
  # mutate(y = c(-0.25, 0.25), y_min = c(-0.5, 0), y_max = c(0, 0.5)) %>%
  mutate(index = 1:n(),
         x = ifelse(1:n() == 1, list(c(-0.5, -0.5, 0.5)), list(c(-0.5, 0.5, 0.5))), 
         y = ifelse(1:n() == 1, list(c(-0.5, 0.5, 0.5)), list(c(-0.5, 0.5, -0.5)))) %>% 
  inner_join(check_numbers)
driver_paired <- driver_paired %>% full_join(add_back)
# Set order of values
driver_paired$wes.id <- factor(driver_paired$wes.id, levels = order_wes$wes.id)
driver_paired$subject <- factor(driver_paired$subject, levels = unique(order_wes$subject))
order_genes <- driver_paired %>% group_by(SYMBOL) %>% summarise(n_pts = length(unique(subject))) %>% arrange(-n_pts)
driver_paired$SYMBOL <- factor(driver_paired$SYMBOL, levels = rev(order_genes$SYMBOL))
driver_paired <- driver_paired %>% arrange(arm, resp2, desc(SYMBOL), wes.id)
driver_paired$subject <- factor(driver_paired$subject, levels = unique(driver_paired$subject))
driver_paired$wes.id <- factor(driver_paired$wes.id, levels = unique(driver_paired$wes.id))
#
order_wes$subject <- factor(order_wes$subject, levels = levels(driver_paired$subject))
order_wes$wes.id <- factor(order_wes$wes.id, levels = levels(driver_paired$wes.id))
wes_x <- order_wes %>% filter(!is.na(subject)) %>% group_by(subject) %>% summarise(x = mean(as.numeric(wes.id)))
#
blank <- filter(order_wes, !is.na(subject)) %>% select(subject, wes.id) %>% cross_join(order_genes)
blank$SYMBOL <- factor(blank$SYMBOL, levels = levels(driver_paired$SYMBOL))
blank$subject <- factor(blank$subject, levels = levels(driver_paired$subject))
blank$wes.id <- factor(blank$wes.id, levels = levels(driver_paired$wes.id))
blank <- blank[complete.cases(blank),]
heatmap <- ggplot(driver_paired, aes(x = as.numeric(wes.id), y = SYMBOL, fill = tumorVAF)) +
  geom_tile(data = blank, fill = 'grey90', colour = 'white') +
  geom_tile(data = filter(driver_paired, is.na(n_hits)), colour = 'white', size = 0.1) +
  geom_polygon(data = filter(driver_paired, !is.na(n_hits)) %>% unnest(c(x, y)),
               aes(x = as.numeric(wes.id)+x, y = as.numeric(SYMBOL)+y, group = interaction(wes.id, Mutation)),
               colour = 'white', size = 0.1) +
  geom_point(data = filter(driver_paired, is.na(n_hits), is.na(tumorVAF)), aes(shape = "Not detected"), size = 1) +
  geom_point(data = filter(driver_paired, !is.na(n_hits), is.na(tumorVAF)), 
             aes(x = ifelse(index == 1, as.numeric(wes.id)-0.25, as.numeric(wes.id)+0.25), 
                 y = ifelse(index == 1, as.numeric(SYMBOL)+0.25, as.numeric(SYMBOL)-0.25),
                 shape = "Not detected"), size = 0.5) +
  scale_x_continuous(breaks = wes_x$x, labels = wes_x$subject) +
  # scale_x_continuous(breaks = 1:length(levels(driver_paired$wes.id)), labels = levels(driver_paired$wes.id)) +
  scale_fill_viridis_c(limits = c(-0.05,1.05), option = 'mako', na.value = NA) +
  scale_shape_manual(values = c(4), na.value = NA) +
  coord_cartesian(expand = 0) + 
  facet_grid(. ~ subject, scales = 'free_x', space = 'free_x') +
  labs(y = "Driver gene", shape = "Detection", fill = "VAF") +
  guides(shape = guide_legend(title.position = 'top'),
         fill = guide_colorbar(title.position = "top", ticks.colour = 'black', ticks.linewidth = 0.1, frame.colour = 'black', frame.linewidth = 0.1)) +
  theme_bw() +
  theme(title = element_text(size = 7), legend.title = element_text(size = 7), legend.text = element_text(size = 5),
        axis.text.y = element_text(size = 5, face = 'italic'), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5),
        axis.ticks.length = unit(2, 'pt'),
        strip.text = element_blank(), axis.title.x = element_blank(),
        legend.position = 'none')# + theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1))

# Metadata
meta <- ggplot(filter(order_wes, !is.na(subject)), aes(x = wes.id)) +
  geom_tile(colour = 'white', size = 0.1, aes(y = 1, fill = arm)) +
  geom_tile(colour = 'white', size = 0.1, aes(y = 2, fill = resp)) +
  geom_tile(colour = 'white', size = 0.1, aes(y = 3, fill = factor(cut(wes.sequenza_purity, c(0, 0.25, 0.5, 0.75, 1))))) +
  geom_tile(colour = 'white', size = 0.1, aes(y = 4, fill = timepoint)) +
  scale_y_reverse(breaks = 1:4, labels = c('arm','bor','purity','timepoint')) +
  scale_fill_manual(values = as.character(c(arm_colors, resp_colors, viridis::viridis(5, option = "mako"), timepoint_colors[c(3,1,2)])), drop = FALSE,
                    breaks = c(names(arm_colors), names(resp_colors), "(0,0.25]","(0.25,0.5]","(0.5,0.75]","(0.5,0.75]","(0.75,1]", names(timepoint_colors[c(3,1,2)]))) +
  facet_grid(. ~ subject, scales = 'free_x', space = 'free_x') +
  coord_cartesian(expand = 0) + 
  theme_bw() +
  theme(text = element_text(size = 5), title = element_text(size = 7),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title = element_blank(), 
        strip.text = element_blank(),
        # strip.background = element_blank(), strip.text = element_text(size = 5, angle = 90, vjust = ), strip.clip = 'off',
        legend.position = 'none')# + theme(legend.position = 'bottom', axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 1))

meta + 
  heatmap + 
  plot_layout(nrow = 2, heights = c(4, length(levels(driver_paired$SYMBOL)))) &
  theme(panel.spacing.x = unit(2, 'pt'), plot.margin = unit(c(1,1,1,1), 'pt'), axis.ticks = element_line(size = 0.1),
        panel.border = element_rect(linewidth = 0.1, colour = 'black'))
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figure1E.pdf", height = 7.4/3, width = 7.2, unit = 'in', family="ArialMT")

# Legends
fill_ <- c('arm','resp',"factor(cut(wes.sequenza_purity, c(0, 0.25, 0.5, 0.75, 1)))",'timepoint')
colors_ <- list(arm_colors, resp_colors, viridis::viridis(5, option = "mako"), timepoint_colors[c(3,1,2)])
breaks_ <- list(names(arm_colors), names(resp_colors), c("(0,0.25]","(0.25,0.5]","(0.5,0.75]","(0.5,0.75]","(0.75,1]"), names(timepoint_colors[c(3,1,2)]))
title_ <- c('S1616 Arm','BOR','Purity (Sequenza)','Timepoint')
legends <- lapply(1:4, function(i){
  plot <- ggplot(filter(order_wes, !is.na(subject)), aes(x = wes.id)) +
    geom_tile(colour = 'white', size = 0.1, aes_string(y = "1", fill = fill_[i])) +
    scale_fill_manual(values = colors_[[i]], breaks = breaks_[[i]]) +
    labs(fill = title_[i]) +
    guides(fill = guide_legend(nrow = 2, override.aes = list(colour = 'black', size = 0.25)))+
    theme(legend.key.width = unit(0.5, 'line'), legend.key.height = unit(0.5, 'line'), 
          legend.justification = c(0,1),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5),
          plot.margin = unit(c(1,1,1,1), 'pt'))
  get_legend(plot) %>% as_ggplot()
})
hm_legend <- heatmap + theme(legend.text = element_text(size = 5),
                             legend.title = element_text(size = 5),
                             legend.position = 'bottom', legend.justification = c(0,1))
hm_legend <- get_legend(hm_legend) %>% as_ggplot()
wrap_plots(legends) + hm_legend + plot_layout(nrow = 1, widths = c(1,1,1.5,1.5,3,1))
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figure1Elegend.pdf", height = 0.7, width = 7.2, unit = 'in', family="ArialMT")


##### Oncoprint ###### 
get_wes <- wes_rep %>% select(wes.id, arm, resp2) %>%
  mutate(group = ifelse(arm == "Ipilimumab", "Ipilimumab", paste0(arm, " ", resp2)))
get_wes$group <- factor(get_wes$group, levels = c("Combination CR/PR","Combination SD/PD","Ipilimumab"))
#
get_cn <- cn_rep %>% filter(description %in% c('CN>4','Deep del'), hgnc_symbol %in% c(ckb_genes, 'B2M'))
get_cn$description <- factor(get_cn$description, levels = c('GOF','LOF','Other mutation','CN>4','Deep del'))
#
get_vars <- var_rep %>% filter(SYMBOL %in% c(ckb_genes, 'NF1'), tumorVAF>0.05)
get_vars <- get_vars %>% arrange(in_ckb, description, -tumorVAF) %>% group_by(wes.id, SYMBOL) %>% slice_head(n=1)
cngenes <- get_cn %>% select(hgnc_symbol) %>% distinct
get_genes <- get_vars %>% group_by(SYMBOL) %>% summarise(n = length(unique(wes.id))) %>% 
  full_join(cngenes, by = c('SYMBOL' = 'hgnc_symbol')) %>%
  filter(n>0.10*nrow(get_wes) | SYMBOL %in% c('JAK1','JAK2','B2M')) %>%
  mutate(goi = case_when(SYMBOL == "NRAS" ~ 1,
                         SYMBOL == "BRAF" ~ 2,
                         SYMBOL == "NF1" ~ 3,
                         SYMBOL == "KIT" ~ 4,
                         SYMBOL %in% c('JAK1','JAK2','B2M') ~ 5,
                         TRUE ~ 6)) %>%
  mutate(goi_next = ifelse(SYMBOL %in% c('BRAF','NRAS','NF1','KIT'), 1, ifelse(SYMBOL %in% c('JAK1','JAK2','B2M'), 2, 3)))
get_vars <- get_vars %>% inner_join(get_vars) %>% mutate(description = ifelse(is.na(description), 'Other mutation', description))
get_vars %>% group_by(description) %>% count()
get_vars$description <- factor(get_vars$description, levels = c('GOF','LOF','Other mutation','CN>4','Deep del'))
get_vars %>% group_by(description) %>% count()
#
all_vars <- get_vars %>% full_join(get_cn, by = c(intersect(colnames(get_vars), colnames(get_cn)), "SYMBOL" = "hgnc_symbol")) %>%
  select(wes.id, SYMBOL, description) %>% mutate(i = as.numeric(description))
var_i <- all_vars %>% arrange(i) %>% group_by(wes.id, SYMBOL) %>% summarise(i = min(i))
var_i <- get_wes %>% cross_join(get_genes) %>% left_join(var_i) %>% mutate(i = ifelse(is.na(i) | i == 3, 10, i))
var_i <- var_i %>% 
  arrange(goi_next, -n)
order_genes <- rev(unique(var_i$SYMBOL))
#
var_i$SYMBOL <- factor(var_i$SYMBOL, levels = order_genes)
# var_i$wes.id <- factor(var_i$wes.id, levels = order_wes)
var_i <- var_i %>% mutate(use_order = ifelse(i<10, 1, 2)) %>% arrange(group, use_order, desc(SYMBOL), i)
order_wes <- unique(var_i$wes.id)
var_i$wes.id <- factor(var_i$wes.id, levels = order_wes)
#
for_vars <- get_vars %>% select(wes.id, SYMBOL, description, in_ckb) %>%
  left_join(get_wes) %>% left_join(get_genes)
for_vars$SYMBOL <- factor(for_vars$SYMBOL, levels = order_genes)
for_vars$wes.id <- factor(for_vars$wes.id, levels = order_wes)
#
for_cn <- get_cn %>% select(wes.id, hgnc_symbol, description) %>%
  left_join(get_wes) %>% left_join(get_genes, by = c("hgnc_symbol" = "SYMBOL")) %>%
  filter(hgnc_symbol %in% order_genes)
for_cn$SYMBOL <- factor(for_cn$hgnc_symbol, levels = order_genes)
for_cn$wes.id <- factor(for_cn$wes.id, levels = order_wes)
# impact_colors <- c(flatui_aussie$`pure apple`, flatui_aussie$`june bud`,
#                    flatui_aussie$`carmine pink`, flatui_aussie$`pink glamour`,
#                    flatui_aussie$`hint of ice`, flatui_aussie$`coastal breeze`,
#                    flatui_aussie$`exodus fruit`,
#                    # flatui_aussie$`deep cove`, #
#                    flatui_aussie$heliotrope,
#                    flatui_aussie$`turbo`)
onco <- ggplot(var_i, aes(x = as.numeric(wes.id), y = as.numeric(SYMBOL),
                  fill = description)) +
  geom_tile(fill = 'grey90', colour = 'white', linewidth = 0.1) +
  geom_tile(data = for_cn, colour = 'white', linewidth = 0.1) +
  #
  geom_tile(data = filter(for_vars, !is.na(in_ckb), !is.na(goi_next)), aes(colour = 'in_ckb'), linewidth = 0.3, width = 0.5, height = 0.5) +
  geom_tile(data = filter(for_vars, is.na(in_ckb), !is.na(goi_next)), colour = NA, width = 0.5, height = 0.5) +
  #
  scale_fill_manual(breaks = levels(get_vars$description),
                    values = as.character(c(flatui_aussie$`pure apple`, flatui_aussie$`carmine pink`,
                                            flatui_aussie$`exodus fruit`,
                                            flatui_aussie$turbo, flatui_aussie$heliotrope)),
                    drop = FALSE) +
  scale_colour_manual(values = c('black'), labels = 'Annotated in\nClinical Knowledgebase', na.translate = FALSE) +
  #
  scale_y_continuous(breaks = 1:length(order_genes), labels = order_genes,
                     sec.axis = sec_axis(~., breaks = 1:length(order_genes), labels = order_genes)) +
  #
  coord_cartesian(expand = 0) +
  facet_grid(goi_next ~ group, scales = 'free', space = 'free') +
  #
  labs(colour = NULL, y = "Gene", fill = "Alteration type") +
  guides(colour = guide_legend(override.aes = list(fill = 'white')),
         fill = guide_legend(override.aes = list(colour = 'white', size = 0.1))) +
  theme_bw() +
  theme(text = element_text(size = 5), title = element_text(size = 7),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        strip.background = element_blank(), strip.text.y = element_blank(), strip.clip = 'off',
        axis.text.y = element_text(size = 5, face = 'italic'),
        panel.spacing.x = unit(2, 'pt'), plot.margin = unit(c(1,1,1,1), 'pt'), axis.ticks = element_line(size = 0.1),
        panel.border = element_rect(linewidth = 0.1, colour = 'black'))

for_wes <- wes_rep %>%
  mutate(group = ifelse(arm == "Ipilimumab", "Ipilimumab", paste0(arm, " ", resp2)))
for_wes$group <- factor(for_wes$group, levels = c("Combination CR/PR","Combination SD/PD","Ipilimumab"))
for_muttype <- wes_ %>%
  mutate(group = ifelse(arm == "Ipilimumab", "Ipilimumab", paste0(arm, " ", resp2)))
for_muttype$muttype <- factor(for_muttype$muttype, levels = rev(c("muttype_C>T","muttype_CC>TT","muttype_T[T>A]T (UV, SBS7c)","muttype_Other")))
for_muttype$group <- factor(for_muttype$group, levels = c("Combination CR/PR","Combination SD/PD","Ipilimumab"))
for_wes$wes.id <- factor(for_wes$wes.id, levels = levels(var_i$wes.id))
for_muttype$wes.id <- factor(for_muttype$wes.id, levels = levels(var_i$wes.id))

plot1 <- ggplot(for_wes, aes(x = as.numeric(wes.id), y = TMB)) +
  geom_tile(aes(y = 5), fill = NA, colour = NA) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 10, colour = 'red', alpha = 0.2) +
  coord_cartesian(expand = 0) +
  scale_y_log10(breaks = c(10, 100)) +
  facet_grid(. ~ group, scales = 'free_x', space = 'free_x') +
  labs(y = "TMB (Mut/Mb)") +
  theme_bw() +
  theme(text = element_text(size = 5), title = element_text(size = 7),
        panel.spacing.x = unit(2, 'pt'), plot.margin = unit(c(1,1,1,1), 'pt'), axis.ticks = element_line(size = 0.1),
        panel.border = element_rect(linewidth = 0.1, colour = 'black'),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        strip.text = element_blank())

plot_subjects <- ggplot(for_wes, aes(x = as.numeric(wes.id))) +
  geom_tile(colour = 'black', linewidth = 0.1, aes(y = 1, fill = factor(as.numeric(arm)))) +
  geom_tile(colour = 'black', linewidth = 0.1, aes(y = 2, fill = factor(as.numeric(resp)+10))) +
  geom_tile(colour = 'black', linewidth = 0.1, aes(y = 3, fill = factor(as.numeric(subtype)+20))) +
  geom_tile(colour = 'black', linewidth = 0.1, aes(y = 4, fill = factor(as.numeric(prior_mapki)+30))) +
  coord_cartesian(expand = 0) +
  scale_y_reverse(breaks = 1:4, labels = c('Arm','BOR','Histology','Prior MAPKi')) +
  scale_fill_manual(values = as.character(c(arm_colors, resp_colors, subtype_colors, prior_colors)), na.value = 'white') +
  labs(y = "Subject") +
  facet_grid(. ~ group, scales = 'free_x', space = 'free_x') +
  theme_bw() +
  theme(text = element_text(size = 5), title = element_text(size = 7),
        panel.spacing.x = unit(2, 'pt'), plot.margin = unit(c(1,1,1,1), 'pt'), axis.ticks = element_line(size = 0.1),
        legend.key.width = unit(0.5, 'line'),
        legend.key.height = unit(0.5, 'line'),
        panel.border = element_rect(linewidth = 0.1, colour = 'black'),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        strip.text = element_blank(), legend.position = 'none')

# Legends
fill_ <- c('arm','resp',"subtype",'prior_mapki')
colors_ <- list(arm_colors, resp_colors, subtype_colors, prior_colors)
breaks_ <- list(names(arm_colors), names(resp_colors), names(subtype_colors), names(prior_colors))
title_ <- c('S1616 Arm','BOR','Histology','Prior MAPKi')
legends <- lapply(1:4, function(i){
  plot <- ggplot(for_wes, aes(x = wes.id)) +
    geom_tile(colour = 'white', size = 0.1, aes_string(y = "1", fill = fill_[i])) +
    scale_fill_manual(values = colors_[[i]], breaks = breaks_[[i]]) +
    labs(fill = title_[i]) +
    guides(fill = guide_legend(nrow = 2, override.aes = list(colour = 'black', size = 0.25)))+
    theme(legend.key.width = unit(0.5, 'line'), legend.key.height = unit(0.5, 'line'), 
          legend.justification = c(0,1),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5),
          plot.margin = unit(c(1,1,1,1), 'pt'))
  get_legend(plot) %>% as_ggplot()
})
hm_legend <- onco +
  guides(fill = guide_legend(nrow = 3)) +
  theme(legend.key.width = unit(0.5, 'line'), legend.key.height = unit(0.5, 'line'), 
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.justification = c(0,1), legend.direction = 'vertical', legend.position = 'bottom')
hm_legend <- get_legend(hm_legend) %>% as_ggplot()
wrap_plots(legends) + hm_legend + plot_layout(nrow = 1, widths = c(1,1,1.5,1.5,3,1))
ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/extfigure2a_legend.pdf", height = 1, width = 7, unit = 'in', family="ArialMT")


plot1 + (onco + theme(legend.position = 'none')) + plot_subjects +
  plot_layout(ncol = 1, heights = c(length(order_genes)/5,length(order_genes),4)) &
  theme(plot.margin = unit(c(2,2,2,2), 'pt'))
ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/extfigure2a.pdf", height = 5, width = 7, unit = 'in', family="ArialMT")


##### END #####
