library(tidyverse)
library(ggpubr)
library(broom)
library(patchwork)
library(ggalluvial)

options(future.globals.maxSize = 24000 * 1024^2)

setwd("~/manuscript/Biopsy analysis of trial S1616/tables")
load("colors.Rda")
load("~/dat/shared/palettes.Rda")

# Plot aesthetics
ae <- list(theme(text = element_text(size = 5), title = element_text(size = 7),
  axis.ticks = element_line(size = 0.1), panel.border = element_rect(linewidth = 0.1, colour = 'black'),
  legend.key.height = unit(0.25, 'line'), legend.key.width = unit(0.25, 'line')
))

# Load in patient and sample annotation
patients <- data.table::fread("table_s1.tsv")
samples <- data.table::fread("table_s2.tsv")

# Set factor levels
patients <- patients %>%
  mutate(group = ifelse(arm == "Ipilimumab", "Ipi", ifelse(resp2 == "CR/PR", "Comb. CR/PR", "Comb. SD/PD")))
patients$group <- factor(patients$group, levels = c("Comb. CR/PR", "Comb. SD/PD", "Ipi"))
patients$arm <- factor(patients$arm, levels = c("Combination", "Ipilimumab"))
patients$resp <- factor(patients$resp, levels = c("CR","PR","SD","PD"))
patients$resp2 <- factor(ifelse(patients$resp %in% c('CR','PR'), 'CR/PR', 'SD/PD'), levels = c('CR/PR','SD/PD'))
samples$timepoint <- factor(samples$timepoint, levels = c('baseline','prior to pd1','ontx c2'))

# Annotation of RNA
rna <- samples %>% filter(!is.na(rnaseq.id)) %>% left_join(patients)
rna %>% group_by(arm, resp2, timepoint) %>% count %>% spread(timepoint, n)
rna %>% group_by(arm, resp2, timepoint, subject) %>% summarise %>% 
  group_by(arm, resp2, timepoint) %>%
  count %>% spread(timepoint, n)

# Load RNAseq results
deg <- readRDS("differentially_expressed_genes.Rds")
table(deg$name)
gsea <- readRDS("gene_set_enrichment_analysis.Rds")
gsva <- readRDS("single_sample_gsea.Rds")
tcr <- readRDS("tcr.Rds")
bcr <- readRDS("bcr.Rds")

# Responders vs. nonresponders at baseline
name_ = "Base, Combo: CR/PR v SD/PD"
labels_ <- rna %>% filter(timepoint == "baseline", arm == "Combination") %>% group_by(arm, group, resp2) %>% count %>%
  rowwise %>% summarise(label = paste0(group, " (n=", n, ")"))
labels_ <- labels_$label
deg_ <- deg %>% filter(name == name_)
gsea_ <- gsea %>% filter(name == name_)
deg_ %>% group_by(name, sig) %>% count
gsea_ %>% group_by(name, sig) %>% count
deg_ %>% group_by(name, sig, dir_i) %>% count
gsea_ %>% group_by(name, sig, dir_i) %>% count

## Volcano plots: # Downsample the non significant results
plot_deg <- deg_ %>% filter(is.na(sig)) %>% slice_sample(prop = 0.10, replace = FALSE) %>%
  full_join(filter(deg_, sig))
plot_deg_ <- ggplot(plot_deg, aes(x = log2FoldChange, y = -log10(pvalue), fill = interaction(sig, dir_i), colour = interaction(sig, dir_i))) +
  geom_point(shape = 21, stroke = 0.1) + 
  ggrepel::geom_text_repel(aes(label = ifelse(sig, Gene, NA)), min.segment.length = 0.1, force = 5,
                           size = 5*0.36, max.overlaps = 30, colour = 'black', show.legend = FALSE) +
  geom_text(data = tibble(n = c(107, 531), x = c(-8, 8), hjust = c(0,1)), aes(x = x, y = 0, label = n, hjust = hjust), 
            inherit.aes = FALSE, size = 5*0.36) +
  scale_fill_manual(values = rev(as.character(resp2_colors)), na.value = 'grey95',
                    labels = c(paste0("FDR<0.05, ", rev(labels_)), "n.s.")) +
  scale_colour_manual(values = c('black','black'), na.value = 'grey50',
                      labels = c(paste0("FDR<0.05, ", rev(labels_)), "n.s.")) +
  coord_cartesian(xlim = c(-8, 8)) +
  labs(fill = "Significance", colour = "Significance", x = "Log2 Fold Change", y = "-log10 (p value)") +
  theme_bw() + ae
plot_gsea <- gsea_ %>% filter(is.na(sig)) %>% 
  full_join(filter(gsea_, sig)) %>% mutate(lab = ifelse(is.na(sig), NA, gsub("(HALLMARK_)|(KEGG_)", "", pathway)))
plot_gsea_ <- ggplot(plot_gsea, aes(x = NES, y = -log10(pval), fill = interaction(sig, dir_i), colour = interaction(sig, dir_i))) +
  geom_point(shape = 21, stroke = 0.1) + 
  ggrepel::geom_text_repel(aes(label = lab), min.segment.length = 0.1, force = 5,
                           size = 5*0.36, max.overlaps = 30, colour = 'black', show.legend = FALSE, direction = 'y', hjust = 1) +
  geom_text(data = tibble(n = c(1, 170), x = c(-4, 4), hjust = c(0,1)), aes(x = x, y = 0, label = n, hjust = hjust), 
            inherit.aes = FALSE, size = 5*0.36) +
  scale_fill_manual(values = rev(as.character(resp2_colors)), na.value = 'grey95',
                    labels = c(paste0("FDR<0.05, ", rev(labels_)), "n.s.")) +
  scale_colour_manual(values = c('black','black'), na.value = 'grey50',
                      labels = c(paste0("FDR<0.05, ", rev(labels_)), "n.s.")) +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(fill = "Significance", colour = "Significance", x = "NES", y = "-log10 (p value)") +
  theme_bw() +
  ae
plot_deg_ + plot_gsea_ +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom')
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figure3AB.pdf", height = 3, width = 5, unit = 'in', family="ArialMT")

select_gsea <- gsea_ %>% filter(sig) %>% arrange(-NES) %>% slice_head(n = 20) %>%
  select(pathway)
select_geasm <- select_gsea %>% inner_join(gsva, by = c("pathway" = "gene_set")) %>%
  filter(timepoint == "baseline", arm == "Combination") %>%
  pivot_wider(id_cols = 'rnaseq.id', names_from = 'pathway', values_from = gsva) %>%
  column_to_rownames("rnaseq.id") %>% data.matrix
M = round(cor(select_geasm, method = 'spearman'), 2)
P = ggcorrplot::cor_pmat(select_geasm, method = "spearman")
cp <- ggcorrplot::ggcorrplot(M,
                 hc.order = TRUE,
                 type = "upper", p.mat = P, show.diag = TRUE
)
cp <- cp$data %>% arrange(-value) %>%
  mutate(label = ifelse(pvalue<0.001, "p<0.001", ifelse(pvalue<0.01, "p<0.01", ifelse(pvalue<0.05, "p<0.05", NA))))
cp$label <- factor(cp$label, levels = c("p<0.001", "p<0.01", "p<0.05"))
ggplot(cp, aes(x = Var1, y = Var2, fill = value, shape = label)) +
  geom_tile(colour = 'black', stroke = 0.1) + 
  geom_point(size = 0.25) +
  coord_flip() +
  scale_fill_distiller(palette = 'RdBu', limits = c(-1.1, 1.1), breaks = c(-1, 0, 1)) +
  scale_shape_manual(values = c(3, 1, 16), 
                     drop = TRUE, na.translate = FALSE) +
  labs(shape = "Significance", fill = "Spearman rho") +
  guides(fill = guide_colorbar(frame.linewidth = 0.1, frame.colour = 'black',
                               ticks.colour = 'black', ticks.linewidth = 0.1)) +
  theme_minimal() +
  theme(text = element_text(size = 5), title = element_text(size = 7),
        axis.ticks = element_line(size = 0.1),
  ) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 3, hjust = 1, vjust = 0.5, angle = 90),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title = element_blank(),
        legend.justification = 'left', legend.key.height = unit(0.25, 'line'), legend.key.width = unit(0.25, 'line'))
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figure3C.pdf", height = 3, width = 2.2, unit = 'in', family="ArialMT")

# "Combo, CR/PR: Ontx v Base"
name_ = "Combo, CR/PR: Ontx v Base"
labels_ <- rna %>% 
  filter(resp2 == "CR/PR", arm == "Combination", timepoint %in% c('baseline','ontx c2')) %>% 
  group_by(arm, group, resp2, timepoint) %>% count %>%
  rowwise %>% summarise(label = paste0(group, ", ", timepoint," (n=", n, ")"))
rna %>% 
  filter(resp2 == "CR/PR", arm == "Combination", timepoint %in% c('baseline','ontx c2')) %>%
  group_by(subject, timepoint) %>% count %>% spread(timepoint, n) %>%
  group_by(baseline, `ontx c2`) %>% count
labels_ <- labels_$label
deg_ <- deg %>% filter(name == name_)
deg_base <- deg %>% filter(name == "Base, Combo: CR/PR v SD/PD")
gsea_ <- gsea %>% filter(name == name_)
gsea_base <- gsea %>% filter(name == "Base, Combo: CR/PR v SD/PD")
deg_ %>% group_by(name, sig) %>% count
gsea_ %>% group_by(name, sig) %>% count
deg_ %>% group_by(name, sig, dir_i) %>% count
gsea_ %>% group_by(name, sig, dir_i) %>% count

## Compare to things enriched at baseline
deg_ %>% left_join(deg_base, by = "Gene", suffix = c("", "_base")) %>%
  filter(sig) %>%
  group_by(name, sig, dir_i, sig_base, dir_i_base) %>% count
gsea_ %>% left_join(gsea_base, by = "pathway", suffix = c("", "_base")) %>%
  # filter(sig) %>%
  group_by(name, sig, dir_i, sig_base, dir_i_base) %>% count
comp_gsea <- gsea_ %>% full_join(gsea_base, by = "pathway", suffix = c("", "_base")) %>%
  mutate(sig = ifelse(is.na(sig), '2nottime', '1time'), sig_base = ifelse(is.na(sig_base), '2notbase', "1base"))
ggplot(comp_gsea, aes(x = NES_base, y = NES, 
                      fill = interaction(sig_base, sig))) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.1) +
  geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.1) +
  geom_point(shape = 21, stroke = 0.1, size = 4) +
  geom_text(data = filter(comp_gsea, NES>0 & NES_base>0 & sig_base != "2notbase"), 
                            aes(label = gsub("HALLMARK_|KEGG_", "", pathway)), hjust = 0,
                            size = 8*0.36, max.overlaps = 5, show.legend = FALSE) +
  coord_cartesian(xlim = c(0,4), ylim = c(0,2.5)) +
  # ggrepel::geom_label_repel(data = filter(comp_gsea, NES>0 & NES_base>0 & sig_base != "2notbase"), 
  #                          aes(label = gsub("HALLMARK_|KEGG_", "", pathway)), hjust = 0,
  #                          seed = 123, alpha = 0.5, colour = 'white',  min.segment.length = 0.1,
  #                          size = 5*0.36, max.overlaps = 5, show.legend = FALSE) +
  # ggrepel::geom_label_repel(data = filter(comp_gsea, NES>0 & NES_base>0 & sig_base != "2notbase"), 
  #                           aes(label = gsub("HALLMARK_|KEGG_", "", pathway)), hjust = 0,
  #                           seed = 123, fill = NA, min.segment.length = 0.1,
  #                           size = 5*0.36, max.overlaps = 5, show.legend = FALSE) +
  scale_fill_manual(values = c(as.character(c(flatui_aussie$`spiced nectarine`, timepoint_colors$`ontx c2`, resp2_colors$`CR/PR`)), 'white'),
                     labels = c('Significant in both',
                                'FDR<0.05, Comb. CR/PR, ontx vs baseline',
                                'FDR<0.05, Comb. CR/PR vs. SD/PD, baseline',
                                'n.s. in both')) +
  # coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4)) +
  labs(x = "NES, Comb. CR/PR vs. SD/PD, baseline",
       y = "NES, Comb. CR/PR, ontx vs baseline",
       fill = "Significance") +
  guides(fill = guide_legend(nrow = 2)) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  ae
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figure3D.pdf", height = 3, width = 2.5, unit = 'in', family="ArialMT")

deg_ %>% filter(sig, dir_i == "+") %>% arrange(-log2FoldChange)
net <- gsea_ %>% filter(NES>0, sig) %>% unnest(leadingEdge) %>%
  left_join(deg_, by = c("name", "dir_i", "leadingEdge" = "Gene"), suffix = c("", "_deg")) %>%
  filter(sig_deg, dir_i=="+") %>% arrange(-NES, -log2FoldChange) %>% mutate(n = 1)
net$pathway <- factor(net$pathway, levels = unique(net$pathway))
net$leadingEdge <- factor(net$leadingEdge, levels = unique(net$leadingEdge))
ggplot(net, aes(axis2 = leadingEdge, axis1 = pathway, y = n)) +
  geom_alluvium(aes(fill = pathway), width = 1/8, alpha = 0.6) +
  geom_stratum(width = 1/8, fill = 'white', alpha = 0.9, color = "black", linewidth = 0.25) +
  geom_text(stat = "stratum", size = 5*0.36,
            aes(label = after_stat(stratum), 
                hjust = ifelse(after_stat(stratum) %in% net$pathway, 0, 0), 
                angle = ifelse(after_stat(stratum) %in% net$pathway, -30, 30))) +
  scale_fill_manual(values = as.character(flatui_aussie)[c(1:5,7:9,11)]) +
  coord_flip() +
  theme_void() +
  theme(legend.position = 'none')
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figure3E.pdf", height = 0.75, width = 2.5, unit = 'in', family="ArialMT")

y_ <- c("n", "h", "clonality")
label_ <- c("N clones", "Shannon Diversity (H)", "Clonality")
lapply(1:3, function(i){
  ggplot(tcr, aes(x = timepoint, y = ))
})

tcr_ <- tcr %>% filter(group == "Comb. CR/PR", timepoint %in% c('baseline','ontx c2'), Gene %in% c('TRA','TRB'))
#   stat_summary(aes(y = Score, group = Emotion), fun.y = mean, geom="line", size = 2.2, alpha = 1.2, width = 0.25, colour = 'gray48') +
fig3f <- ggplot(tcr_, aes(x = timepoint, y = h, fill = timepoint)) +
  geom_boxplot(outlier.colour = NA, linewidth = 0.25) +
  stat_summary(aes(group = subject), fun.y = median, geom="line", size = 0.1, colour = 'black', alpha = 0.8) +
  ggbeeswarm::geom_beeswarm(shape = 21, stroke = 0.1) +
  stat_compare_means(comparisons = list(c('baseline','ontx c2')), size = 5*0.36) +
  scale_fill_manual(values = as.character(timepoint_colors)) +
  scale_x_discrete(labels = c('baseline','ontx')) +
  facet_wrap(~ Gene) +
  labs(y = "Shannon Diversity (H)", x = "Timepoint") +
  theme_bw() +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 7, face = 'italic'),
        legend.position = 'none') + ae
fig3g <- ggplot(tcr_, aes(x = timepoint, y = clonality, fill = timepoint)) +
  geom_boxplot(outlier.colour = NA, linewidth = 0.25) +
  stat_summary(aes(group = subject), fun.y = median, geom="line", size = 0.1, colour = 'black', alpha = 0.8) +
  ggbeeswarm::geom_beeswarm(shape = 21, stroke = 0.1) +
  stat_compare_means(comparisons = list(c('baseline','ontx c2')), size = 5*0.36) +
  scale_fill_manual(values = as.character(timepoint_colors)) +
  scale_x_discrete(labels = c('baseline','ontx')) +
  facet_wrap(~ Gene) +
  labs(y = "Clonality", x = "Timepoint") +
  theme_bw() +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 7, face = 'italic'),
        legend.position = 'none') + ae
fig3f + fig3g &
  theme(plot.margin = unit(c(2,2,2,2), 'pt'))
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figure3FG.pdf", height = 1.75, width = 4.75, unit = 'in', family="ArialMT")


# "Combo, SD/PD: Ontx v Base"
name_ = "Combo, SD/PD: Ontx v Base"
labels_ <- rna %>% 
  filter(resp2 == "SD/PD", arm == "Combination", timepoint %in% c('baseline','ontx c2')) %>% 
  group_by(arm, group, resp2, timepoint) %>% count %>%
  rowwise %>% summarise(label = paste0(group, ", ", timepoint," (n=", n, ")"))
rna %>% 
  filter(resp2 == "SD/PD", arm == "Combination", timepoint %in% c('baseline','ontx c2')) %>%
  group_by(subject, timepoint) %>% count %>% spread(timepoint, n) %>%
  group_by(baseline, `ontx c2`) %>% count
labels_ <- labels_$label
deg_ <- deg %>% filter(name == name_)
deg_base <- deg %>% filter(name == "Combo, SD/PD: Ontx v Base")
gsea_ <- gsea %>% filter(name == name_)
gsea_base <- gsea %>% filter(name == "Combo, SD/PD: Ontx v Base")
deg_ %>% group_by(name, sig) %>% count
gsea_ %>% group_by(name, sig) %>% count
deg_ %>% group_by(name, sig, dir_i) %>% count
gsea_ %>% group_by(name, sig, dir_i) %>% count

## Volcano plots: # Downsample the non significant results
base_genes <- filter(deg, sig, dir_i == "+", name == "Base, Combo: CR/PR v SD/PD")$Gene
pick_genes <- c('CCL21','IGHG3','IGHG4','IGHG2','CXCL13','IL2RA','IL21R')
plot_deg <- deg_ %>% filter(is.na(sig)) %>% slice_sample(prop = 0.10, replace = FALSE) %>%
  full_join(filter(deg_, sig)) %>%
  mutate(base = ifelse(is.na(sig) | ! Gene %in% base_genes, "2", "1")) %>%
  mutate(label = ifelse(Gene %in% pick_genes, "1", "2")) %>%
  arrange(desc(label), desc(sig), desc(base))
plot_deg_ <- ggplot(plot_deg, aes(x = log2FoldChange, y = -log10(pvalue),
                     fill = interaction(sig, dir_i), colour = base)) +
  geom_point(shape = 21, stroke = 0.25) + 
  ggrepel::geom_label_repel(aes(label = ifelse(sig & label == "1", Gene, NA)), min.segment.length = 0.1, force = 5,
                           size = 5*0.36, max.overlaps = 30, colour = 'black', show.legend = FALSE, fill = 'white') +
  geom_text(data = tibble(n = c(29, 146), x = c(-8, 8), hjust = c(0,1)), aes(x = x, y = 0, label = n, hjust = hjust), 
            inherit.aes = FALSE, size = 5*0.36) +
  scale_fill_manual(values = as.character(timepoint_colors), na.value = 'grey95',
                    labels = c(paste0("FDR<0.05, ", rev(labels_)), "n.s.")) +
  scale_colour_manual(values = c('black','grey80'), na.value = 'grey80',
                      labels = c("Sig. higher in Comb. CR/PR, base (Fig 2)","")) +
  guides(colour = guide_legend(override.aes = list(size = 3, stroke = 0.5)),
         fill = guide_legend(override.aes = list(stroke = 0, size = 3))) +
  coord_cartesian(xlim = c(-8, 8)) +
  labs(fill = "Significance", colour = "Baseline analysis", x = "Log2 Fold Change", y = "-log10 (p value)") +
  theme_bw() +
  theme(panel.border = element_rect(linewidth = 0.1, colour = 'black'),
        plot.margin = unit(c(1,1,1,1), 'pt'),
        legend.key.size = unit(c(0.5,0.5), 'line'))
base_path <- filter(gsea, sig, dir_i == "+", name == "Base, Combo: CR/PR v SD/PD")$pathway[1:15]
plot_gsea <- gsea_ %>% filter(is.na(sig)) %>% 
  full_join(filter(gsea_, sig)) %>% 
  mutate(base = ifelse(is.na(sig) | ! pathway %in% base_path, "2", "1")) %>%
  mutate(lab = ifelse(is.na(sig) | base == "2", NA, gsub("(HALLMARK_)|(KEGG_)", "", pathway))) %>%
  arrange(desc(sig), desc(base))
plot_gsea_ <- ggplot(plot_gsea, aes(x = NES, y = -log10(pval), fill = interaction(sig, dir_i), colour = base)) +
  geom_point(shape = 21, stroke = 0.1) + 
  ggrepel::geom_text_repel(aes(label = lab), min.segment.length = 0.1, force = 5,
                           size = 5*0.36, max.overlaps = 30, colour = 'black', show.legend = FALSE, direction = 'y', hjust = 1) +
  geom_text(data = tibble(n = c(1, 87), x = c(-3, 3), hjust = c(0,1)), aes(x = x, y = 0, label = n, hjust = hjust), 
            inherit.aes = FALSE, size = 5*0.36) +
  scale_fill_manual(values = as.character(timepoint_colors), na.value = 'grey95',
                    labels = c(paste0("FDR<0.05, ", rev(labels_)), "n.s.")) +
  scale_colour_manual(values = c('black','grey80'), na.value = 'grey80',
                      labels = c("Sig. higher in Comb. CR/PR, base (Fig 2)","")) +
  guides(colour = guide_legend(override.aes = list(size = 3, stroke = 0.5)),
         fill = guide_legend(override.aes = list(stroke = 0, size = 3))) +
  coord_cartesian(xlim = c(-3, 3)) +
  labs(fill = "Significance", colour = "Baseline analysis", x = "NES", y = "-log10 (p value)") +
  theme_bw() +
  theme(panel.border = element_rect(linewidth = 0.1, colour = 'black'),
        plot.margin = unit(c(1,1,1,10), 'pt'),
        legend.key.size = unit(c(0.5,0.5), 'line'))
plot_deg_ + plot_gsea_ +
  plot_layout(guides = 'collect') &
  theme(text = element_text(size = 5), title = element_text(size = 7),
        legend.position = 'right', legend.direction = 'vertical',
        legend.margin = margin(2,2,2,2,'pt'), 
        legend.title.position = 'top', legend.justification = c(0,1))
ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figure4ab.pdf", height = 2.2, width = 7.2, unit = 'in', family="ArialMT")


base_paths <- filter(gsea, sig, dir_i == "+", name == "Base, Combo: CR/PR v SD/PD")$pathway
sdpd_paths <- filter(gsea, sig, dir_i == "+", name == "Combo, SD/PD: Ontx v Base")$pathway
rel_gsva <- gsva %>% filter(gene_set %in% intersect(base_paths, sdpd_paths)) %>%
  filter(grepl('Comb', group), timepoint %in% c('baseline','ontx c2')) %>%
  group_by(gene_set) %>% mutate(gsva = scale(gsva)[,1])
rel_gsva <- rel_gsva %>% group_by(rnaseq.id, subject, timepoint, group) %>% 
  summarise(gsva = mean(gsva))

ggplot(rel_gsva, aes(x = timepoint, y = gsva, fill = group, alpha = timepoint)) +
  geom_boxplot(linewidth = 0.25) +
  stat_summary(geom = 'line', fun = median, aes(group = subject), linewidth = 0.1, alpha =1) +
  ggbeeswarm::geom_beeswarm(shape = 21, size = 1, alpha = 1, stroke = 0.1) +
  scale_fill_manual(values = as.character(resp2_colors)) +
  scale_alpha_manual(values = c(1,0.5))+
  facet_wrap(~ group) +
  labs(y = "Relative expression (ssGSEA)") +
  theme_bw() +
  theme(text = element_text(size = 5), title = element_text(size = 7),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(linewidth = 0.1, colour = 'black'),
        plot.margin = margin(1,1,1,1,'pt'),
        strip.background = element_blank())
ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figure4c.pdf", height = 2.2, width = 1.5, unit = 'in', family="ArialMT")



plot_deg %>% group_by(sig, dir_i, base) %>% count
base_path <- filter(gsea, sig, dir_i == "+", name == "Base, Combo: CR/PR v SD/PD")$pathway
plot_gsea <- gsea_ %>% filter(is.na(sig)) %>% 
  full_join(filter(gsea_, sig)) %>% 
  mutate(base = ifelse(is.na(sig) | ! pathway %in% base_path, "2", "1")) %>%
  mutate(lab = ifelse(is.na(sig) | base == "2", NA, gsub("(HALLMARK_)|(KEGG_)", "", pathway))) %>%
  arrange(desc(sig), desc(base))
plot_gsea %>% group_by(sig, dir_i, base) %>% count
deg_upnr <- plot_deg %>% filter(sig, dir_i == "+", base == "2")
gsea_unr <- plot_gsea %>% filter(sig, dir_i == "+", base == "2")
upnr <- gsea_unr %>% unnest(leadingEdge) %>% inner_join(deg_upnr, by = c("leadingEdge" = "Gene"))

# panel.border = element_rect(linewidth = 0.1, colour = 'black'),
# plot.margin = unit(c(1,1,1,1), 'pt'),
# legend.key.size = unit(c(0.5,0.5), 'line')
tcrbcr <- tcr %>% filter(group == "Comb. SD/PD", timepoint %in% c('baseline','ontx c2'), Gene %in% c('TRA','TRB','IGH','IGL'))
tcrbcr <- bcr %>% filter(group == "Comb. SD/PD", timepoint %in% c('baseline','ontx c2'), Gene %in% c('TRA','TRB','IGH','IGL')) %>%
  full_join(tcrbcr)
tcrbcr$Gene <- factor(tcrbcr$Gene, levels = c('TRA','TRB','IGL','IGH'))
#   stat_summary(aes(y = Score, group = Emotion), fun.y = mean, geom="line", size = 2.2, alpha = 1.2, width = 0.25, colour = 'gray48') +
fig4d <- ggplot(tcrbcr, aes(x = timepoint, y = h, fill = timepoint)) +
  geom_boxplot(outlier.colour = NA, linewidth = 0.25) +
  stat_summary(aes(group = subject), fun.y = median, geom="line", size = 0.1, colour = 'black', alpha = 0.8) +
  ggbeeswarm::geom_beeswarm(shape = 21, stroke = 0.1) +
  stat_compare_means(comparisons = list(c('baseline','ontx c2')), size = 5*0.36) +
  scale_fill_manual(values = as.character(timepoint_colors)) +
  scale_x_discrete(labels = c('baseline','ontx')) +
  facet_wrap(~ Gene, nrow = 1) +
  labs(y = "Shannon Diversity (H)", x = "Timepoint") +
  theme_bw() +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 7, face = 'italic'), axis.text.x = element_blank(),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = 'none') + ae
fig4e <- ggplot(tcrbcr, aes(x = timepoint, y = clonality, fill = timepoint)) +
  geom_boxplot(outlier.colour = NA, linewidth = 0.25) +
  stat_summary(aes(group = subject), fun.y = median, geom="line", size = 0.1, colour = 'black', alpha = 0.8) +
  ggbeeswarm::geom_beeswarm(shape = 21, stroke = 0.1) +
  stat_compare_means(comparisons = list(c('baseline','ontx c2')), size = 5*0.36) +
  scale_fill_manual(values = as.character(timepoint_colors)) +
  scale_x_discrete(labels = c('baseline','ontx')) +
  facet_wrap(~ Gene, nrow = 1) +
  labs(y = "Clonality", x = "Timepoint") +
  theme_bw() +
  theme(strip.background = element_blank(), 
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        legend.position = 'none') + ae
fig4d + fig4e +
  plot_layout(ncol = 1) &
  theme(plot.margin = unit(c(2,2,2,2), 'pt'))
ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figure4de.pdf", height = 3, width = 3.25, unit = 'in', family="ArialMT")

# "Ipi, SD/PD: Ontx v Base"
name_ = "Ipi, SD/PD: Ontx v Base"
labels_ <- rna %>% 
  filter(resp2 == "SD/PD", arm == "Ipilimumab", timepoint %in% c('baseline','ontx c2')) %>% 
  group_by(arm, group, resp2, timepoint) %>% count %>%
  rowwise %>% summarise(label = paste0(group, ", ", timepoint," (n=", n, ")"))
rna %>% 
  filter(resp2 == "SD/PD", arm == "Ipilimumab", timepoint %in% c('baseline','ontx c2')) %>% 
  group_by(subject, timepoint) %>% count %>% spread(timepoint, n) %>%
  group_by(baseline, `ontx c2`) %>% count
labels_ <- labels_$label
deg_ <- deg %>% filter(name == name_)
gsea_ <- gsea %>% filter(name == name_)
deg_ %>% group_by(name, sig) %>% count
gsea_ %>% group_by(name, sig) %>% count
deg_ %>% group_by(name, sig, dir_i) %>% count
gsea_ %>% group_by(name, sig, dir_i) %>% count

base_genes <- filter(deg, sig, dir_i == "+", name == "Base, Combo: CR/PR v SD/PD")$Gene
plot_deg <- deg_ %>% filter(is.na(sig)) %>% slice_sample(prop = 0.10, replace = FALSE) %>%
  full_join(filter(deg_, sig)) %>%
  mutate(base = ifelse(is.na(sig) | ! Gene %in% base_genes, "2", "1")) %>%
  arrange(desc(sig), desc(base))
plot_deg_ <- ggplot(plot_deg, aes(x = log2FoldChange, y = -log10(pvalue),
                                  fill = interaction(sig, dir_i), colour = base)) +
  geom_point(shape = 21, stroke = 0.25) + 
  geom_text(data = tibble(n = c(41, 282), x = c(-8, 8), hjust = c(0,1)), aes(x = x, y = 0, label = n, hjust = hjust), 
            inherit.aes = FALSE, size = 5*0.36) +
  scale_fill_manual(values = as.character(timepoint_colors), na.value = 'grey95',
                    labels = c(paste0("FDR<0.05, ", rev(labels_)), "n.s.")) +
  scale_colour_manual(values = c('black','grey80'), na.value = 'grey80',
                      labels = c("Sig. higher in Comb. CR/PR, base (Fig 2)","")) +
  guides(colour = guide_legend(override.aes = list(size = 3, stroke = 0.5)),
         fill = guide_legend(override.aes = list(stroke = 0, size = 3))) +
  coord_cartesian(xlim = c(-8, 8)) +
  labs(fill = "Significance", colour = "Baseline analysis", x = "Log2 Fold Change", y = "-log10 (p value)") +
  theme_bw() +
  theme(panel.border = element_rect(linewidth = 0.1, colour = 'black'),
        plot.margin = unit(c(1,1,1,1), 'pt'),
        legend.key.size = unit(c(0.5,0.5), 'line'))
base_path <- filter(gsea, sig, dir_i == "+", name == "Base, Combo: CR/PR v SD/PD")$pathway[1:15]
plot_gsea <- gsea_ %>% filter(is.na(sig)) %>% 
  full_join(filter(gsea_, sig)) %>% 
  mutate(base = ifelse(is.na(sig) | ! pathway %in% base_path, "2", "1")) %>%
  mutate(lab = ifelse(is.na(sig) | base == "2", NA, gsub("(HALLMARK_)|(KEGG_)", "", pathway))) %>%
  arrange(desc(sig), desc(base))
plot_gsea_ <- ggplot(plot_gsea, aes(x = NES, y = -log10(pval), fill = interaction(sig, dir_i), colour = base)) +
  geom_point(shape = 21, stroke = 0.1) + 
  ggrepel::geom_text_repel(aes(label = lab), min.segment.length = 0.1, force = 5,
                           size = 5*0.36, max.overlaps = 30, colour = 'black', show.legend = FALSE, direction = 'y', hjust = 1) +
  geom_text(data = tibble(n = c(57, 50), x = c(-3, 3), hjust = c(0,1)), aes(x = x, y = 0, label = n, hjust = hjust), 
            inherit.aes = FALSE, size = 5*0.36) +
  scale_fill_manual(values = as.character(timepoint_colors), na.value = 'grey95',
                    labels = c(paste0("FDR<0.05, ", rev(labels_)), "n.s.")) +
  scale_colour_manual(values = c('black','grey80'), na.value = 'grey80',
                      labels = c("Sig. higher in Comb. CR/PR, base (Fig 2)","")) +
  guides(colour = guide_legend(override.aes = list(size = 3, stroke = 0.5)),
         fill = guide_legend(override.aes = list(stroke = 0, size = 3))) +
  coord_cartesian(xlim = c(-3, 3)) +
  labs(fill = "Significance", colour = "Baseline analysis", x = "NES", y = "-log10 (p value)") +
  theme_bw() +
  theme(panel.border = element_rect(linewidth = 0.1, colour = 'black'),
        plot.margin = unit(c(1,1,1,10), 'pt'),
        legend.key.size = unit(c(0.5,0.5), 'line'))
plot_deg_ + plot_gsea_ +
  plot_layout(guides = 'collect') &
  theme(text = element_text(size = 5), title = element_text(size = 7),
        legend.position = 'right', legend.direction = 'vertical',
        legend.margin = margin(2,2,2,2,'pt'), 
        legend.title.position = 'top', legend.justification = c(0,1))
ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figure5ab.pdf", height = 2.2, width = 7.2, unit = 'in', family="ArialMT")


# filter(gsea, sig, dir_i == "+", name == "Base, Combo: CR/PR v SD/PD")
base_paths <- filter(gsea, sig, dir_i == "+", name == "Base, Combo: CR/PR v SD/PD")$pathway
ipi_paths <- filter(gsea, sig, dir_i == "+", name == "Ipi, SD/PD: Ontx v Base")$pathway
rel_gsva <- gsva %>% filter(gene_set %in% base_paths) %>% #intersect(base_paths, ipi_paths)) %>%
  filter(grepl('Ipi', group) | group == "Comb. CR/PR", timepoint %in% c('baseline','ontx c2')) %>%
  group_by(gene_set) %>% mutate(gsva = scale(gsva)[,1])
rel_gsva <- rel_gsva %>% group_by(rnaseq.id, subject, timepoint, group) %>% 
  summarise(gsva = mean(gsva))

group_colors <- as.character(c(flatui_aussie$`pure apple`, flatui_aussie$`carmine pink`, flatui_aussie$`deep koamaru`))
ggplot(rel_gsva, aes(x = timepoint, y = gsva, fill = group, alpha = timepoint)) +
  geom_boxplot(linewidth = 0.25, outlier.colour = NA) +
  stat_summary(geom = 'line', fun = median, aes(group = subject), linewidth = 0.1, alpha =1) +
  ggbeeswarm::geom_beeswarm(shape = 21, size = 1, alpha = 1, stroke = 0.1) +
  stat_summary(data = filter(rel_gsva, subject == "PT0550"), geom = 'line', fun = median, aes(group = subject, colour = "PT0550"), linewidth = 0.1, alpha =1) +
  geom_point(data = filter(rel_gsva, subject == "PT0550"), aes(colour = "PT0550"), shape = 23, size = 1.25, alpha = 1, stroke = 0.1, fill = 'red3') +
  geom_point(data = filter(rel_gsva, subject == "PT0765"), aes(colour = "PT0765"), shape = 23, size = 1.25, alpha = 1, stroke = 0.1, fill = 'orange2') +
  scale_fill_manual(values = group_colors[c(1,3)]) +
  scale_alpha_manual(values = c(1, 0.2)) +
  scale_colour_manual(values = c("red3","orange2")) +
  facet_wrap(~ group) +
  labs(y = "Relative expression (ssGSEA)", colour = "") +
  guides(fill = 'none', alpha = 'none', colour = guide_legend(nrow = 2, override.aes = list(shape = 23, size = 1.25))) +
  theme_bw() +
  theme(text = element_text(size = 5), title = element_text(size = 7),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = 'bottom', legend.justification = c(0,1), legend.margin = margin(0,0,0,0,'pt'),
        legend.key.size = unit(c(0.5,0.5),'line'), legend.spacing = unit(c(0,0,0,0),'pt'),
        legend.box.margin = margin(0,0,0,0,'pt'), 
        panel.border = element_rect(linewidth = 0.1, colour = 'black'),
        plot.margin = margin(1,1,1,1,'pt'),
        strip.background = element_blank())
ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figure5c.pdf", height = 2.2, width = 1.25, unit = 'in', family="ArialMT")


# Merge TCR/BCR summary metrics with Table S2
# table_s2 <- data.table::fread("table_s2.tsv")


#### END ####
