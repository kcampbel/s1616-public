library(data.table)
library(tidyverse)
library(SPOT)
library(ggpubr)
library(scales)
library(patchwork)
library(ggdendro)
library(spatstat)
library(sf)

options(future.globals.maxSize = 24000 * 1024^2)

setwd("~/manuscript/Biopsy analysis of trial S1616/tables")
load("colors.Rda")
load("~/dat/shared/palettes.Rda")

color_sorted <- c("carmine pink", "greenland green", "turbo", "pure apple", "steel pink",
                  "quince jelly", "exodus fruit", "deep koamaru", "soaring eagle", "hint of ice",
                  "heliotrope","june bud","wizard hat", "beekeeper")
color_sorted <- as.character(flatui_aussie[color_sorted])
group_colors <- as.character(c(flatui_aussie$`pure apple`, flatui_aussie$`carmine pink`, flatui_aussie$`deep koamaru`))
ae <- list(theme(text = element_text(size = 5), title = element_text(size = 7),
           axis.ticks = element_line(size = 0.1), panel.border = element_rect(linewidth = 0.1, colour = 'black'),
           legend.key.height = unit(0.25, 'line'), legend.key.width = unit(0.25, 'line')
))

# Load in patient and sample annotation
patients <- data.table::fread("table_s1.tsv")
samples <- data.table::fread("table_s2.tsv")
samples %>% filter(!is.na(mibi.id)) %>% nrow

# Set factor levels
patients$arm <- factor(patients$arm, levels = c("Combination", "Ipilimumab"))
patients$resp <- factor(patients$resp, levels = c("CR","PR","SD","PD"))
patients$resp2 <- factor(ifelse(patients$resp %in% c('CR','PR'), 'CR/PR', 'SD/PD'), levels = c('CR/PR','SD/PD'))
samples$timepoint <- factor(samples$timepoint, levels = c('baseline','prior to pd1','ontx c2'))

# Annotation of MIBI
mibi <- samples %>% filter(!is.na(mibi.id)) %>% left_join(patients)
mibi <- mibi %>%
  mutate(group = ifelse(arm == "Ipilimumab", "Ipi PD", ifelse(resp2 == "CR/PR", "Comb. CR/PR", "Comb. PD")))
mibi$group <- factor(mibi$group, levels = c("Comb. CR/PR", "Comb. PD", "Ipi PD"))

mibi %>% group_by(timepoint) %>% count
mibi %>% group_by(arm, resp2, timepoint) %>% count %>% spread(timepoint, n)
mibi %>% group_by(arm, resp2, timepoint, subject) %>% summarise %>% 
  group_by(arm, resp2, timepoint) %>%
  count %>% spread(timepoint, n)
mibi %>% group_by(arm, resp2, subject) %>% summarise %>% 
  group_by(arm, resp2) %>% count

mibi_ <- mibi %>% select(mibi.id, subject, timepoint, group)
order_x <- mibi %>% arrange(group, timepoint) %>% group_by(group, timepoint) %>% summarise %>%
  ungroup %>% mutate(x = 1:n())
missing <- mibi_ %>% group_by(subject, group, timepoint) %>% summarise(tmp = 1) %>% 
  spread(timepoint, tmp) %>% gather(timepoint, missing, -subject, -group) %>% 
  filter(is.na(missing)) %>% left_join(order_x)

# Segments
rois <- fread("mibi_roi_map.tsv")
segs <- fread("mibi_segments.tsv")
segs <- segs %>% left_join(rois) %>% left_join(mibi_)

# Cell type densities
tissue_dens <- segs %>% group_by(subject, timepoint, group) %>% summarise(area_px = sum(area))
tissue_dens <- tissue_dens %>% ungroup %>% mutate(area_mm2 = area_px * (.800/2048)**2)

# type order
types <- tibble(type = c('Tumor','CD8 T cell','CD4 T cell','Treg','B cell',
           'NK cell','NKT cell','Dendritic','Neutrophil','Plasma cell',
           'Macrophage','Monocyte','Endothelial','Fibroblast'))
#
type_dens <- segs %>% 
  filter(type %in% types$type) %>%
  group_by(subject, timepoint, group, type) %>% count %>% 
  left_join(tissue_dens) %>% rowwise %>% mutate(dens = n/area_mm2)
type_dens <- tissue_dens %>% cross_join(types) %>%
  full_join(type_dens) %>%
  mutate(dens = ifelse(is.na(dens), 0, dens)) %>%
  left_join(order_x)
  # filter(type %in% c('Tumor','CD8 T cell','CD4 T cell','Treg'))
type_dens$type <- factor(type_dens$type, levels = types$type)

dens_max <- type_dens %>% group_by(type) %>% summarise(y = max(dens))
comparisons <- tibble(group1 = as.character(c(1,3,5,1)), 
                      group2 = as.character(c(2,4,6,3)),
                      plus = 0.075*(c(1,2,1,3)))
stats <- compare_means(dens ~ x, data = type_dens, group.by = 'type') %>%
  inner_join(comparisons) %>% left_join(dens_max) %>%
  mutate(p.format = signif(p, 2))

ggplot(type_dens, aes(x = x,
                      y = dens, 
                      fill = group)) +
  geom_blank(data = dens_max, aes(x = NULL, y = y*0.35, fill = NULL)) +
  geom_boxplot(aes(alpha = timepoint), outlier.colour = NA, size = 0.25) +
  geom_line(aes(group = subject), size = 0.25) +
  ggbeeswarm::geom_beeswarm(shape = 21, colour = 'black', stroke = 0.1, size = 1) +
  #
  geom_text(data = stats, aes(x = (as.numeric(group1)+as.numeric(group2))/2,
                              y = y+y*(plus+0.01),
                              label = p.format),
            inherit.aes = FALSE, size = 5*0.36, vjust = 0) +
  geom_segment(data = stats, aes(x = as.numeric(group1), xend = as.numeric(group2),
                                 y = y+y*(plus), yend = y+y*(plus)),
               inherit.aes = FALSE, linewidth = 0.1, lineend = 'square') +
  geom_segment(data = stats, aes(x = as.numeric(group1), xend = as.numeric(group1),
                                 y = y+y*(plus-0.01), yend = y+y*(plus)),
               inherit.aes = FALSE, linewidth = 0.1, lineend = 'square') +
  geom_segment(data = stats, aes(x = as.numeric(group2), xend = as.numeric(group2),
                                 y = y+y*(plus-0.01), yend = y+y*(plus)),
               inherit.aes = FALSE, linewidth = 0.1, lineend = 'square') +
  #
  scale_x_continuous(breaks = c(1.5, 3.5, 5.5), labels = gsub(" ", "\n", levels(type_dens$group))) +
  scale_alpha_manual(values = c(0.9, 0.45)) +
  scale_fill_manual(values = group_colors) +
  guides(alpha = guide_legend(override.aes = list(fill = 'black'))) +
  labs(y = "Cell Type Density (cells per mm2 tissue)",
       fill = "Patient subset", alpha = "Timepoint") +
  facet_wrap(~ type, scales = 'free_y', nrow = 3) +
  theme_bw() +
  theme(strip.background = element_blank(),
        text = element_text(size = 5), title = element_text(size = 7),
        legend.key.height = unit(0.25, 'line'), legend.key.width = unit(0.25, 'line'),
        axis.title.x = element_blank()) +
  ae
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figures3.pdf", height = 5, width = 7.2, unit = 'in', family="ArialMT")
#

select_dens <- type_dens %>% 
  filter(type %in% c('Tumor','CD8 T cell','CD4 T cell','Treg','Monocyte'),
         group == "Comb. CR/PR") %>%
  group_by(type) %>% mutate(z = scale(dens)[,1])
change_dens <- select_dens %>% pivot_wider(id_cols = c('subject','group','type'),
                            names_from = c('x'), values_from = 'dens') %>%
  rowwise %>% mutate(change = log2(`2`/`1`)) %>%
  arrange(type, desc(`2`), change, desc(`1`))
select_dens <- select_dens %>% left_join(change_dens)
missing_ <- missing %>% filter(x %in% c(1, 2))
left <- ggplot(select_dens, aes(x = x,
                        y = factor(subject, levels = rev(unique(change_dens$subject)))
                        )) +
  geom_point(aes(fill = z, size = cut(dens, 
                                      breaks = c(0, 250, 500, 750,
                                                 1000, 5000,10000))),
             shape = 21, colour = 'black', stroke = 0.1) +
  geom_point(data = missing_, aes(shape = "Paired sample N/A"),
             colour = 'black', stroke = 0.1, size = 2) +
  scale_x_continuous(limits = c(0.5, 2.5), breaks = c(1,2), labels = c('baseline','ontx c2')) +
  scale_fill_distiller(palette = 'RdBu', limits = c(-1.05, 1.05), oob = squish) +
  scale_shape_manual(values = c(4)) +
  scale_size_manual(values = c(0.5,1,1.5,3,4,5),
                    labels = c('<250','250-500','500-750','750-1,000','1,000-5,000','>5,000')) +
  facet_grid(. ~ type, scales = 'free_x', space = 'free_x') +
  labs(shape = '', size = "Cell density (cells/mm2)", fill = "Relative cell\ndensity (Z)",
       x = 'Timepoint', y = 'Patient (Comb. CR/PR)') +
  guides(fill = guide_colorbar(#title.position = "top", 
                               ticks.colour = 'black', ticks.linewidth = 0.1, 
                               frame.colour = 'black', frame.linewidth = 0.1, ),
         size = guide_legend(ncol = 2)) +
  theme_bw() + ae +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
right <- ggplot(select_dens, aes(x = change,
                        y = factor(subject, levels = rev(unique(change_dens$subject))),
                        )) +
  geom_vline(xintercept = 0, linetype = 'dotted', colour = 'grey80') +
  geom_point(aes(fill = type), shape = 23, colour = 'black', stroke = 0.1, size = 1.5) +
  geom_point(data = missing_, aes(x = 0, shape = "Paired sample N/A"),
             colour = 'black', stroke = 0.1, size = 2) +
  scale_x_continuous(limits = c(-7.5, 9)) +
  scale_shape_manual(values = c(4)) +
  scale_fill_manual(values = color_sorted) +
  labs(x = "Log2 Fold Change", shape = '') +
  theme_bw() + ae +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
left + right +
  plot_layout(guides = 'collect', widths = c(2,1)) & 
  theme(legend.position = 'none', plot.margin = unit(c(1,1,1,1), 'pt'), panel.spacing = unit(1, 'pt'))
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/_figure3a.pdf", height = 2.5, width = 3.75, unit = 'in', family="ArialMT")

legends <- lapply(list(left, right), function(gg){
  gg <- gg + theme(legend.position = 'bottom')
  get_legend(gg) %>% as_ggplot()
})
wrap_plots(legends, nrow = 1)
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/_figure3alegend.pdf", height = 2.5, width = 3.75, unit = 'in', family="ArialMT")

# MIBI cluster cleanup: Remove groups that were patient-specific
percluster <- segs %>% filter(!is.na(type_group)) %>%
  group_by(subject, leiden, type_group) %>% count
percluster <- split(x = percluster, f = percluster$type_group)
percluster <- lapply(percluster, function(tg){ 
  confmat <- tg %>% spread(leiden, n, fill = 0) %>% ungroup %>% select(-type_group) %>% 
    column_to_rownames("subject") %>% data.matrix
  confmat_p <- confmat/rowSums(confmat)
  
  confmat_zs <- lapply(1:ncol(confmat_p), function(col_i){
    p = confmat[,col_i]/rowSums(confmat)
    margin <- sqrt(mean(p)*(1-mean(p))/length(p))
    z_scores <- (p-mean(p)) / margin # evaluate for outliers
    return(z_scores)
  })
  confmat_zs <- enframe(confmat_zs, name = "label", value = "z") %>% 
    mutate(label = colnames(confmat_p), group = list(rownames(confmat_p))) %>%
    unnest(c(group, z))
  })
z_plots <- lapply(c('CD8','CD4','Myeloid'), function(tg){
  df <- percluster[[tg]] %>%
    group_by(label) %>% 
    mutate(keep = ifelse(any(z>5), "Remove", "Keep"))
  ggplot(df, aes(x = z, y = label, fill = keep)) +
    ggridges::geom_density_ridges(scale = 1, linewidth = 0.1) +
    geom_vline(xintercept = 5, linetype = 'dotted', colour = 'grey80') +
    geom_label(data = filter(df, z>5 & keep == "Remove"),
              aes(label = group), vjust = 0, size = 5*0.36, alpha = 0.5, colour = NA) +
    geom_label(data = filter(df, z>5 & keep == "Remove"),
               aes(label = group), vjust = 0, size = 5*0.36, fill = NA) +
    scale_x_continuous(limits = c(-5, 20)) +
    scale_fill_manual(values = c('grey','red')) +
    guides(fill = guide_legend(nrow = 1)) +
    labs(fill = "Keep subpopulation", y = "Leiden cluster", x = "Z", title = tg) +
    theme_bw() +
    theme(plot.margin = unit(c(2,2,2,2), 'pt'),legend.position = 'bottom', legend.direction = 'vertical',
          text = element_text(size = 7), title = element_text(size = 7),
          axis.ticks = element_line(size = 0.1), panel.border = element_rect(linewidth = 0.1, colour = 'black'))
})
wrap_plots(z_plots) +
  plot_layout(guides = 'collect', ncol = 1) &
  theme(legend.position = 'bottom')
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figures4left.pdf", height = 8.8, width = 2, unit = 'in', family="ArialMT")

keepgroups <- percluster %>% enframe(name = "type_group") %>% unnest(value) %>% 
  group_by(label) %>% mutate(keep = ifelse(any(z>5), FALSE, TRUE))
keep <- percluster %>% enframe(name = "type_group") %>% unnest(value) %>% 
  group_by(label) %>% summarise(keep = ifelse(any(z>5), FALSE, TRUE))
segs_tg <- segs %>% inner_join(keepgroups, by = c("subject" = "group", "leiden" = "label", "type_group"))
umaps <- lapply(c('CD8','CD4','Myeloid'), function(tg){
  df <- filter(segs_tg, type_group == tg, keep == TRUE) %>%
    arrange(leiden)
  means <- df %>% group_by(leiden) %>% 
    summarise(umap1 = mean(umap1), umap2 = mean(umap2)) %>%
    mutate(label = gsub(paste0(tg, "_"), "", leiden))
  ggplot(df, aes(x = umap1, y = umap2, colour = leiden, fill = leiden)) +
    ggrastr::geom_point_rast(alpha = 0.2, size = 0.25) +
    geom_label(data = means, aes(label = label), size = 5*0.36, colour = 'white') +
    scale_colour_manual(values = color_sorted) +
    scale_fill_manual(values = color_sorted) +
    labs(x = "UMAP1", y = "UMAP2", title = paste0(tg, " (", prettyNum(nrow(df), big.mark = ","), " cells)")) +
    theme_void() +
    theme(legend.position = 'none',
          axis.title.x = element_text(size = 5), 
          axis.title.y = element_text(size = 5, angle = 90), 
          plot.title = element_text(size = 5, hjust = 0.5))
  # ggsave(paste0("~/manuscript/Biopsy analysis of trial S1616/figures/_figure3", tg, ".pdf"), height = 1, width = 1, unit = 'in', family="ArialMT")
})

rgg <- fread("mibi_leidenrankedfeatures.tsv")
# write.table(rgg, file = "table_s7.tsv", sep = '\t', quote = F, row.names = F)
rgg <- rgg %>% inner_join(keep, by = c("leiden" = "label")) %>% filter(keep == TRUE)
features <- split(x = rgg, f = rgg$type_group)
pf <- grep("neighbors", unique(rgg$names), value = T, invert = TRUE)
nhf <- grep("neighbors", unique(rgg$names), value = T, invert = FALSE)
change_names <- tibble(names = c("HLA class 1 A, B, and C, Na-K-ATPase"),
                       label = c("HLA I/NaK-ATPase"))
#
exp <- fread("mibi_leideninputmeasurements.tsv")
segs_exp <- segs %>% select(roi, segment, type_group, leiden) %>% 
  inner_join(keep, by = c("leiden" = "label")) %>% filter(keep == TRUE) %>%
  inner_join(exp, by = c("roi" = "id", "segment"))
leiden_expz <- segs_exp %>% 
  group_by(type_group, leiden) %>% summarise_at(c(pf, nhf), mean) %>%
  group_by(type_group) %>% mutate_at(c(pf, nhf), ~(scale(.) %>% as.vector))
leiden_expz <- split(x = leiden_expz, f = leiden_expz$type_group)

hms <- lapply(c('CD8','CD4','Myeloid'), function(tg) {
  f_ <- leiden_expz[[tg]] %>% column_to_rownames("leiden") %>%
    ungroup %>% select(all_of(pf), all_of(nhf)) %>% data.matrix
  f_hc <- hclust(dist(f_, method = 'euclidean'), method = 'ave')
  row_order <- labels(as.dendrogram(f_hc))
  right_dendro = ggdendrogram(f_hc) +
    coord_flip(xlim = c(0.5, max(dendro_data(f_hc)$segments$x)+0.5),
               ylim = c(0, max(dendro_data(f_hc)$segments$y)*1.05),
               expand = 0) +
    theme_void()
  #
  pf_ <- leiden_expz[[tg]] %>%
    column_to_rownames("leiden") %>%
    ungroup %>% select(all_of(pf)) %>% data.matrix
  pf_hc <- hclust(dist(t(pf_), method = 'euclidean'), method = 'ave')
  pf_order <- labels(as.dendrogram(pf_hc))
  nhf_ <- leiden_expz[[tg]] %>%
    column_to_rownames("leiden") %>%
    ungroup %>% select(all_of(nhf)) %>% data.matrix
  nhf_hc <- hclust(dist(t(nhf_), method = 'euclidean'), method = 'ave')
  nhf_order <- labels(as.dendrogram(nhf_hc))
  df <- leiden_expz[[tg]] %>% gather(names, exp, -type_group, -leiden)
  df <- features[[tg]] %>% left_join(df) %>% mutate(feat = ifelse(names %in% pf, "Protein expression", "Neighboring cell types"))
  xorder <- tibble(names = c(pf_order, nhf_order)) %>%
    left_join(change_names) %>%
    mutate(label = ifelse(is.na(label), names, label),
           label = gsub(".*fracneighbors_", "", label))
  ggplot(df, aes(y = factor(leiden, levels = row_order),
                 x = factor(names, levels = c(pf_order, nhf_order)),
                 colour = ifelse(pvals_adj<0.05, "1", "2"),
                 fill = exp, size = cut(abs(logfoldchanges), c(0, 1, 2, 5, 10, 25)))) +
    geom_point(shape = 21, stroke = 0.1) +
    scale_fill_distiller(palette = 'RdBu', limits = c(-1.4, 1.4), oob = squish) +
    scale_size_manual(values = c(0.75,1.25, 1.75, 2, 3), 
                      labels = c("<1", "1-2", "2-5", "5-10", ">10")) +
    scale_colour_manual(values = c('red3','black'), labels = c("padj<0.05","n.s.")) +
    scale_x_discrete(breaks = xorder$names, labels = xorder$label) +
    facet_grid(. ~ feat, scales = 'free_x', space = 'free_x') +
    labs(x = "Feature", y = "Leiden cluster", fill = "Relative expression",
         size = "Abs. Log2 Fold Change", colour = "Significance") +
    theme_bw() + ae +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          strip.background = element_blank(), legend.position = 'bottom')
})
wrap_plots(hms) + plot_layout(ncol = 1, guides = 'collect') &
  theme(legend.position = 'none', axis.title = element_text(size = 5),
        plot.margin = unit(c(2,2,2,2), 'pt'), axis.ticks.length = unit(2, 'pt'),
        strip.text = element_text(size = 7), panel.grid = element_blank())
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figures4right.pdf", height = 8.10, width = 5.4, unit = 'in', family="ArialMT")

hm_legs <- hms[[1]] + 
  guides(fill = guide_colorbar(ticks.colour = 'black', ticks.linewidth = 0.1, 
                               frame.colour = 'black', frame.linewidth = 0.1,
                               theme = theme(legend.key.width = unit(5, 'line'), 
                                             legend.key.height = unit(0.75, 'line'))), 
         colour = guide_legend(nrow = 2),
         size = guide_legend(nrow = 3)) +
  theme(legend.justification = 'left', legend.position = 'bottom', legend.title.position = 'top')
get_legend(hm_legs) %>% as_ggplot()
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figures4rightlegend.pdf", height = 0.75, width = 5.4, unit = 'in', family="ArialMT")

# Spot analysis
prep_spot <- segs %>% left_join(keep, by = c("leiden" = "label")) %>%
  filter((keep == TRUE) | (type %in% c('Endothelial','Fibroblast','B cell','Plasma cell','Tumor'))) %>%
  mutate(type = ifelse(is.na(leiden), type, leiden))
filters <- list(quote(group == "Comb. CR/PR"),
                quote(group == "Comb. PD"),
                quote(timepoint == "baseline" & grepl("Comb.", group)),
                quote(timepoint == "ontx c2" & grepl("Comb.", group)))
results <- list(quote(ifelse(timepoint == "ontx c2", TRUE, FALSE)),
                quote(ifelse(timepoint == "ontx c2", TRUE, FALSE)),
                quote(ifelse(group == "Comb. CR/PR", TRUE, FALSE)),
                quote(ifelse(group == "Comb. CR/PR", TRUE, FALSE)))
comps <- c('Combo, CR/PR: Ontx v Base',
                'Combo, PD: Ontx v Base',
                'Base, Combo: CR/PR v PD',
                'Ontx, Combo: CR/PR v PD')
# Prepare comparisons for spot
prep_comparisons <- lapply(1:length(filters), function(i){
  prep_spot %>% filter(!! filters[[i]]) %>%
    mutate(x = x_um, y = y_um,
           result = !! results[[i]],
           PID = subject, id = roi) %>%
    select(PID, id, timepoint, group, x, y, type, result)
})

run_spot <- function(dat) {
  print("Choosing radii...")
  smallest.dim <- dat %>% group_by(id) %>%
    summarise(x.range = abs(max(x) - min(x)),
              y.range = abs(max(y) - min(y)),
              min = min(x.range, y.range))
  radii.ripleys.rule <- seq(0, 0.25 * min(smallest.dim$min), length.out = 100)
  
  # Create all combinations of immune cell types
  cell.types <- unique(dat %>% select(type) %>% unlist())
  cell.type.combos <- combn(cell.types, 2, simplify = FALSE)
  print(paste0("Number of cell types: ", length(cell.types)))
  print(paste0("Number of cell type combinations: ", length(cell.type.combos)))
  
  # Initialize a matrix to store the results
  colocalization.matrix <- 
    matrix(nrow = length(cell.types), ncol = length(cell.types),
           dimnames = list(cell.types, cell.types))
  
  # Save the metadata
  colocalization.metadata <-
    matrix(list(), nrow = length(cell.types), ncol = length(cell.types),
           dimnames = list(cell.types, cell.types))
  
  progress <- seq(1, length(cell.type.combos), 10)
  # Iterate through the cell type combinations
  for (i in 1:length(cell.type.combos)) {
    if(i %in% progress){
      print(paste0(signif(i*100/length(cell.type.combos)), "% complete."))
    }
    svMisc::progress(i/(length(cell.type.combos)/100))
    
    # Save the combination
    cell.type.i <- cell.type.combos[[i]]
    
    # Save the row and column index
    row.ind <- which(cell.types == cell.type.i[1])
    col.ind <- which(cell.types == cell.type.i[2])
    
    # Run SPOT
    res <- spot(data = dat, 
                radii = radii.ripleys.rule,
                outcome = "result", # outcome cannot be "outcome" lol
                use.K = FALSE, K.diff = FALSE, 
                model.type = "logistic", cell.type = cell.type.i,
                marked = TRUE, pick.roi = "average",
                print.progress = FALSE) 
    # pick.roi- patient may have been counted multiple times
    # try pick.roi = "average"
    
    # Save the result
    colocalization.matrix[row.ind, col.ind] <- res$overall.pval
    colocalization.metadata[row.ind, col.ind][[1]] <- res
  }
  
  # Colocalization with tumor cells
  colocalization.pvals <- gdata::unmatrix(colocalization.matrix)
  colocalization.pvals <- colocalization.pvals[!is.na(colocalization.pvals)]
  
  # Create a table with the results
  colocalization.matrix.df <- reshape2::melt(colocalization.matrix, na.rm = TRUE)
  colnames(colocalization.matrix.df) <- c("Cell Type 1", "Cell Type 2", "SPOT P-Value")
  colocalization.matrix.df$`SPOT Q-Value` <- p.adjust(colocalization.matrix.df$`SPOT P-Value`, method = "fdr")
  
  pvals <- lapply(cell.type.combos, function(i){ 
    if( is.null(colocalization.metadata[i[1], i[2]][[1]]$pval.df) ) {
      return(NULL)
    } else {
      df <- colocalization.metadata[i[1], i[2]][[1]]$pval.df %>%
        mutate(type1 = i[1], type2 = i[2])
      return(df)
    }
  }) %>%
    enframe %>% unnest(value) %>%
    mutate(pval.neg.log10 = -log10(pval))
  print("Done!")
  return(list("radii" = radii.ripleys.rule, "cell_type_combos" = cell.type.combos, "p_vals" = colocalization.matrix.df, "res" = pvals))
}
# res_spot <- lapply(prep_comparisons, run_spot)
# full_spot <- lapply(res_spot, function(sp) { sp$res }) %>% enframe(name = "comparison") %>% mutate(comparison = comps) %>%
#   unnest(value) %>% mutate(name = NULL) %>% arrange(comparison, pval)
# write.table(full_spot, file = "table_s8.tsv", sep = '\t', quote = FALSE, row.names = FALSE)

full_spot <- fread("table_s8.tsv")
full_spot %>% group_by(comparison, type1, type2) %>% summarise %>% 
  group_by(comparison) %>% count
full_spot %>% filter(pval<0.05) %>%
  group_by(comparison, type1, type2) %>% count %>%
  group_by(comparison) %>% count
full_spot %>% filter(pval<0.05) %>%
  gather(i, type, type1, type2) %>%
  group_by(comparison, type) %>% count %>%
  group_by(comparison) %>% count
#
report_spot <- full_spot %>% filter(comparison %in% comps[3:4])
report_spot <- split(report_spot, report_spot$comparison)
plot_spot <- lapply(names(report_spot), function(sp){
  dat <- report_spot[[sp]] %>% arrange(pval.neg.log10)
  labs <- dat %>% filter(pval<0.05) %>% mutate(dir = ifelse(coef>0, 1, -1)) %>%
    arrange(-abs(coef)) %>% group_by(dir, type1, type2) %>% slice_head(n=1) %>%
    group_by(dir) %>% slice_head(n=3) %>% 
    rowwise %>% mutate(label = paste0(type1, ":", type2))
  ns <- dat %>% filter(pval<0.05) %>% mutate(dir = ifelse(coef>0, 3, -3)) %>%
    group_by(dir, type1, type2) %>% count %>% group_by(dir) %>% count %>%
    mutate(hjust = ifelse(dir >0, 1, 0))
  ggplot(dat, aes(x = coef, y = pval.neg.log10,
                  alpha = factor(ifelse(pval<0.05 & coef>0, "1", ifelse(pval<0.05 & coef<0, "2", "3"))),
                  fill = factor(ifelse(pval<0.05 & coef>0, "1", ifelse(pval<0.05 & coef<0, "2", "3"))),
                  colour = factor(ifelse(pval<0.05 & coef>0, "1", ifelse(pval<0.05 & coef<0, "2", "3"))))) +
    ggrastr::geom_point_rast(shape = 21, size = 0.5, stroke = 0.1) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dotted', colour = 'black', size = 0.25) +
    scale_fill_manual(breaks = c("1", "2", "3"), values = c(as.character(resp2_colors), 'grey80'),
                      labels = c("p<0.05, Comb. CR/PR", "p<0.05, Comb.SD/PD", 'ns'), drop = FALSE) +
    scale_colour_manual(breaks = c("1", "2", "3"), values = c('black', 'black', 'grey80'),
                        labels = c("p<0.05, Comb. CR/PR", "p<0.05, Comb.SD/PD", 'ns'), drop = FALSE) +
    scale_alpha_manual(breaks = c("1", "2", "3"), values = c(0.8, 0.8, 0.2),
                        labels = c("p<0.05, Comb. CR/PR", "p<0.05, Comb.SD/PD", 'ns'), drop = FALSE) +
    ggrepel::geom_text_repel(data = labs, aes(label = label, hjust = ifelse(coef>0, 0, 1)), size = 5*0.36, alpha = 1, show.legend = FALSE) +
    geom_text(data = ns, aes(x = dir, y = 0.05, hjust = hjust, label = n), size = 5*0.36, colour = 'black', inherit.aes = NULL) +
    coord_cartesian(xlim = c(-3, 3), ylim = c(0, 1.75)) +
    labs(x = "Coefficient", y = "-log10 (p value)",
         fill = "Significance", colour = "Significance", alpha = "Significance") +
    theme_bw() + ae +
    theme(legend.key.width = unit(0.5, 'line'), legend.key.height = unit(0.5, 'line'))
})
wrap_plots(plot_spot) +
  plot_layout(nrow = 1) &
  theme(title = element_text(size = 5), legend.position = 'none',
        plot.margin = unit(c(2,2,2,2), 'pt'))
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/_figure3EF.pdf", height = 1.5, width = 3.5, unit = 'in', family="ArialMT")
get_legend(plot_spot[[2]] + guides(fill = guide_legend(override.aes = list(size = 1))) +theme(legend.title = element_text(size = 5), legend.position = 'bottom')) %>% as_ggplot()
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/_figure3EFlegend.pdf", height = 0.25, width = 3.5, unit = 'in', family="ArialMT")

## Networks in responding biopsies
net <- report_spot$`Ontx, Combo: CR/PR v PD` %>%
  filter(coef>0, pval<0.05) %>% group_by(type1, type2) %>% 
  summarise(n = n(), across(c('radius','coef'), list(median = median, mean = mean, sd = sd)))

# Convert the data frame to a graph object
g <- igraph::graph_from_data_frame(
  d = net,
  directed = FALSE
)

# Assign radius_mean as the edge weights
igraph::E(g)$weight <- net$radius_mean
# Perform clustering using the Louvain method
clusters <- igraph::cluster_louvain(g)
igraph::V(g)$cluster <- igraph::membership(clusters)

# Sort vertices by cluster manually
sorted_vertices <- order(igraph::V(g)$cluster)
g <- igraph::permute(g, sorted_vertices)

# Create custome graph plot
pts <- data.frame(v = igraph::V(g), cluster = igraph::V(g)$cluster) %>% rownames_to_column("type") %>% 
  mutate(network = paste0("Network ", cluster)) %>% select(type, network)
pts$network <- factor(pts$network, levels = paste0("Network ", c(2,1,4,5,3)))
#
edges <- net %>% ungroup %>% mutate(edge = 1:n()) %>%
  left_join(pts, by = c("type1" = "type")) %>% left_join(pts, by = c("type2" = "type"), suffix = c("1","2"))
horizontal_edges <- edges %>% filter(network1 == network2) %>%
  gather(i, type, type1, type2) %>% left_join(pts) %>% group_by(network, type) %>%
  count(name = "n_horizontal")
vertical_edges <- edges %>% filter(network1 != network2) %>%
  gather(i, type, type1, type2) %>% left_join(pts) %>% group_by(network, type) %>%
  count(name = "n_vertical")
weight_edges <- horizontal_edges %>% full_join(vertical_edges) %>% 
  mutate(n_vertical = ifelse(is.na(n_vertical), 0, n_vertical))
#
org_edges <- weight_edges %>% arrange(network, -n_horizontal, -n_vertical) %>%
  group_by(network) %>% mutate(x = (1:n()-1)*(-1)^(1:n())) %>%
  arrange(network, x) %>% mutate(x = 1:n()) %>% select(network, type, x)
push <- edges %>% filter(network1 != network2) %>% 
  left_join(org_edges, by = c("type1" = "type", "network1" = "network")) %>% 
  left_join(org_edges, by = c("type2" = "type", "network2" = "network"), suffix = c("1", "2")) %>%
  mutate(xdiff = x2-x1) %>% select(matches('type'), xdiff) %>% gather(i, type, -xdiff) %>%
  left_join(pts) %>% group_by(type, network) %>% summarise(push = mean(xdiff)) %>%
  group_by(network) %>% summarise(push = mean(push))
org_edges <- org_edges %>% left_join(push) %>% mutate(x = x + abs(push))
#
plot_edges <- edges %>% select(matches('type'), n, matches('network')) %>%
  left_join(org_edges, by = c("type1" = "type", "network1" = "network")) %>% 
  left_join(org_edges, by = c("type2" = "type", "network2" = "network"), suffix = c("1", "2"))
plot_edges$network1 <- factor(plot_edges$network1, levels = levels(org_edges$network))
plot_edges$network2 <- factor(plot_edges$network2, levels = levels(org_edges$network))
#
get_labs <- rgg %>%
  filter(scores>5, pvals_adj<0.05) %>%
  group_by(type_group, leiden) %>% slice_head(n = 3) %>%
  group_by(type_group, leiden) %>% summarise(n = n(), label = paste0(names, collapse = ", ")) %>%
  mutate(label = gsub("HLA class 1 A, B, and C, Na-K-ATPase", "HLA I", label)) %>%
  mutate(label = gsub("HLA DDQDR", "HLA II", label)) %>%
  mutate(label = gsub("fracneighbors_", "", label)) %>%
  mutate(label = gsub("PD-1", "PD1", label))
#
org_edges_ <- org_edges %>% left_join(get_labs, by = c("type" = "leiden"))

rects <- org_edges %>% mutate(y = as.numeric(network)) %>%
  group_by(network, y) %>% summarise(minx = min(x), maxx = max(x))

org_edges_$type_group <- factor(org_edges_$type_group, levels = c('CD8','CD4','Myeloid'))

## GRAPH DIAGRAM THIS ONE
ggplot(data = NULL, aes(y = as.numeric(network))) +
  # rectangles
  geom_text(data = rects, aes(x = minx-0.55, y = y, label = paste0("Louvain ", network)),
            hjust = 1.05, size = 5*0.36) +
  geom_rect(data = rects, aes(xmin = minx - 0.5, xmax = maxx + 0.5, ymin = y-0.375, ymax = y+0.375, fill = factor(network, levels = rev(levels(rects$network)))),
            inherit.aes = FALSE, colour = NA, alpha = 0.2) +
  # edges
  geom_curve(data = filter(plot_edges, network1 != network2),
             aes(y = as.numeric(network1), yend = as.numeric(network2), x = x1, xend = x2), 
             curvature = 0.1, alpha = 0.2, size = 0.1) +
  geom_curve(data = filter(plot_edges, network1 == network2),
             aes(y = as.numeric(network1), yend = as.numeric(network2), x = x1, xend = x2), 
             curvature = 0.2, alpha = 0.6, size = 0.1) +
  # vertices
  geom_point(data = org_edges_, aes(x = x, fill = factor(network, levels = rev(levels(rects$network))), 
                                    shape = type_group), size = 5, stroke = 0.05) +
  # cell labels
  geom_label(data = org_edges, aes(x = x, y = as.numeric(network)+0.1, label = gsub("Treg", "", type)),
             vjust = 0, size = 5*0.36, label.size = 0.05) +
  # descriptors
  geom_label(data = org_edges_, aes(x = x, y = as.numeric(network)-0.1, label = str_wrap(label, 20)), # 
             vjust = 1, size = 5*0.36, fill = 'white', alpha = 1, label.size = 0.05) +
  scale_y_continuous(breaks = 1:length(levels(pts$network)), labels = levels(pts$network)) +
  scale_x_continuous(limits = c(min(rects$minx-1.5), max(rects$maxx+0.5))) +
  scale_shape_manual(values = c(21, 22, 23), na.value = 25, labels = c("CD8 T cell subsets", "CD4 T cell subsets", "Myeloid subsets", "Other cell types")) +
  scale_fill_manual(values = color_sorted) +
  labs(shape = "Cell populations", fill = "Network\n(Louvain cluster)") +
  # guides(fill = guide_legend(override.aes = list(alpha = 1, shape = NA))) +
  guides(fill = 'none',
         shape = guide_legend(override.aes = list(stroke = 0.25))) +
  theme_void() +
  theme(text = element_text(size = 5), title = element_text(size = 5),
        legend.key.size = unit(c(0.5,0.5),'line'),
        legend.position = c(1.1,0.10), legend.direction = 'vertical', legend.justification = c(1,0),
        plot.margin = unit(c(2,2,2,2), 'pt'))
# ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/_figure3G.pdf", height = 2.5, width = 7.2, unit = 'in', family="ArialMT")

## Find examples of Louvain Networks 3 and 4
pops <- pts %>% filter(network == "Network 5")
net_ <- net %>% filter(type1 %in% pops$type & type2 %in% pops$type)
pops_ <- segs %>% 
  mutate(type = ifelse(is.na(leiden), type, leiden)) %>%
  filter(type %in% pops$type,
         timepoint == "ontx c2",
         group == "Comb. CR/PR")
pops_ <- split(pops_, pops_$roi)
nets_ <- lapply(pops_, function(p_){
  t_ <- split(p_, p_$type)
  pp_ <- lapply(t_, function(t_i) { ppp(t_i$x_um, t_i$y_um, window = owin(c(0, 800), c(0, 800))) })
  apply(net_, 1, function(i) {
    if ( ! (i['type1'] %in% names(pp_) & i['type2'] %in% names(pp_) ) ) {
      return(NULL)
    } else if ( nrow(t_[[i['type1']]]) > 0 & nrow(t_[[i['type2']]])>0 ) {
      cd <- crossdist(pp_[[i['type1']]], pp_[[i['type2']]])
      # ind <- which(cd<i['radius_mean'], arr.ind = TRUE)
      ind <- which(cd<30, arr.ind=TRUE)
      # if (i['type1'] == "Myeloid_5" | i['type2'] == "Myeloid_5") {
      #   ind <- which(cd<30, arr.ind=TRUE)
      # } else {
      #   ind <- which(cd<15, arr.ind=TRUE) 
      # }
      t1 <- t_[[i['type1']]] %>% mutate(i1 = 1:n()) %>% select(segment, i1)
      t2 <- t_[[i['type2']]] %>% mutate(i2 = 1:n()) %>% select(segment, i2)
      t12 <- tibble(type1 = i['type1'], type2 = i['type2'], i1 = ind[,'row'], i2 = ind[,'col']) %>%
        left_join(t1, by = c("i1")) %>% left_join(t2, by = c("i2"), suffix = c("1", "2"))
    } else {
      return(NULL)
    }
  }) %>% enframe %>% unnest(value)
})
network3 <- nets_ %>% enframe(name = "roi") %>% unnest(value)
network3 %>% group_by(roi, type1, type2) %>% count %>% unite(types, type1, type2, sep = ":") %>%
  spread(types, n, fill = 0)
network3 %>% group_by(roi) %>% count %>% arrange(-n)

# roi_ = "PT0524_060118-BX_ROI1"#"PT0688_091019-BX-1_ROI1"#"PT0524_060118-BX_ROI1"#"PT0679_082119-BX-2_ROI3"#PT0683_090919-BX-1_ROI3"#
# choose_ = c("PT0688_091019-BX-1_ROI1","PT0524_060118-BX_ROI1","PT0679_082119-BX-2_ROI3","PT0683_090919-BX-1_ROI3")
# c("PT0694_101419-BX-3-3_ROI2","PT0689_092419-BX-1_ROI4","PT0683_090919-BX-1_ROI3","PT0524_060118-BX_ROI1")
# unique(network3$roi)
# print_rois <- 
lapply(c("PT0621_022519-BXA1_ROI6"), function(roi_){
  roi_name = str_replace(roi_, "_", "-")
  outline_path <- "~/Google Drive/My Drive/Morrison/MORRISON ptII- Cohort 3/morrison-cohort-3-mibi/structures/"
  outlines <- list.files(path = outline_path, pattern = "dsDNA.+geojson")
  tumors <- list.files(path = outline_path, pattern = "SOX10.+geojson")
  endothelials <- list.files(path = outline_path, pattern = "CD31.+geojson")
  roi_outline <- paste0(outline_path, grep(str_replace(roi_, "_", "-"), outlines, value = TRUE))
  roi_tumor <- paste0(outline_path, grep(str_replace(roi_, "_", "-"), tumors, value = TRUE))
  roi_endoth <- paste0(outline_path, grep(str_replace(roi_, "_", "-"), endothelials, value = TRUE))
  outline_geom <- st_read(roi_outline) %>% st_geometry()
  tumor_geom <- st_read(roi_tumor) %>% st_geometry()
  endoth_geom <- st_read(roi_endoth) %>% st_geometry()
  #
  outline_sfc <- lapply(outline_geom, st_polygon) %>% st_as_sfc()
  tumor_sfc <- lapply(tumor_geom, st_polygon) %>% st_as_sfc()
  # tumor_sfc <- tumor_sfc[st_is_valid(tumor_sfc)]
  tumor_inoutline <- st_intersection(tumor_sfc, outline_sfc)#tumor_sfc#
  tumor_inoutline <- st_cast(tumor_inoutline, "MULTIPOLYGON")
  endoth_sfc <- lapply(endoth_geom, st_polygon) %>% st_as_sfc()
  endoth_sfc <- endoth_sfc[st_is_valid(endoth_sfc)]
  endoth_sfc <- st_cast(endoth_sfc, "POLYGON")
  endoth_inoutline <- st_intersection(endoth_sfc, outline_sfc)
  endoth_inoutline <- st_cast(endoth_inoutline, "MULTIPOLYGON")
  #
  outline_poly <- lapply(outline_sfc, function(g){
    outline <- st_coordinates(g)
    ## reverse Y coords
    outline[,"Y"] <- outline[,"Y"]*-1
    return(outline)
  })
  tumor_poly <- lapply(tumor_inoutline, function(g){
    tumor <- st_coordinates(g)
    ## reverse Y coords
    tumor[,"Y"] <- tumor[,"Y"]*-1
    return(tumor)
  })
  endoth_poly <- lapply(endoth_inoutline, function(g){
    endoth <- st_coordinates(g)
    ## reverse Y coords
    endoth[,"Y"] <- endoth[,"Y"]*-1
    return(endoth)
  })

  plot_networks <- network3 %>% filter(roi == roi_) %>% mutate(network3 = TRUE) %>%
    left_join(segs, by = c("roi", "segment1" = "segment")) %>%
    left_join(segs, by = c("roi", "segment2" = "segment"), suffix = c("_1", "_2"))
  plot_network_cells <- network3 %>% filter(roi == roi_) %>% mutate(network3 = TRUE) %>%
    gather(i, segment, segment1, segment2) %>% select(roi, segment) %>% distinct %>%
    inner_join(segs) %>% mutate(type = ifelse(is.na(leiden), type, leiden))
  # plot_others <- segs %>% filter(roi == roi_) %>% filter(! segment %in% c(plot_networks$segment1, plot_networks$segment2))

  tumor_polygg <- lapply(tumor_poly, function(poly) {
    geom_polygon(data = poly, aes(x = X, y = Y, fill = '_A'), colour = NA)
  })
  endoth_polygg <- lapply(endoth_poly, function(poly) {
    geom_polygon(data = poly, aes(x = X, y = Y, fill = '_B'), colour = NA)
  })
  outline_polygg <- lapply(outline_poly, function(poly) {
    geom_polygon(data = poly, aes(x = X, y = Y), fill = NA, colour = 'black', linewidth = 0.1)
  })
  p <- ggplot(data = NULL) +
    tumor_polygg +
    endoth_polygg +
    outline_polygg +
    geom_segment(data = plot_networks, aes(x = x_um_1, xend = x_um_2, y = -y_um_1, yend = -y_um_2),
                 size = 0.1) +
    geom_point(data = plot_network_cells, aes(x = x_um, y = -y_um, fill = type),
               size = 0.5, shape = 21, colour = 'black', stroke = 0.05) +
    scale_fill_manual(values = c('grey90','#efbcb8',color_sorted[c(1:4,11,8)]),
                      breaks = c('_A','_B',sort(pops$type)),
                      labels = c("Tumor region", "Endothelial region",
                                 sort(pops$type))) +
    scale_y_continuous(breaks = seq(0, -800, -200),
                       labels = c("0", "", "400", "", "800"),
                       sec.axis = sec_axis(~ ., breaks = seq(0,-800,-200), labels = c("0", "", "400", "", "800"))) +
    scale_x_continuous(breaks = seq(0, 800, 200), labels = paste0(seq(0, 800, 200), " um"),
                       sec.axis = sec_axis(~ ., breaks = seq(0, 800, 200))) +
    coord_fixed(xlim = c(0, 800), ylim = c(-800, 0), expand = 0) +
    labs(fill = "Annotation", x = roi_name) +
    guides(fill = guide_legend(nrow = 2)) +
    theme_bw() +
    theme(text = element_text(size = 5), title = element_text(size = 5),
          axis.title.y = element_blank(),
          axis.title.x = element_text(hjust = 0),
          axis.text.y.left = element_text(hjust = 1.3), axis.text.x = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks = element_line(size = 0.05),
          panel.border = element_rect(linewidth = 0.1, colour = 'black'),
          panel.grid = element_blank(), plot.margin = unit(c(1,1,1,1), 'pt'),
          legend.position = 'none')
  # legend.position = 'bottom', legend.title.position = 'top')
  p
  ggsave(paste0("../figures/3h/n5_figure3.",roi_,".pdf"), height = 1.5, width = 1.5, unit = 'in', family="ArialMT")
  return(p)
})

# 

print_legend <- lapply(c("PT0524_060118-BX_ROI1"), function(roi_){
  roi_name = str_replace(roi_, "_", "-")
  outline_path <- "~/Google Drive/My Drive/Morrison/MORRISON ptII- Cohort 3/morrison-cohort-3-mibi/structures/"
  outlines <- list.files(path = outline_path, pattern = "dsDNA.+geojson")
  tumors <- list.files(path = outline_path, pattern = "SOX10.+geojson")
  endothelials <- list.files(path = outline_path, pattern = "CD31.+geojson")
  roi_outline <- paste0(outline_path, grep(str_replace(roi_, "_", "-"), outlines, value = TRUE))
  roi_tumor <- paste0(outline_path, grep(str_replace(roi_, "_", "-"), tumors, value = TRUE))
  roi_endoth <- paste0(outline_path, grep(str_replace(roi_, "_", "-"), endothelials, value = TRUE))
  outline_geom <- st_read(roi_outline) %>% st_geometry()
  tumor_geom <- st_read(roi_tumor) %>% st_geometry()
  endoth_geom <- st_read(roi_endoth) %>% st_geometry()
  #
  outline_sfc <- lapply(outline_geom, st_polygon) %>% st_as_sfc()
  tumor_sfc <- lapply(tumor_geom, st_polygon) %>% st_as_sfc()
  # tumor_sfc <- tumor_sfc[st_is_valid(tumor_sfc)]
  tumor_inoutline <- st_intersection(tumor_sfc, outline_sfc)#tumor_sfc#
  tumor_inoutline <- st_cast(tumor_inoutline, "MULTIPOLYGON")
  endoth_sfc <- lapply(endoth_geom, st_polygon) %>% st_as_sfc()
  endoth_sfc <- endoth_sfc[st_is_valid(endoth_sfc)]
  endoth_sfc <- st_cast(endoth_sfc, "POLYGON")
  endoth_inoutline <- st_intersection(endoth_sfc, outline_sfc)
  endoth_inoutline <- st_cast(endoth_sfc, "MULTIPOLYGON")
  #
  outline_poly <- lapply(outline_sfc, function(g){
    outline <- st_coordinates(g)
    ## reverse Y coords
    outline[,"Y"] <- outline[,"Y"]*-1
    return(outline)
  })
  tumor_poly <- lapply(tumor_inoutline, function(g){
    tumor <- st_coordinates(g)
    ## reverse Y coords
    tumor[,"Y"] <- tumor[,"Y"]*-1
    return(tumor)
  })
  endoth_poly <- lapply(endoth_inoutline, function(g){
    endoth <- st_coordinates(g)
    ## reverse Y coords
    endoth[,"Y"] <- endoth[,"Y"]*-1
    return(endoth)
  })
  
  plot_networks <- network3 %>% filter(roi == roi_) %>% mutate(network3 = TRUE) %>%
    left_join(segs, by = c("roi", "segment1" = "segment")) %>%
    left_join(segs, by = c("roi", "segment2" = "segment"), suffix = c("_1", "_2"))
  plot_network_cells <- network3 %>% filter(roi == roi_) %>% mutate(network3 = TRUE) %>%
    gather(i, segment, segment1, segment2) %>% select(roi, segment) %>% distinct %>%
    inner_join(segs) %>% mutate(type = ifelse(is.na(leiden), type, leiden))
  # plot_others <- segs %>% filter(roi == roi_) %>% filter(! segment %in% c(plot_networks$segment1, plot_networks$segment2))
  
  tumor_polygg <- lapply(tumor_poly, function(poly) { 
    geom_polygon(data = poly, aes(x = X, y = Y, fill = '_A'), colour = NA) 
  })
  endoth_polygg <- lapply(endoth_poly, function(poly) { 
    geom_polygon(data = poly, aes(x = X, y = Y, fill = '_B'), colour = NA) 
  })
  outline_polygg <- lapply(outline_poly, function(poly) { 
    geom_polygon(data = poly, aes(x = X, y = Y), fill = NA, colour = 'black', linewidth = 0.1)
  })
  p <- ggplot(data = NULL) +
    tumor_polygg +
    endoth_polygg +
    outline_polygg +
    geom_segment(data = plot_networks, aes(x = x_um_1, xend = x_um_2, y = -y_um_1, yend = -y_um_2), 
                 size = 0.1) +
    geom_point(data = plot_network_cells, aes(x = x_um, y = -y_um, fill = type), 
               size = 0.5, shape = 21, colour = 'black', stroke = 0.05) +
    scale_fill_manual(values = c('grey90','#efbcb8',color_sorted[c(1:4,11,8)]),
                      breaks = c('_A','_B',sort(pops$type)),
                      labels = c("Tumor region", "Endothelial region",
                                 sort(pops$type))) +
    scale_y_continuous(breaks = seq(0, -800, -200), 
                       labels = c("0", "", "400", "", "800"),
                       sec.axis = sec_axis(~ ., breaks = seq(0,-800,-200), labels = c("0", "", "400", "", "800"))) +
    scale_x_continuous(breaks = seq(0, 800, 200), labels = paste0(seq(0, 800, 200), " um"),
                       sec.axis = sec_axis(~ ., breaks = seq(0, 800, 200))) +
    coord_fixed(xlim = c(0, 800), ylim = c(-800, 0), expand = 0) +
    labs(fill = "Annotation", x = roi_name) +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size = 3))) +
    theme_bw() +
    theme(text = element_text(size = 5), title = element_text(size = 5),
          axis.title.y = element_blank(),
          axis.title.x = element_text(hjust = 0),
          axis.text.y.left = element_text(hjust = 1.3), axis.text.x = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks = element_line(size = 0.05), 
          panel.border = element_rect(linewidth = 0.1, colour = 'black'),
          panel.grid = element_blank(), plot.margin = unit(c(1,1,1,1), 'pt'),
          legend.key.size = unit(c(0.5,0.5), 'line'), legend.title.position = 'top')
  return(p)
})

# get_legend((print_legend[[1]] + theme(legend.justification = c(0,1)))) %>% as_ggplot()
# ggsave(paste0("../figures/_figure3.legend.pdf"), height = 1.5, width = 1.25, unit = 'in', family="ArialMT")

# export <- network3 %>% filter(roi == "PT0683_090919-BX-1_ROI3") %>% select(segment1, type1)
# export <- network3 %>% filter(roi == "PT0683_090919-BX-1_ROI3") %>% select(segment2, type2) %>%
#   full_join(export, by = c("segment2" = "segment1", "type2" = "type1")) %>% distinct
# write.table(export, file = "~/Downloads/export.csv", sep = ',', quote = F, col.names = F, row.names = F)

## Print out all series of networks in images

## Print all networks
network_segs <- lapply(unique(pts$network), function(netname){
  pops <- pts %>% filter(network == netname)
  net_ <- net %>% filter(type1 %in% pops$type & type2 %in% pops$type)
  pops_ <- segs %>% 
    mutate(type = ifelse(is.na(leiden), type, leiden)) %>%
    filter(type %in% pops$type,
           timepoint == "ontx c2",
           group == "Comb. PD")
  # ,
  # timepoint == "ontx c2",
  # group == "Comb. CR/PR"
  pops_ <- split(pops_, pops_$roi)
  nets_ <- lapply(pops_, function(p_){
    t_ <- split(p_, p_$type)
    pp_ <- lapply(t_, function(t_i) { ppp(t_i$x_um, t_i$y_um, window = owin(c(0, 800), c(0, 800))) })
    apply(net_, 1, function(i) {
      if ( ! (i['type1'] %in% names(pp_) & i['type2'] %in% names(pp_) ) ) {
        return(NULL)
      } else if ( nrow(t_[[i['type1']]]) > 0 & nrow(t_[[i['type2']]])>0 ) {
        cd <- crossdist(pp_[[i['type1']]], pp_[[i['type2']]])
        ind <- which(cd < (ceiling(as.numeric(i['radius_mean']) / 5) * 5), arr.ind = TRUE)
        t1 <- t_[[i['type1']]] %>% mutate(i1 = 1:n()) %>% select(segment, i1)
        t2 <- t_[[i['type2']]] %>% mutate(i2 = 1:n()) %>% select(segment, i2)
        t12 <- tibble(type1 = i['type1'], type2 = i['type2'], i1 = ind[,'row'], i2 = ind[,'col']) %>%
          left_join(t1, by = c("i1")) %>% left_join(t2, by = c("i2"), suffix = c("1", "2"))
      } else {
        return(NULL)
      }
    }) %>% enframe %>% unnest(value)
  })
  network <- nets_ %>% enframe(name = "roi") %>% unnest(value) %>% mutate(value = NULL)
  return(network)
})
network_segs <- network_segs %>% enframe(name = "network") %>% mutate(network = unique(pts$network)) %>%
  unnest(value)


net_order <- paste0("Network ", c(3, 5, 4, 1, 2))

# specials <- c("PT0688_091019-BX-1_ROI1")
# unique(network_segs$roi)
tmp <- sig_results %>% filter(timepoint == "ontx c2", group == "Comb. CR/PR")
# print_rois <- lapply(tmp$id, function(roi_){
# lapply(c("PT0694_101419-BX-3-3_ROI2"), function(roi_) {
lapply(c("PT0621_022519-BXA1_ROI6"), function(roi_) {
  print(roi_)
  roi_name = str_replace(roi_, "_", "-")
  outline_path <- "~/Google Drive/My Drive/Morrison/MORRISON ptII- Cohort 3/morrison-cohort-3-mibi/structures/"
  outlines <- list.files(path = outline_path, pattern = "dsDNA.+geojson")
  tumors <- list.files(path = outline_path, pattern = "SOX10.+geojson")
  endothelials <- list.files(path = outline_path, pattern = "CD31.+geojson")
  roi_outline <- paste0(outline_path, grep(str_replace(roi_, "_", "-"), outlines, value = TRUE))
  roi_tumor <- paste0(outline_path, grep(str_replace(roi_, "_", "-"), tumors, value = TRUE))
  roi_endoth <- paste0(outline_path, grep(str_replace(roi_, "_", "-"), endothelials, value = TRUE))
  outline_geom <- st_read(roi_outline) %>% st_geometry()
  tumor_geom <- st_read(roi_tumor) %>% st_geometry()
  endoth_geom <- st_read(roi_endoth) %>% st_geometry()
  #
  outline_sfc <- lapply(outline_geom, st_polygon) %>% st_as_sfc()
  tumor_sfc <- lapply(tumor_geom, st_polygon) %>% st_as_sfc()
  # tumor_sfc <- tumor_sfc[st_is_valid(tumor_sfc)]
  tumor_inoutline <- st_intersection(tumor_sfc, outline_sfc)#tumor_sfc#
  tumor_inoutline <- st_cast(tumor_inoutline, "MULTIPOLYGON")
  endoth_sfc <- lapply(endoth_geom, st_polygon) %>% st_as_sfc()
  endoth_sfc <- endoth_sfc[st_is_valid(endoth_sfc)]
  endoth_sfc <- st_cast(endoth_sfc, "POLYGON")
  endoth_inoutline <- st_intersection(endoth_sfc, outline_sfc)
  endoth_inoutline <- st_cast(endoth_inoutline, "MULTIPOLYGON")
  #
  outline_poly <- lapply(outline_sfc, function(g){
    outline <- st_coordinates(g)
    ## reverse Y coords
    outline[,"Y"] <- outline[,"Y"]*-1
    return(outline)
  })
  tumor_poly <- lapply(tumor_inoutline, function(g){
    tumor <- st_coordinates(g)
    ## reverse Y coords
    tumor[,"Y"] <- tumor[,"Y"]*-1
    return(tumor)
  })
  endoth_poly <- lapply(endoth_inoutline, function(g){
    endoth <- st_coordinates(g)
    ## reverse Y coords
    endoth[,"Y"] <- endoth[,"Y"]*-1
    return(endoth)
  })
  
  plot_networks <- network_segs %>% filter(roi == roi_) %>% 
    left_join(segs, by = c("roi", "segment1" = "segment")) %>%
    left_join(segs, by = c("roi", "segment2" = "segment"), suffix = c("_1", "_2"))
  plot_network_cells <- network_segs %>% filter(roi == roi_) %>% 
    gather(i, segment, segment1, segment2) %>% select(roi, segment, network) %>% distinct %>%
    inner_join(segs)
  plot_plasma <- segs %>%
    filter(roi == roi_, type == "Plasma cell")
  
  tumor_polygg <- lapply(tumor_poly, function(poly) {
    geom_polygon(data = poly, aes(x = X, y = Y, fill = '_A'), colour = NA)
  })
  endoth_polygg <- lapply(endoth_poly, function(poly) {
    geom_polygon(data = poly, aes(x = X, y = Y, fill = '_B'), colour = NA)
  })
  outline_polygg <- lapply(outline_poly, function(poly) {
    geom_polygon(data = poly, aes(x = X, y = Y), fill = NA, colour = 'black', linewidth = 0.1)
  })
  p <- ggplot(data = NULL) +
    tumor_polygg +
    endoth_polygg +
    outline_polygg +
    geom_segment(data = plot_networks, aes(x = x_um_1, xend = x_um_2, y = -y_um_1, yend = -y_um_2),
                 size = 0.1) +
    geom_point(data = plot_network_cells, aes(x = x_um, y = -y_um, fill = network),
               size = 0.5, shape = 21, colour = 'black', stroke = 0.05) +
    # geom_point(data = plot_plasma, aes(x = x_um, y = -y_um, fill = type), 
    #            size = 0.75, shape = 21, colour = 'black', stroke = 0.10) +
    scale_fill_manual(values = c('grey90','#efbcb8',color_sorted, "red3"),
                      breaks = c('_A','_B',net_order, "Plasma cell"),
                      labels = c("Tumor region", "Endothelial region",
                                 net_order, "Plasma cell"),
                      drop = FALSE) +
    scale_y_continuous(breaks = seq(0, -800, -200),
                       labels = c("0", "", "400", "", "800"),
                       sec.axis = sec_axis(~ ., breaks = seq(0,-800,-200), labels = c("0", "", "400", "", "800"))) +
    scale_x_continuous(breaks = seq(0, 800, 200), labels = paste0(seq(0, 800, 200), " um"),
                       sec.axis = sec_axis(~ ., breaks = seq(0, 800, 200))) +
    coord_fixed(xlim = c(0, 800), ylim = c(-800, 0), expand = 0) +
    labs(fill = "Annotation", x = roi_name) +
    guides(fill = guide_legend(nrow = 2)) +
    theme_bw() +
    theme(text = element_text(size = 5), title = element_text(size = 5),
          axis.title.y = element_blank(),
          axis.title.x = element_text(hjust = 0),
          axis.text.y.left = element_text(hjust = 1.3), axis.text.x = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks = element_line(size = 0.05),
          panel.border = element_rect(linewidth = 0.1, colour = 'black'),
          panel.grid = element_blank(), plot.margin = unit(c(1,1,1,1), 'pt'),
          legend.position = 'none')
  p
  # ggsave(paste0("../figures/s5/figures5_plasma.",roi_,".pdf"), height = 1.5, width = 1.5, unit = 'in', family="ArialMT")
  ggsave(paste0("../figures/3h/alln_every.",roi_,".pdf"), height = 1.5, width = 1.5, unit = 'in', family="ArialMT")
  return(p)
})

(print_rois[[1]] +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = 'bottom', legend.margin = margin(0,0,0,0,'pt'),
        legend.justification = c(0,1), legend.key.size = unit(c(0.5,0.5), 'line'))) %>%
  get_legend() %>% as_ggplot()
# ggsave(paste0("../figures/s5/legend.pdf"), height = 1.5, width = 4, unit = 'in', family="ArialMT")

# 

print_legend <- lapply(c("PT0524_060118-BX_ROI1"), function(roi_){
  roi_name = str_replace(roi_, "_", "-")
  print(roi_)
  roi_name = str_replace(roi_, "_", "-")
  outline_path <- "~/Google Drive/My Drive/Morrison/MORRISON ptII- Cohort 3/morrison-cohort-3-mibi/structures/"
  outlines <- list.files(path = outline_path, pattern = "dsDNA.+geojson")
  tumors <- list.files(path = outline_path, pattern = "SOX10.+geojson")
  endothelials <- list.files(path = outline_path, pattern = "CD31.+geojson")
  roi_outline <- paste0(outline_path, grep(str_replace(roi_, "_", "-"), outlines, value = TRUE))
  roi_tumor <- paste0(outline_path, grep(str_replace(roi_, "_", "-"), tumors, value = TRUE))
  roi_endoth <- paste0(outline_path, grep(str_replace(roi_, "_", "-"), endothelials, value = TRUE))
  outline_geom <- st_read(roi_outline) %>% st_geometry()
  tumor_geom <- st_read(roi_tumor) %>% st_geometry()
  endoth_geom <- st_read(roi_endoth) %>% st_geometry()
  #
  outline_sfc <- lapply(outline_geom, st_polygon) %>% st_as_sfc()
  tumor_sfc <- lapply(tumor_geom, st_polygon) %>% st_as_sfc()
  # tumor_sfc <- tumor_sfc[st_is_valid(tumor_sfc)]
  tumor_inoutline <- tumor_sfc#st_intersection(tumor_sfc, outline_sfc)#tumor_sfc#
  tumor_inoutline <- st_cast(tumor_inoutline, "MULTIPOLYGON")
  endoth_sfc <- lapply(endoth_geom, st_polygon) %>% st_as_sfc()
  endoth_sfc <- endoth_sfc[st_is_valid(endoth_sfc)]
  endoth_sfc <- st_cast(endoth_sfc, "POLYGON")
  endoth_inoutline <- st_intersection(endoth_sfc, outline_sfc)
  endoth_inoutline <- st_cast(endoth_sfc, "MULTIPOLYGON")
  #
  outline_poly <- lapply(outline_sfc, function(g){
    outline <- st_coordinates(g)
    ## reverse Y coords
    outline[,"Y"] <- outline[,"Y"]*-1
    return(outline)
  })
  tumor_poly <- lapply(tumor_inoutline, function(g){
    tumor <- st_coordinates(g)
    ## reverse Y coords
    tumor[,"Y"] <- tumor[,"Y"]*-1
    return(tumor)
  })
  endoth_poly <- lapply(endoth_inoutline, function(g){
    endoth <- st_coordinates(g)
    ## reverse Y coords
    endoth[,"Y"] <- endoth[,"Y"]*-1
    return(endoth)
  })
  
  plot_networks <- network_segs %>% filter(roi == roi_) %>% 
    left_join(segs, by = c("roi", "segment1" = "segment")) %>%
    left_join(segs, by = c("roi", "segment2" = "segment"), suffix = c("_1", "_2"))
  plot_network_cells <- network_segs %>% filter(roi == roi_) %>% 
    gather(i, segment, segment1, segment2) %>% select(roi, segment, network) %>% distinct %>%
    inner_join(segs)
  
  tumor_polygg <- lapply(tumor_poly, function(poly) {
    geom_polygon(data = poly, aes(x = X, y = Y, fill = '_A'), colour = NA)
  })
  endoth_polygg <- lapply(endoth_poly, function(poly) {
    geom_polygon(data = poly, aes(x = X, y = Y, fill = '_B'), colour = NA)
  })
  outline_polygg <- lapply(outline_poly, function(poly) {
    geom_polygon(data = poly, aes(x = X, y = Y), fill = NA, colour = 'black', linewidth = 0.1)
  })
  p <- ggplot(data = NULL) +
    tumor_polygg +
    endoth_polygg +
    outline_polygg +
    geom_segment(data = plot_networks, aes(x = x_um_1, xend = x_um_2, y = -y_um_1, yend = -y_um_2),
                 size = 0.1) +
    geom_point(data = plot_network_cells, aes(x = x_um, y = -y_um, fill = network),
               size = 0.5, shape = 21, colour = 'black', stroke = 0.05) +
    scale_fill_manual(values = c('grey90','#efbcb8',color_sorted),
                      breaks = c('_A','_B',net_order),
                      labels = c("Tumor region", "Endothelial region",
                                 net_order),
                      drop = FALSE) +
    scale_y_continuous(breaks = seq(0, -800, -200),
                       labels = c("0", "", "400", "", "800"),
                       sec.axis = sec_axis(~ ., breaks = seq(0,-800,-200), labels = c("0", "", "400", "", "800"))) +
    scale_x_continuous(breaks = seq(0, 800, 200), labels = paste0(seq(0, 800, 200), " um"),
                       sec.axis = sec_axis(~ ., breaks = seq(0, 800, 200))) +
    coord_fixed(xlim = c(0, 800), ylim = c(-800, 0), expand = 0) +
    labs(fill = "Annotation", x = roi_name) +
    guides(fill = guide_legend(override.aes = list(size = 3))) +
    theme_bw() +
    theme(text = element_text(size = 5), title = element_text(size = 5),
          axis.title.y = element_blank(),
          axis.title.x = element_text(hjust = 0),
          axis.text.y.left = element_text(hjust = 1.3), axis.text.x = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks = element_line(size = 0.05),
          panel.border = element_rect(linewidth = 0.1, colour = 'black'),
          panel.grid = element_blank(), plot.margin = unit(c(1,1,1,1), 'pt'),
          legend.key.size = unit(c(0.5,0.5), 'line'), legend.title.position = 'top', legend.justification = 'left')
  return(p)
})

get_legend(print_legend[[1]]) %>% as_ggplot()
# ggsave(paste0("../figures/_figure3.legend.pdf"), height = 1.5, width = 1.25, unit = 'in', family="ArialMT")

# export <- network3 %>% filter(roi == "PT0683_090919-BX-1_ROI3") %>% select(segment1, type1)
# export <- network3 %>% filter(roi == "PT0683_090919-BX-1_ROI3") %>% select(segment2, type2) %>%
#   full_join(export, by = c("segment2" = "segment1", "type2" = "type1")) %>% distinct
# write.table(export, file = "~/Downloads/export.csv", sep = ',', quote = F, col.names = F, row.names = F)


# network_bar <- network_segs %>% 
#   gather(i, segment, segment1, segment2) %>% select(roi, segment, network) %>% distinct %>%
#   full_join(segs)
net_order <- paste0("Network ", c(3, 5, 4, 1, 2))
network_bar <- segs %>% mutate(type = ifelse(is.na(leiden), type, leiden)) %>%
  left_join(pts) %>% 
  mutate(istumor = ifelse(type == "Tumor", TRUE, NA)) %>%
  mutate(network = ifelse(is.na(network), "Other cell populations", as.character(network)))
network_bar_ <- network_bar %>% group_by(roi, subject, timepoint, group, network, istumor) %>% count %>% # include/exclude roi
  group_by(roi, subject, timepoint) %>% mutate(frac = n/sum(n))
network_bar_ <- network_bar_ %>% filter(subject == "PT0550") #filter(grepl("Ipi", group)) # Switch for combo
network_bar_$network <- factor(network_bar_$network, levels = c("Other cell populations", rev(net_order)))
#
order_sub <- network_bar_ %>% filter(istumor) %>% arrange(group, timepoint, -frac)
# network_bar_$subject <- factor(network_bar_$subject, levels = unique(order_sub$subject))
network_bar_$roi <- factor(network_bar_$roi, levels = unique(order_sub$roi))

ggplot(network_bar_, # Filter to on-therapy
       aes(x = roi, y = n, fill = network, alpha = ifelse(is.na(istumor), FALSE, TRUE), colour = istumor)) +
  geom_col(position = 'fill', size = 0) +
  geom_col(position = 'fill', size = 0.5, fill = NA) +
  scale_fill_manual(values = c('grey90',color_sorted[5:1])) +
  scale_colour_manual(values = c('black'), na.value = NA, na.translate = FALSE,
                      labels = "SOX10+ melanoma") +
  scale_alpha_manual(values = c(1, 1), guide = FALSE) +
  scale_x_discrete(breaks = levels(network_bar_$roi), labels = gsub(".+(ROI\\d)", "\\1", levels(network_bar_$roi))) +
  facet_grid(. ~ paste0("PT0550- ", timepoint), scales = 'free_x', space = 'free_x') +
  guides(fill = guide_legend(reverse = TRUE, override.aes = list(colour = 'white', linewidth = 0.1)),
         colour = guide_legend(override.aes = list(fill = NA))) +
  labs(colour = "Tumor cells", fill = "Louvain Network",
       y = "Proportion of all cells", x = "ROI") +
  theme_bw() +
  theme(text = element_text(size = 5), title = element_text(size = 5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.border = element_rect(linewidth = 0.1, colour = 'black'),
        strip.background = element_blank(),
        legend.key.size = unit(c(0.5,0.5),'line'),
    axis.ticks = element_line(size = 0.1), plot.margin = unit(c(1,4,1,1), 'pt'), legend.justification = 'top', legend.margin = margin(0,0,0,0),
    legend.spacing = unit(0,'pt'), legend.position = 'none', strip.clip = "off", panel.spacing = unit(1,'pt')
    )
# ggsave("../figures/figure5d.pdf", height = 1.5, width = 1.5, family="ArialMT")
# ggsave("../figures/figure5f.pdf", height = 1.5, width = 2.6, family="ArialMT")

# ae <- list(theme(text = element_text(size = 5), title = element_text(size = 7),
#                  axis.ticks = element_line(size = 0.1), panel.border = element_rect(linewidth = 0.1, colour = 'black'),
#                  legend.key.height = unit(0.25, 'line'), legend.key.width = unit(0.25, 'line')
# ))
# Try densities instead
network_segs <- lapply(unique(pts$network), function(netname){
  pops <- pts %>% filter(network == netname)
  net_ <- net %>% filter(type1 %in% pops$type & type2 %in% pops$type)
  pops_ <- segs %>% 
    mutate(type = ifelse(is.na(leiden), type, leiden)) %>%
    filter(type %in% pops$type)
  pops_ <- split(pops_, pops_$roi)
  nets_ <- lapply(pops_, function(p_){
    t_ <- split(p_, p_$type)
    pp_ <- lapply(t_, function(t_i) { ppp(t_i$x_um, t_i$y_um, window = owin(c(0, 800), c(0, 800))) })
    apply(net_, 1, function(i) {
      if ( ! (i['type1'] %in% names(pp_) & i['type2'] %in% names(pp_) ) ) {
        return(NULL)
      } else if ( nrow(t_[[i['type1']]]) > 0 & nrow(t_[[i['type2']]])>0 ) {
        cd <- crossdist(pp_[[i['type1']]], pp_[[i['type2']]])
        ind <- which(cd < (ceiling(as.numeric(i['radius_mean']) / 5) * 5), arr.ind = TRUE)
        t1 <- t_[[i['type1']]] %>% mutate(i1 = 1:n()) %>% select(segment, i1)
        t2 <- t_[[i['type2']]] %>% mutate(i2 = 1:n()) %>% select(segment, i2)
        t12 <- tibble(type1 = i['type1'], type2 = i['type2'], i1 = ind[,'row'], i2 = ind[,'col']) %>%
          left_join(t1, by = c("i1")) %>% left_join(t2, by = c("i2"), suffix = c("1", "2"))
      } else {
        return(NULL)
      }
    }) %>% enframe %>% unnest(value)
  })
  network <- nets_ %>% enframe(name = "roi") %>% unnest(value) %>% mutate(value = NULL)
  return(network)
})
network_segs <- network_segs %>% enframe(name = "network") %>% mutate(network = unique(pts$network)) %>%
  unnest(value)

# Cell type densities
roi_dens <- segs %>% group_by(subject, roi, timepoint, group) %>% summarise(area_px = sum(area))
roi_dens <- roi_dens %>% ungroup %>% mutate(area_mm2 = area_px * (.800/2048)**2)
roi_1s <- roi_dens %>% group_by(subject, timepoint, group) %>% count(name = "nrois") %>% filter(nrois==1)
roi_dens <- rbind(roi_1s, roi_1s) %>% mutate(roi_i = 1:n()) %>% full_join(roi_dens)
roiXnet <- network_segs %>% select(network) %>% distinct() %>% cross_join(roi_dens)
# tissue_dens <- segs %>% group_by(subject, timepoint, group) %>% summarise(area_px = sum(area))
# tissue_dens <- tissue_dens %>% ungroup %>% mutate(area_mm2 = area_px * (.800/2048)**2)
# Cells actually tied to these networks
network_dens <- network_segs %>%
  gather(i, segment, segment1, segment2) %>% select(roi, segment, network) %>% distinct %>%
  full_join(segs)
network_dens_ <- network_dens %>% filter(type != "Tumor") %>%
  group_by(roi, subject, timepoint, group, network) %>% count %>%
  full_join(roiXnet) %>% rowwise %>% mutate(dens = ifelse(is.na(n), 0, n/area_mm2))
network_dens_ <- network_dens_ %>% filter(grepl("Comb", group), timepoint == "ontx c2")
network_dens_ <- network_dens_ %>% filter(!is.na(network))
network_dens_$subject <- factor(network_dens_$subject, levels = rev(unique(order_sub$subject)))

ggplot(network_dens_, 
       aes(y = subject, x = dens, fill = network, size = roi_i)) +
  ggridges::geom_density_ridges(alpha = 0.5, scale = TRUE, rel_min_height = 0.01) +
  # ggridges::geom_density_ridges(stat = "identity", alpha = 0.5) +
  geom_point(shape = 21) +
  scale_fill_manual(values = color_sorted) +
  facet_wrap( ~ group, scales = 'free') +
  # scale_x_log10() +
  theme_bw()





#
networkassoc_dens <- segs %>% mutate(type = ifelse(is.na(leiden), type, leiden)) %>%
  left_join(pts) %>% 
  mutate(istumor = ifelse(type == "Tumor", TRUE, NA)) %>%
  mutate(network = ifelse(is.na(network), "Other cell populations", as.character(network)))
network_bar_ <- network_bar %>% group_by(subject, timepoint, group, network, istumor) %>% count %>%
  group_by(subject, timepoint) %>% mutate(frac = n/sum(n))


##### Comb. PD #####
# Null distribution model to identify samples with plasma cells
plasma <- segs %>% filter(type == "Plasma cell")
plasma_ <- split(plasma, plasma$roi)
cd4 <- segs %>% filter(type_group == "CD4")
cd4_ <- split(cd4, cd4$roi)
cd4_4 <- segs %>% filter(type_group == "CD4", leiden == "CD4_4")
cd4_4_ <- split(cd4_4, cd4_4$roi)
tumor <- segs %>% filter(type == "Tumor")
tumor_ <- split(tumor, tumor$roi)


permutation_model <- function(id, plasma_cells, cd4_cells, 
                              distance_threshold = 50, n_permutations = 5000){
  # Create spatial point patterns
  plasma_pp <- ppp(plasma_cells$x_um, plasma_cells$y_um, 
                   window = owin(c(0, 800), c(0, 800)))
  cd4_pp <- ppp(cd4_cells$x_um, cd4_cells$y_um,
                window = owin(c(0, 800), c(0, 800)))
  
  if (nrow(plasma_cells)>1 & nrow(cd4_cells)>1) {
    # Calculate observed metric
    observed_distances <- crossdist(plasma_pp, cd4_pp)
    observed_neighbors <- observed_distances <= distance_threshold & observed_distances>0# Number of pairs within distance
    # observed_metric <- mean((rowSums(observed_neighbors) > 0)) # Average proportion of plasma cells with a CD4
    # observed_metric <- mean(rowSums(observed_neighbors))

    # Initialize storage for null model metrics
    null_metrics <- numeric(n_permutations)

    # Null model: shuffle locations and calculate metric
    for (i in 1:n_permutations) {
      # Shuffle plasma cell coordinates
      shuffled_plasma <- plasma_cells
      shuffled_plasma$x_um <- runif(nrow(plasma_cells), min = 0, max = 800)
      shuffled_plasma$y_um <- runif(nrow(plasma_cells), min = 0, max = 800)

      # Shuffle CD4 T cell coordinates
      shuffled_cd4 <- cd4_cells
      shuffled_cd4$x_um <- runif(nrow(cd4_cells), min = 0, max = 800)
      shuffled_cd4$y_um <- runif(nrow(cd4_cells), min = 0, max = 800)

      # Create shuffled point patterns
      shuffled_plasma_pp <- ppp(shuffled_plasma$x_um, shuffled_plasma$y_um, window = owin(c(0, 800), c(0, 800)))
      shuffled_cd4_pp <- ppp(shuffled_cd4$x_um, shuffled_cd4$y_um, window = owin(c(0, 800), c(0, 800)))

      # Calculate distances and neighbors
      shuffled_distances <- crossdist(shuffled_plasma_pp, shuffled_cd4_pp)
      shuffled_neighbors <- shuffled_distances <= distance_threshold

      # Metric: Proportion of plasma cells with neighboring CD4 T cells
      null_metrics[i] <- mean(rowSums(shuffled_neighbors) > 0)
      # null_metrics[i] <- mean(rowSums(shuffled_neighbors))
    }

    # Compare observed metric to null distribution
    p_value <- mean(null_metrics >= observed_metric)

  } else {
    observed_metric <- NA
    null_metrics <- NA
    p_value <- NA
  }
  results <- data.frame(
    id = id,
    n_1 = nrow(plasma_cells),
    n_2 = nrow(cd4_cells),
    observed = observed_metric,
    mean_null = mean(null_metrics),
    p_value = p_value
  )
  return(results)
}

roi_dens <- segs %>% group_by(subject, roi, timepoint, group) %>% summarise(area_px = sum(area), total_cells = n())
roi_dens <- roi_dens %>% ungroup %>% mutate(area_mm2 = area_px * (.800/2048)**2)

# plasm_cd4_4 <- lapply(intersect(names(plasma_), names(cd4_4_)),
#                       function(roi) {
#   permutation_model(roi, plasma_[[roi]], cd4_4_[[roi]])
# })
# plasm_cd4_4 <- plasm_cd4_4 %>% enframe %>% unnest(value)
# plasm_cd4_4 <- plasm_cd4_4 %>% full_join(rois, by = c("id" = "roi")) %>% 
#   left_join(roi_dens, by = c("id" = "roi")) %>%
#   left_join(mibi_) %>%
#   mutate(n_1 = ifelse(is.na(n_1), 0, n_1), n_2 = ifelse(is.na(n_2), 0, n_2),
#          p_value = ifelse(is.na(p_value), 1, p_value),
#          n1_dens = n_1/area_mm2, n1_prop = n_1/total_cells,
#          n2_dens = n_2/area_mm2, n2_prop = n_2/total_cells)

# plasm_cd4 <- lapply(intersect(names(plasma_), names(cd4_)),
#                       function(roi) {
#                         permutation_model(roi, plasma_[[roi]], cd4_[[roi]])
#                       })
# plasm_cd4 <- plasm_cd4 %>% enframe %>% unnest(value)
# plasm_cd4 <- plasm_cd4 %>% full_join(rois, by = c("id" = "roi")) %>% 
#   left_join(roi_dens, by = c("id" = "roi")) %>%
#   left_join(mibi_) %>%
#   mutate(n_1 = ifelse(is.na(n_1), 0, n_1), n_2 = ifelse(is.na(n_2), 0, n_2),
#          p_value = ifelse(is.na(p_value), 1, p_value),
#          n1_dens = n_1/area_mm2, n1_prop = n_1/total_cells,
#          n2_dens = n_2/area_mm2, n2_prop = n_2/total_cells)

# plasm_tumor <- lapply(intersect(names(plasma_), names(tumor_)),
#                     function(roi) {
#                       permutation_model(roi, plasma_[[roi]], tumor_[[roi]])
#                     })
# plasm_tumor <- plasm_tumor %>% enframe %>% unnest(value)
# plasm_tumor <- plasm_tumor %>% full_join(rois, by = c("id" = "roi")) %>% 
#   left_join(roi_dens, by = c("id" = "roi")) %>%
#   left_join(mibi_) %>%
#   mutate(n_1 = ifelse(is.na(n_1), 0, n_1), n_2 = ifelse(is.na(n_2), 0, n_2),
#          p_value = ifelse(is.na(p_value), 1, p_value),
#          n1_dens = n_1/area_mm2, n1_prop = n_1/total_cells,
#          n2_dens = n_2/area_mm2, n2_prop = n_2/total_cells)

# tumor_plasm <- lapply(intersect(names(plasma_), names(tumor_)),
#                       function(roi) {
#                         permutation_model(roi, tumor_[[roi]], plasma_[[roi]])
#                       })
# tumor_plasm <- tumor_plasm %>% enframe %>% unnest(value)
# tumor_plasm <- tumor_plasm %>% full_join(rois, by = c("id" = "roi")) %>% 
#   left_join(roi_dens, by = c("id" = "roi")) %>%
#   left_join(mibi_) %>%
#   mutate(n_1 = ifelse(is.na(n_1), 0, n_1), n_2 = ifelse(is.na(n_2), 0, n_2),
#          p_value = ifelse(is.na(p_value), 1, p_value),
#          n1_dens = n_1/area_mm2, n1_prop = n_1/total_cells,
#          n2_dens = n_2/area_mm2, n2_prop = n_2/total_cells)

# plasm_plasm <- lapply(intersect(names(plasma_), names(plasma_)),
#                       function(roi) {
#                         permutation_model(roi, plasma_[[roi]], plasma_[[roi]])
#                       })
# plasm_plasm <- plasm_plasm %>% enframe %>% unnest(value)
# plasm_plasm <- plasm_plasm %>% full_join(rois, by = c("id" = "roi")) %>% 
#   left_join(roi_dens, by = c("id" = "roi")) %>%
#   left_join(mibi_) %>%
#   mutate(n_1 = ifelse(is.na(n_1), 0, n_1), n_2 = ifelse(is.na(n_2), 0, n_2),
#          p_value = ifelse(is.na(p_value), 1, p_value),
#          n1_dens = n_1/area_mm2, n1_prop = n_1/total_cells,
#          n2_dens = n_2/area_mm2, n2_prop = n_2/total_cells)

# permutation_results <- lapply(list(plasm_cd4_4, plasm_cd4, plasm_plasm), function(i) { mutate(i, name = NULL) }) %>%
#   enframe %>% mutate(name = "permutation model, 50um, 5000 permutations") %>%
#   mutate(type_1 = c("Plasma cell", "Plasma cell", "Plasma cell"), type_2 = c("CD4_4", "CD4 T cell", "Plasma cell")) %>%
#   select(name, type_1, type_2, value) %>% unnest(value)
# write.table(permutation_results, file = "table_s9_pre.tsv", sep = '\t', quote = FALSE, row.names = FALSE)
# permutation_results <- lapply(list(plasm_plasm), function(i) { mutate(i, name = NULL) }) %>%
#   enframe %>% mutate(name = "permutation model, 50um, 5000 permutations") %>%
#   mutate(type_1 = c("Plasma cell"), type_2 = c("Plasma cell")) %>%
#   select(name, type_1, type_2, value) %>% unnest(value)
# write.table(permutation_results, file = "table_s9.tsv", sep = '\t', quote = FALSE, row.names = FALSE)
permutation_results <- fread("table_s9.tsv")
permutation

# Plot for Fig 3F
# Plasma cell densities
roi_dens <- segs %>% group_by(subject, roi, timepoint, group) %>% summarise(area_px = sum(area), total_cells = n())
roi_dens <- roi_dens %>% ungroup %>% mutate(area_mm2 = area_px * (.800/2048)**2)
tissue_dens <- segs %>% group_by(subject, timepoint, group) %>% summarise(area_px = sum(area), total_cells = n())
tissue_dens <- tissue_dens %>% ungroup %>% mutate(area_mm2 = area_px * (.800/2048)**2)
roi_plasma <- segs %>% filter(type == "Plasma cell") %>% group_by(roi) %>% count(name = "n_plasma") %>%
  full_join(rois) %>% mutate(n_plasma = ifelse(is.na(n_plasma), 0, n_plasma)) %>%
  full_join(roi_dens) %>% mutate(plasma_dens = n_plasma/area_mm2) %>%
  mutate(y = ifelse(group == "Comb. PD", subject, as.character(group)),
         facet = case_when(group == "Comb. CR/PR" ~ "2", group == "Comb. PD" ~ "1", TRUE ~ "3"))
tissue_plasma <- roi_plasma %>% group_by(facet, y, subject, timepoint, group) %>% summarise(n_plasma = sum(n_plasma)) %>%
  full_join(tissue_dens) %>% mutate(plasma_dens = n_plasma/area_mm2)
order_pts <- tissue_plasma %>% filter(group == "Comb. PD") %>% arrange(timepoint, -plasma_dens)
tissue_plasma$y <- factor(tissue_plasma$y, levels = c(unique(order_pts$y), 'Comb. CR/PR','Ipi PD'))
roi_plasma$y <- factor(roi_plasma$y, levels = c(unique(order_pts$y), 'Comb. CR/PR','Ipi PD'))
missing_ <- missing %>% 
  mutate(y = ifelse(group == "Comb. PD", subject, as.character(group)),
         facet = case_when(group == "Comb. CR/PR" ~ "2", group == "Comb. PD" ~ "1", TRUE ~ "3")) %>%
  filter(facet == "1")
missing_$y <- factor(missing_$y, levels = levels(tissue_plasma$y))

sig_results <- permutation_results %>% filter(n_1>4 & (p_value < 0.05 | n1_dens > 100)) %>%
  group_by(id) %>% summarise(sig = TRUE) %>%
  inner_join(roi_plasma, by = c("id" = "roi"))
nonsig_results <- roi_plasma %>% anti_join(sig_results, by = c("roi" = "id"))
ggplot(data = NULL, aes(x = plasma_dens+1, y = as.numeric(y))) +
  geom_vline(xintercept = 100, linetype = 'dashed', colour = 'grey50', size = 0.25) +
  geom_segment(data = tissue_plasma, aes(xend = plasma_dens+1, y = as.numeric(y)+0.5, yend = as.numeric(y)-0.5),
               colour = NA) +
  geom_point(data = nonsig_results, aes(fill = timepoint), shape = 21, size = 1, stroke = 0.1) +
  geom_point(data = sig_results, aes(fill = timepoint, colour = sig), shape = 21, size = 1, stroke = 0.25) +
  geom_segment(data = tissue_plasma, aes(xend = plasma_dens+1, y = as.numeric(y)+0.25, yend = as.numeric(y)-0.25),
               colour = 'black', size = 0.5) +
  geom_point(data = missing_, aes(shape = "Missing"), x = 0, size = 1, stroke = 0.5) +
  scale_x_continuous(trans = 'log10') +
  scale_y_reverse(breaks = 1:length(levels(tissue_plasma$y)), labels = levels(tissue_plasma$y)) +
  coord_cartesian(expand = 0, xlim = c(0.5, 5000)) +
  scale_shape_manual(labels = "Paired sample N/A", values = 4) +
  scale_fill_manual(values = as.character(timepoint_colors)) +
  scale_colour_manual(values = "red3", na.value = 'black', labels = c("Monte Carlo p<0.05")) +
  facet_grid(facet ~ timepoint, scales = 'free_y', space = 'free_y') +
  labs(y = "Subject (Comb. PD)", x = "Plasma cell density (cells per mm2)", shape = "",
       fill = "Timepoint", color = "Permutation test") +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y = element_blank(),
        text = element_text(size = 5), title = element_text(size = 5),
        legend.key.size = unit(c(0.5,0.5),'line'), legend.position = 'right',
        legend.justification = c(1,1), legend.margin = margin(2,2,2,2,'pt'),legend.spacing = unit(0,'pt'),
        plot.margin = unit(c(2,2,2,2), 'pt'),
        axis.title.y = element_text(hjust = 0.65))
ggsave("../figures/figure4f.pdf", height = 2, width = 5.2, family="ArialMT")

# perm_ <- permutation_results %>% 
#   mutate(y = ifelse(group == "Comb. PD", subject, as.character(group)),
#          facet = case_when(group == "Comb. CR/PR" ~ "2", group == "Comb. PD" ~ "1", TRUE ~ "3"))
# perm_$y <- factor(perm_$y, levels = levels(tissue_plasma$y))
# ggplot(data = NULL, aes(x = p_value, y = as.numeric(y))) +
#   geom_segment(data = tissue_plasma, aes(x = 0, xend = 0, y = as.numeric(y)+0.5, yend = as.numeric(y)-0.5),
#                colour = NA) +
#   geom_point(data = perm_, aes(fill = timepoint), shape = 21, stroke = 0.1) +
#   scale_y_reverse(breaks = 1:length(levels(tissue_plasma$y)), labels = levels(tissue_plasma$y)) +
#   scale_x_reverse() +
#   facet_grid(facet ~ type_2, scales = 'free_y', space = 'free_y') +
#   labs(y = "Subject") +
#   theme_bw() +
#   theme(strip.background = element_blank(), strip.text.y = element_blank())
# tmp <- fread("table_s9_pre.tsv")
# tmp <- tmp %>% filter(type_2 == "CD4 T cell")
prep_perm <- permutation_results %>% select(type_2) %>% distinct
prep_perm <- prep_perm %>% cross_join(mibi_)
permutation_results %>% filter(group == "Comb. PD", timepoint == "baseline") %>% summary
# plasma cell density has an average of 109 across all samples
# plasma cell number median is 5 for these samples
summ_perm <- permutation_results %>% filter(n_1>4 & (p_value < 0.05 | n1_dens > 100)) %>% 
  group_by(type_2, subject, timepoint, group) %>%
  summarise(n_rois = n(), sig = TRUE) %>%
  full_join(prep_perm) %>% mutate(n_rois = ifelse(is.na(n_rois), 0, n_rois), sig = ifelse(is.na(sig), FALSE, sig))
ggplot(summ_perm, aes(x = group, fill = interaction(group, sig), alpha = sig)) +
  geom_bar(position = 'fill') +
  facet_wrap(~ timepoint)
summ_perm_ <- split(summ_perm, summ_perm$timepoint, drop = TRUE)
fisher <- lapply(summ_perm_, function(df){
  df <- filter(df, grepl("Comb", group))
  contingency <- table(df$group, df$sig)
  fisher <- fisher.test(contingency)
  signif(fisher$p.value, 2)
}) %>% enframe(name = "timepoint") %>% unnest(value)
fisher$timepoint <- factor(fisher$timepoint)

ggplot(summ_perm, aes(x = group, fill = group, alpha = sig)) +
  geom_bar(position = 'fill') +
  geom_segment(data = NULL, inherit.aes = FALSE, x = 1, xend = 2, y = 1.02, yend = 1.02, linewidth = 0.25) +
  geom_segment(data = NULL, inherit.aes = FALSE, x = 2, xend = 2, y = 1.01, yend = 1.02, linewidth = 0.25) +
  geom_segment(data = NULL, inherit.aes = FALSE, x = 1, xend = 1, y = 1.01, yend = 1.02, linewidth = 0.25) +
  geom_text(data = fisher, aes(x = 1.5, y = 1.03, label = value), vjust = 0, inherit.aes = FALSE, size = 5*0.36) +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_fill_manual(values = group_colors) +
  coord_cartesian(ylim = c(0, 1.035)) +
  facet_wrap(~ timepoint) +
  labs(y = "% of patients with dense CD138+ plasma cells") +
  theme_bw() +
  theme(strip.background = element_blank(), legend.position = 'none', axis.title.x = element_blank(),
        text = element_text(size = 5), title = element_text(size = 5),
        legend.key.size = unit(c(0.5,0.5),'line'),
        legend.justification = c(1,0),
        plot.margin = unit(c(2,2,2,2), 'pt'))
ggsave("../figures/figure4g.pdf", height = 2, width = 2.5, family="ArialMT")

# Plasma cells observed on-therapy in responding biopsies



# Ipi highlight: PT0550


##### LATER #####
# CURRENTLY NOT USING
## Networks in Comb. PD
net <- report_spot$`Ontx, Combo: CR/PR v PD` %>%
  filter(coef<0, pval<0.05) %>% group_by(type1, type2) %>% 
  summarise(n = n(), across(c('radius','coef'), list(median = median, mean = mean, sd = sd)))

# Convert the data frame to a graph object
g <- igraph::graph_from_data_frame(
  d = net,
  directed = FALSE
)

# Assign radius_mean as the edge weights
igraph::E(g)$weight <- net$radius_mean
# Perform clustering using the Louvain method
clusters <- igraph::cluster_louvain(g)
igraph::V(g)$cluster <- igraph::membership(clusters)

# Sort vertices by cluster manually
sorted_vertices <- order(igraph::V(g)$cluster)
g <- igraph::permute(g, sorted_vertices)

# Create custome graph plot
pts <- data.frame(v = igraph::V(g), cluster = igraph::V(g)$cluster) %>% rownames_to_column("type") %>% 
  mutate(network = paste0("Network ", cluster+5)) %>% select(type, network)
pts$network <- factor(pts$network)#, levels = paste0("Network ", c(2,1,4,5,3)))
#
edges <- net %>% ungroup %>% mutate(edge = 1:n()) %>%
  left_join(pts, by = c("type1" = "type")) %>% left_join(pts, by = c("type2" = "type"), suffix = c("1","2"))
horizontal_edges <- edges %>% filter(network1 == network2) %>%
  gather(i, type, type1, type2) %>% left_join(pts) %>% group_by(network, type) %>%
  count(name = "n_horizontal")
vertical_edges <- edges %>% filter(network1 != network2) %>%
  gather(i, type, type1, type2) %>% left_join(pts) %>% group_by(network, type) %>%
  count(name = "n_vertical")
weight_edges <- horizontal_edges %>% full_join(vertical_edges) %>% 
  mutate(n_vertical = ifelse(is.na(n_vertical), 0, n_vertical))
#
org_edges <- weight_edges %>% arrange(network, -n_horizontal, -n_vertical) %>%
  group_by(network) %>% mutate(x = (1:n()-1)*(-1)^(1:n())) %>%
  arrange(network, x) %>% mutate(x = 1:n()) %>% select(network, type, x)
push <- edges %>% filter(network1 != network2) %>% 
  left_join(org_edges, by = c("type1" = "type", "network1" = "network")) %>% 
  left_join(org_edges, by = c("type2" = "type", "network2" = "network"), suffix = c("1", "2")) %>%
  mutate(xdiff = x2-x1) %>% select(matches('type'), xdiff) %>% gather(i, type, -xdiff) %>%
  left_join(pts) %>% group_by(type, network) %>% summarise(push = mean(xdiff)) %>%
  group_by(network) %>% summarise(push = mean(push))
org_edges <- org_edges %>% left_join(push) %>% 
  mutate(push = ifelse(is.na(push), 0, push)) %>%
  mutate(x = x + abs(push))
#
plot_edges <- edges %>% select(matches('type'), n, matches('network')) %>%
  left_join(org_edges, by = c("type1" = "type", "network1" = "network")) %>% 
  left_join(org_edges, by = c("type2" = "type", "network2" = "network"), suffix = c("1", "2"))
plot_edges$network1 <- factor(plot_edges$network1, levels = levels(org_edges$network))
plot_edges$network2 <- factor(plot_edges$network2, levels = levels(org_edges$network))
#
get_labs <- rgg %>%
  filter(scores>5, pvals_adj<0.05) %>%
  group_by(type_group, leiden) %>% slice_head(n = 3) %>%
  group_by(type_group, leiden) %>% summarise(n = n(), label = paste0(names, collapse = ", ")) %>%
  mutate(label = gsub("HLA class 1 A, B, and C, Na-K-ATPase", "HLA I", label)) %>%
  mutate(label = gsub("HLA DDQDR", "HLA II", label)) %>%
  mutate(label = gsub("fracneighbors_", "", label)) %>%
  mutate(label = gsub("PD-1", "PD1", label))
#
org_edges_ <- org_edges %>% left_join(get_labs, by = c("type" = "leiden"))

rects <- org_edges %>% mutate(y = as.numeric(network)) %>%
  group_by(network, y) %>% summarise(minx = min(x), maxx = max(x))

org_edges_$type_group <- factor(org_edges_$type_group, levels = c('CD8','CD4','Myeloid'))

## GRAPH DIAGRAM THIS ONE
ggplot(data = NULL, aes(y = as.numeric(network))) +
  # rectangles
  geom_text(data = rects, aes(x = minx-0.55, y = y, label = paste0("Louvain ", network)),
            hjust = 1.05, size = 5*0.36) +
  geom_rect(data = rects, aes(xmin = minx - 0.5, xmax = maxx + 0.5, ymin = y-0.375, ymax = y+0.375, fill = factor(network, levels = rev(levels(rects$network)))),
            inherit.aes = FALSE, colour = NA, alpha = 0.2) +
  # edges
  geom_curve(data = filter(plot_edges, network1 != network2),
             aes(y = as.numeric(network1), yend = as.numeric(network2), x = x1, xend = x2), 
             curvature = 0.1, alpha = 0.2, size = 0.1) +
  geom_curve(data = filter(plot_edges, network1 == network2),
             aes(y = as.numeric(network1), yend = as.numeric(network2), x = x1, xend = x2), 
             curvature = 0.2, alpha = 0.6, size = 0.1) +
  # vertices
  geom_point(data = org_edges_, aes(x = x, fill = factor(network, levels = rev(levels(rects$network))), 
                                    shape = type_group), size = 5, stroke = 0.05) +
  # cell labels
  geom_label(data = org_edges, aes(x = x, y = as.numeric(network)+0.1, label = gsub("Treg", "", type)),
             vjust = 0, size = 5*0.36, label.size = 0.05) +
  # descriptors
  geom_label(data = org_edges_, aes(x = x, y = as.numeric(network)-0.1, label = str_wrap(label, 16)), # 
             vjust = 1, size = 5*0.36, fill = 'white', alpha = 1, label.size = 0.05) +
  scale_y_continuous(breaks = 1:length(levels(pts$network)), labels = levels(pts$network)) +
  scale_x_continuous(limits = c(min(rects$minx-1.5), max(rects$maxx+0.5))) +
  scale_shape_manual(values = c(21, 22, 23), na.value = 25, labels = c("CD8 T cell subsets", "CD4 T cell subsets", "Myeloid subsets", "Other cell types")) +
  scale_fill_manual(values = color_sorted[c(6:8,11:12)]) +
  labs(shape = "Cell populations", fill = "Network\n(Louvain cluster)") +
  # guides(fill = guide_legend(override.aes = list(alpha = 1, shape = NA))) +
  guides(fill = 'none',
         shape = guide_legend(override.aes = list(size = 4, stroke = 0.25))) +
  theme_void() +
  theme(text = element_text(size = 5), title = element_text(size = 5),
        legend.key.size = unit(c(0.5,0.5),'line'),
        legend.position = c(1.11, 0.1), legend.direction = 'vertical', 
        legend.justification = c(1,0),
        plot.margin = unit(c(2,2,2,2), 'pt'))
ggsave("~/manuscript/Biopsy analysis of trial S1616/figures/figure4e.pdf", height = 3, width = 4.5, unit = 'in', family="ArialMT")

# Print individual regions
network_segs <- lapply(unique(pts$network), function(netname){
  pops <- pts %>% filter(network == netname)
  net_ <- net %>% filter(type1 %in% pops$type & type2 %in% pops$type)
  pops_ <- segs %>% 
    mutate(type = ifelse(is.na(leiden), type, leiden)) %>%
    filter(type %in% pops$type)
  #
  pops_ <- split(pops_, pops_$roi)
  nets_ <- lapply(pops_, function(p_){
    t_ <- split(p_, p_$type)
    pp_ <- lapply(t_, function(t_i) { ppp(t_i$x_um, t_i$y_um, window = owin(c(0, 800), c(0, 800))) })
    apply(net_, 1, function(i) {
      if ( ! (i['type1'] %in% names(pp_) & i['type2'] %in% names(pp_) ) ) {
        return(NULL)
      } else if ( nrow(t_[[i['type1']]]) > 0 & nrow(t_[[i['type2']]])>0 ) {
        cd <- crossdist(pp_[[i['type1']]], pp_[[i['type2']]])
        ind <- which(cd < (ceiling(as.numeric(i['radius_mean']) / 5) * 5), arr.ind = TRUE)
        t1 <- t_[[i['type1']]] %>% mutate(i1 = 1:n()) %>% select(segment, i1)
        t2 <- t_[[i['type2']]] %>% mutate(i2 = 1:n()) %>% select(segment, i2)
        t12 <- tibble(type1 = i['type1'], type2 = i['type2'], i1 = ind[,'row'], i2 = ind[,'col']) %>%
          left_join(t1, by = c("i1")) %>% left_join(t2, by = c("i2"), suffix = c("1", "2"))
      } else {
        return(NULL)
      }
    }) %>% enframe %>% unnest(value)
  })
  network <- nets_ %>% enframe(name = "roi") %>% unnest(value) %>% mutate(value = NULL)
  return(network)
})
network_segs <- network_segs %>% enframe(name = "network") %>% mutate(network = unique(pts$network)) %>%
  unnest(value)
network_segs %>% group_by(network, roi) %>% count %>%
  left_join(rois) %>% left_join(mibi_) %>% arrange(-n) %>%
  filter(! subject %in% c('PT0735','PT0621'), group == "Comb. PD")

net_order <- paste0("Network ", c(9, 8, 7, 6, 10))

lapply(grep("PT0542", unique(network_segs$roi), value = TRUE), function(roi_){
  print(roi_)
  roi_name = str_replace(roi_, "_", "-")
  outline_path <- "~/Google Drive/My Drive/Morrison/MORRISON ptII- Cohort 3/morrison-cohort-3-mibi/structures/"
  outlines <- list.files(path = outline_path, pattern = "dsDNA.+geojson")
  tumors <- list.files(path = outline_path, pattern = "SOX10.+geojson")
  endothelials <- list.files(path = outline_path, pattern = "CD31.+geojson")
  roi_outline <- paste0(outline_path, grep(str_replace(roi_, "_", "-"), outlines, value = TRUE))
  roi_tumor <- paste0(outline_path, grep(str_replace(roi_, "_", "-"), tumors, value = TRUE))
  roi_endoth <- paste0(outline_path, grep(str_replace(roi_, "_", "-"), endothelials, value = TRUE))
  outline_geom <- st_read(roi_outline) %>% st_geometry()
  tumor_geom <- st_read(roi_tumor) %>% st_geometry()
  endoth_geom <- st_read(roi_endoth) %>% st_geometry()
  #
  outline_sfc <- lapply(outline_geom, st_polygon) %>% st_as_sfc()
  tumor_sfc <- lapply(tumor_geom, st_polygon) %>% st_as_sfc()
  # tumor_sfc <- tumor_sfc[st_is_valid(tumor_sfc)]
  tumor_inoutline <- st_intersection(tumor_sfc, outline_sfc)#tumor_sfc#
  tumor_inoutline <- st_cast(tumor_inoutline, "MULTIPOLYGON")
  endoth_sfc <- lapply(endoth_geom, st_polygon) %>% st_as_sfc()
  endoth_sfc <- endoth_sfc[st_is_valid(endoth_sfc)]
  endoth_sfc <- st_cast(endoth_sfc, "POLYGON")
  endoth_inoutline <- st_intersection(endoth_sfc, outline_sfc)
  endoth_inoutline <- st_cast(endoth_sfc, "MULTIPOLYGON")
  #
  outline_poly <- lapply(outline_sfc, function(g){
    outline <- st_coordinates(g)
    ## reverse Y coords
    outline[,"Y"] <- outline[,"Y"]*-1
    return(outline)
  })
  tumor_poly <- lapply(tumor_inoutline, function(g){
    tumor <- st_coordinates(g)
    ## reverse Y coords
    tumor[,"Y"] <- tumor[,"Y"]*-1
    return(tumor)
  })
  endoth_poly <- lapply(endoth_inoutline, function(g){
    endoth <- st_coordinates(g)
    ## reverse Y coords
    endoth[,"Y"] <- endoth[,"Y"]*-1
    return(endoth)
  })
  
  plot_networks <- network_segs %>% filter(roi == roi_) %>% 
    left_join(segs, by = c("roi", "segment1" = "segment")) %>%
    left_join(segs, by = c("roi", "segment2" = "segment"), suffix = c("_1", "_2"))
  plot_network_cells <- network_segs %>% filter(roi == roi_) %>% 
    gather(i, segment, segment1, segment2) %>% select(roi, segment, network) %>% distinct %>%
    inner_join(segs)
  
  tumor_polygg <- lapply(tumor_poly, function(poly) {
    geom_polygon(data = poly, aes(x = X, y = Y, fill = '_A'), colour = NA)
  })
  endoth_polygg <- lapply(endoth_poly, function(poly) {
    geom_polygon(data = poly, aes(x = X, y = Y, fill = '_B'), colour = NA)
  })
  outline_polygg <- lapply(outline_poly, function(poly) {
    geom_polygon(data = poly, aes(x = X, y = Y), fill = NA, colour = 'black', linewidth = 0.1)
  })
  p <- ggplot(data = NULL) +
    tumor_polygg +
    endoth_polygg +
    outline_polygg +
    geom_segment(data = plot_networks, aes(x = x_um_1, xend = x_um_2, y = -y_um_1, yend = -y_um_2),
                 size = 0.1) +
    geom_point(data = plot_network_cells, aes(x = x_um, y = -y_um, fill = network),
               size = 0.5, shape = 21, colour = 'black', stroke = 0.05) +
    # scale_fill_manual(values = color_sorted[c(6:8,11:12)]) +
    scale_fill_manual(values = c('grey90','#efbcb8',color_sorted[c(6:8,11:12)]),
                      breaks = c('_A','_B',net_order),
                      labels = c("Tumor region", "Endothelial region",
                                 net_order),
                      drop = FALSE) +
    scale_y_continuous(breaks = seq(0, -800, -200),
                       labels = c("0", "", "400", "", "800"),
                       sec.axis = sec_axis(~ ., breaks = seq(0,-800,-200), labels = c("0", "", "400", "", "800"))) +
    scale_x_continuous(breaks = seq(0, 800, 200), labels = paste0(seq(0, 800, 200), " um"),
                       sec.axis = sec_axis(~ ., breaks = seq(0, 800, 200))) +
    coord_fixed(xlim = c(0, 800), ylim = c(-800, 0), expand = 0) +
    labs(fill = "Annotation", x = roi_name) +
    guides(fill = guide_legend(nrow = 2)) +
    theme_bw() +
    theme(text = element_text(size = 5), title = element_text(size = 5),
          axis.title.y = element_blank(),
          axis.title.x = element_text(hjust = 0),
          axis.text.y.left = element_text(hjust = 1.3), axis.text.x = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks = element_line(size = 0.05),
          panel.border = element_rect(linewidth = 0.1, colour = 'black'),
          panel.grid = element_blank(), plot.margin = unit(c(1,1,1,1), 'pt'),
          legend.position = 'none')
  p
  ggsave(paste0("../figures/4f/nets.",roi_,".pdf"), height = 1.5, width = 1.5, unit = 'in', family="ArialMT")
  return(p)
})





##### END #####


## ARXIV
#
# outline_geom <- lapply(outline_sfc, function(g){
#   # outline <- st_polygon(g)
#   ## get the geom coordinates as data.frame
#   outline <- st_coordinates(g)
#   ## reverse Y coords
#   outline[,"Y"] <- outline[,"Y"]*-1
#   ## re-create geom
#   outline <- st_as_sfc(st_as_text(st_linestring(outline)))
# })
# tumor_geom <- lapply(tumor_inoutline, function(g){
#   # tumor <- st_polygon(g)
#   ## get the geom coordinates as data.frame
#   tumor <- st_coordinates(g)
#   ## reverse Y coords
#   tumor[,"Y"] <- tumor[,"Y"]*-1
#   ## re-create geom
#   tumor <- st_as_sfc(st_as_text(st_linestring(tumor)))
# })
# # outline_geom <- outline_geom %>% st_as_sfc
# outline_layers <- lapply(outline_geom, function(i) { 
#   geom_sf(data = st_polygon(i), aes(fill = "Tissue"), linewidth = 0.1)
# })
# tumor_layers <- lapply(tumor_geom, function(i) { 
#   geom_sf(data = st_polygon(i), aes(fill = "Tumor"), linewidth = 0.1)
# })


