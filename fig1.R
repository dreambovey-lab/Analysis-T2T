library(stats)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)
library(openxlsx)

### data preprocessing
sample_lvls <- c("D0", "D1-1", "D1-2", "D1-2-1", "D1-2-2", "D1-2-3", "D1-3", "D1-2(10)", "D1-1(100)", "D2-1", "D2-2", "D2-2-1", "D2-2-2", "D2-2-3", "D2-3", "D2-2(10)", "D2-1(100)", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8")
ILMN = read.table("ILMN_performance_stats.tsv", sep = "\t", h = T) %>% 
	mutate(Lab = gsub("_.*", "", Lab)) %>% 
	dplyr::filter(!grepl("Undet", Prefix) & !grepl("NC|_R|S3", Sample)) %>% 
	mutate(
		Lab = factor(Lab), 
		Sample = factor(gsub("_D", "", Sample), levels = sample_lvls)) %>% 
	arrange(Lab, Sample)
MGI = read.table("MGI_performance_stats.tsv", sep = "\t", h = T) %>% 
	mutate(Lab = gsub("_.*", "", Lab)) %>% 
	dplyr::filter(!grepl("NC|_R", Sample)) %>% 
	mutate(
		Lab = factor(Lab), 
		Sample = factor(gsub("_D", "", Sample), levels = sample_lvls)) %>% 
	arrange(Lab, Sample)

sequencers = c("ILMN", "MGI")
sequencer_cols <- c(ILMN = "#D55E00", MGI = "#E7B800")

### consensus plot
p_consensus <- wrap_plots(lapply(sequencers, function(sequencer){
	assign(
		paste0("p_cons_", sequencer), 
		ggplot(eval(sym(sequencer)) %>% 
				select(Lab, Sample, Host_consensus, Non_host_consensus, Multi_labels, Host_reads_ground_truth) %>% 
				unite("Sample", c(Lab, Sample), sep = ": ") %>% 
				mutate(Sample = factor(Sample, levels = unique(Sample))) %>% 
				pivot_longer(c("Non_host_consensus", "Multi_labels", "Host_consensus"), "Label") %>% 
				mutate(Label = factor(Label, levels = c("Non_host_consensus", "Multi_labels", "Host_consensus"))), 
			aes(x = fct_rev(Sample), y = value, fill = Label)) + 
		geom_bar(stat = "identity", position = "stack") + 
		scale_y_continuous(expand = c(0, 0)) + 
		theme_classic() + 
		theme(
			axis.text.x = element_text(angle = 60, hjust = 1), 
			axis.text.y = element_blank(), 
			axis.ticks.y = element_blank(), 
			legend.key.size=unit(0.4, "cm")) + 
		scale_fill_manual(
			labels = c("Non-host consensus", "Multi-labels", "Host consensus"), 
			values = c("#A6CEE3", "grey50", "#FB9A99")) + 
		ggtitle(sequencer) + 
		labs(x = "", y = "", fill = "") + 
		coord_flip(), 
	envir = globalenv())
}), nrow = 1, guides = "collect") & theme(legend.position = 'bottom')
df_avg_host_ratio <- rbind(ILMN %>% select(Sequencer, Host_reads_perc), MGI %>% select(Sequencer, Host_reads_perc)) %>% mutate(Sequencer = factor(Sequencer, levels = c("ILMN", "MGI")))

### avg host ratio plot
p_avg_host_ratio <- ggplot(df_avg_host_ratio, aes(x = Sequencer, y = Host_reads_perc, fill = Sequencer)) + geom_half_violin(adjust = 0.5, trim = T, color = "black", side = "r", scale = "width", alpha = 0.9) + geom_half_point(aes(color = Sequencer), side = "l", transformation = position_jitter(width = 0.05), alpha = 0.9) + scale_fill_manual(values = sequencer_cols) + scale_color_manual(values = sequencer_cols) + labs(x = "", y = "Host reads (%)", fill = "", color = "") + guides(fill = guide_legend(nrow = 1, byrow = TRUE)) + theme_classic() + theme(legend.position = "none")


### performance stats and plot
performance_measures <- c("sens", "spec", "MCC")
ref_genome <- c("T2T", "hg38", "YH")

if(!file.exists("stats_tests/compare_algorithm/")) {dir.create("stats_tests/compare_algorithm/")}
p_values <- unlist(lapply(performance_measures, function(measure) {
	my_wb <- createWorkbook()
	p <- lapply(sequencers, function(sequencer){
		my_wb2 <- createWorkbook()
		df <- eval(sym(sequencer)) %>% select(Lab, Sample, ends_with(measure)) %>% unite("Sample", c(Sample, Lab)) %>% mutate(Sample = factor(Sample, levels = unique(Sample))) %>% pivot_longer(ends_with(measure), "Method") %>% separate(Method, c("Genome", "Algorithm", "Measure"), sep = "__") %>% mutate(Genome = factor(Genome, levels = c("hg38", "T2T", "YH")), Algorithm = factor(gsub("_", " ", Algorithm), levels = c("bowtie2 vf", "bowtie2 vs", "bwa mem", "minimap2", "winnowmap", "STAT", "kraken2")))
		df %>% group_by(Genome) %>% group_walk(~{for(type in c("greater", "less")){sheetname = paste(.y$Genome, type, sep = "_"); addWorksheet(my_wb2, sheetName = sheetname); writeData(my_wb2, sheetname, pairwise.wilcox.test(.x$value, .x$Algorithm, paired = T, p.adj = "bonf", alternative = type)$p.value, rowNames = T)}})
		if(measure == "MCC"){
			df_perf_plot <- df %>% dplyr::filter(Measure == "MCC") %>% pivot_wider(names_from = "Algorithm", values_from = "value") %>% select(-Sample, -Measure)
			num_comparison <- (ncol(df_perf_plot) - 1)*(ncol(df_perf_plot) - 2)/2
			legendBarDiag <- function(data, mapping, ...) {ggally_barDiag(data = data, mapping = mapping, ...) + annotate("text", label = "hg38\nn=67", x = 1, y = Inf, vjust = 2.5) + annotate("text", label = "T2T\nn=67", x = 2, y = Inf, vjust = 2.5) + annotate("text", label = "YH\nn=67", x = 3, y = Inf, vjust = 2.5)}
			ablinePoints <- function(data, mapping, ...) {ggally_points(data = data, mapping = mapping, ...) + geom_abline(slope = 1, intercept = 0)}
			medianPairDiff <- function(x, y) {p1 <- wilcox.test(x, y, alternative = "greater", paired = T)$p.value;p2 <- wilcox.test(x, y, alternative = "less", paired = T)$p.value;return(paste(round(median(x - y), 4), case_when(p1 < 0.001 | p2 < 0.001 ~ "***", p1 < 0.01 | p2 < 0.01 ~ "**", p1 < 0.05 | p2 < 0.05 ~ "*", TRUE ~ "")))}
			lessBreaksHist <- function(data, mapping, ...) {ggally_facethist(data = data, mapping = mapping, ...) + scale_y_continuous(n.breaks = 4)}
			mat_plot <- ggpairs(df_perf_plot , aes(color = Genome), diag = list(discrete = wrap(legendBarDiag, alpha = 0.9), continuous = wrap("densityDiag", alpha = 0.8)), upper = list(continuous = wrap("statistic", text_fn = medianPairDiff, title = "p-Wil", display_grid = F), combo = wrap("box_no_facet", alpha = 0.9)), lower = list(continuous = wrap(ablinePoints, alpha = 0.9), combo = wrap(lessBreaksHist, alpha = 0.9, color = "black", bins = 20))) + theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
			for(row in seq_len(mat_plot$nrow))
				for(col in seq_len(mat_plot$ncol))
					mat_plot[row, col] <- mat_plot[row, col] +
						scale_fill_brewer(palette = "Set2") +
						scale_color_brewer(palette = "Set2")
			assign(paste0("p_", sequencer), mat_plot, envir = globalenv())
		}
		for(type in c("greater", "less")){
			sheetname = paste(sequencer, type, sep = "_")
			addWorksheet(my_wb, sheetName = sheetname)
			writeData(my_wb, sheetname, pairwise.wilcox.test(df$value, interaction(df$Genome, df$Algorithm), paired = T, p.adj = "bonf", alternative = type)$p.value, rowNames = T)
		}
		saveWorkbook(my_wb2, file = paste0("stats_tests/compare_algorithm/", measure, "_", sequencer, ".xlsx"), overwrite = T)
		kw = kruskal.test(df$value, interaction(df$Genome, df$Algorithm))
		return(c(measure, sequencer, kw$p.value))
	})
	saveWorkbook(my_wb, file = paste0("stats_tests/genome_cross_algorithm_", measure, ".xlsx"), overwrite = T)
	return(p)
}), recursive = F)
write.table(data.frame(do.call(rbind, p_values)), "stats_tests/Omnibus_test_Kruskal.txt", sep = "\t", col.names = F, row.names = F)

g1 <- plot_grid(
	plot_grid(
		plot_grid(p_consensus, p_avg_host_ratio, 
		nrow = 1, rel_widths = c(2.5, 1), 
		labels = c("A", "B"), 
		align = "h", 
		axis = "b" 
		), 
		ncol = 1,
		rel_heights = c(1.2, 1),
		labels = c("", "C")
	), 
	plot_grid(ggmatrix_gtable(p_ILMN), ggmatrix_gtable(p_MGI), ncol = 1, rel_heights = c(1, 1), labels = c("D", "E")), 
	ncol = 2, rel_widths = c(1, 1.2)
)

ggsave("Figure1.pdf", g1, width = 18, height = 18)
ggsave("Figure1.tiff", g1, dpi = 320, w = 400, h = 400, units = "mm", compression = "lzw")
ggsave("Figure1.png", g1, w = 400, h = 400, units = "mm")