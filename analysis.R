### data preprocessing
sample_lvls <- c("D0", "D1-1", "D1-2", "D1-2-1", "D1-2-2", "D1-2-3", "D1-3", "D1-2(10)", "D1-1(100)", "D2-1", "D2-2", "D2-2-1", "D2-2-2", "D2-2-3", "D2-3", "D2-2(10)", "D2-1(100)", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8")
ONT <- read.table("ONT_performance_stats.tsv", sep = "\t", h = T) %>% 
	dplyr::filter(!grepl("NC", Sample)) %>% 
	mutate(
		Lab = factor(gsub("POOLING[12]_", "", Lab)), 
		Sample = factor(
			gsub("_D", "", Sample), 
			levels = sample_lvls)) %>% 
	arrange(Lab, Sample)
ILMN <- read.table("ILMN_performance_stats.tsv", sep = "\t", h = T) %>% 
	mutate(Lab = gsub("_.*", "", Lab)) %>% 
	dplyr::filter(!grepl("Undet", Prefix) & !grepl("NC|_R|S3", Sample)) %>% 
	mutate(
		Lab = factor(Lab), 
		Sample = factor(
			gsub("_D", "", Sample), 
			levels = sample_lvls)) %>% 
	arrange(Lab, Sample)
MGI <- read.table("MGI_performance_stats.tsv", sep = "\t", h = T) %>% 
	mutate(Lab = gsub("_.*", "", Lab)) %>% 
	dplyr::filter(!grepl("NC|_R", Sample)) %>% 
	mutate(
		Lab = factor(Lab), 
		Sample = factor(
			gsub("_D", "", Sample), 
			levels = sample_lvls)) %>% 
	arrange(Lab, Sample)

performance_measures = c("sens", "spec", "MCC")
sequencers = c("ILMN", "MGI", "ONT")


### consensus
library(patchwork)
library(gghalves)
sequencer_cols <- c(ILMN = "#D55E00", MGI = "#E7B800", ONT = "#0072B2")
p_consensus <- wrap_plots(unlist(lapply(sequencers, function(sequencer){
	p_label <- ggplot(eval(sym(sequencer)) %>% mutate(Lab = gsub("(Lab.*|[n]?HD)", "Ref_\\1", Lab)) %>% count(Lab) %>% arrange(desc(Lab)) %>% mutate(end = cumsum(n), start = end - n + 1)) + geom_segment(aes(x = 1, xend = 1, y = start - 0.25, yend = end + 0.25), linewidth = 2) + geom_text(aes(x = 0, y = (start + end)/2, label = Lab), angle = 90) + scale_x_continuous(limits = c(-0.5, 1), expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0.5)) + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
	p_barplot <- ggplot(
			eval(sym(sequencer)) %>% 
				select(Lab, Sample, Host_consensus, Non_host_consensus, Multi_labels, Host_reads_ground_truth) %>% 
				mutate(Sample = paste(Lab, Sample, sep = ": ")) %>% 
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
			legend.key.size=unit(0.4, "cm"), 
			plot.margin = margin(0, 0, 0, 0)) + 
		scale_fill_manual(
			labels = c("Non-host consensus", "Multi-labels", "Host consensus"), 
			values = c("#A6CEE3", "grey50", "#FB9A99")) + 
		ggtitle(sequencer) + 
		labs(x = "", y = "", fill = "") + 
		coord_flip()
	return(list(p_label, p_barplot))
}), recursive = F), nrow = 1, guides = "collect", widths = rep(c(1, 8), length(sequencers))) & theme(legend.position = 'bottom')

df_avg_host_ratio <- rbind(
		ILMN %>% select(Sequencer, Host_reads_perc), 
		MGI %>% select(Sequencer, Host_reads_perc), 
		ONT %>% select(Sequencer, Host_reads_perc)) %>% 
	mutate(Sequencer = factor(Sequencer, levels = c("ILMN", "MGI", "ONT")))
p_avg_host_ratio <- ggplot(
		df_avg_host_ratio, 
		aes(
			x = Sequencer, 
			y = Host_reads_perc, 
			fill = Sequencer)) + 
	geom_half_violin(
		adjust = 0.5, 
		trim = T, 
		color = "black", 
		side = "r", 
		scale = "width", 
		alpha = 0.9) + 
	geom_half_point(
		aes(color = Sequencer), 
		side = "l", 
		transformation = position_jitter(width = 0.05), 
		alpha = 0.9) + 
	scale_fill_manual(values = sequencer_cols) + 
	scale_color_manual(values = sequencer_cols) + 
	labs(x = "", y = "Host reads (%)", fill = "", color = "") + 
	guides(fill = guide_legend(nrow = 1, byrow = TRUE)) + 
	theme_classic() + 
	theme(legend.position = "none")


### performance
library(stats)
library(GGally)
library(cowplot)
library(patchwork)
library(openxlsx)
# performance_measures <- c("sens", "spec", "MCC")
# sequencers <- c("ILMN", "MGI", "ONT")
ref_genome <- c("T2T", "hg38", "YH")
#p_ <- GGally::print_if_interactive
#genome_cols <- c(hg38 = "#00AFBB", YH = "#E7B800", T2T = "#FC4E07")

p_values <- unlist(lapply(performance_measures, function(measure) {
	my_wb <- createWorkbook()
	p <- lapply(sequencers, function(sequencer){
		my_wb2 <- createWorkbook()
		df <- eval(sym(sequencer)) %>% select(Lab, Sample, ends_with(measure)) %>% unite("Sample", c(Sample, Lab)) %>% mutate(Sample = factor(Sample, levels = unique(Sample))) %>% pivot_longer(ends_with(measure), "Method") %>% separate(Method, c("Genome", "Algorithm", "Measure"), sep = "__") %>% mutate(Genome = factor(Genome, levels = c("hg38", "T2T", "YH")), Algorithm = factor(gsub("_", " ", Algorithm), levels = c("bowtie2 vf", "bowtie2 vs", "bwa mem", "minimap2", "winnowmap", "STAT", "kraken2")))
		df %>% group_by(Genome) %>% group_walk(~{for(type in c("greater", "less")){sheetname = paste(.y$Genome, type, sep = "_"); addWorksheet(my_wb2, sheetName = sheetname); writeData(my_wb2, sheetname, pairwise.wilcox.test(.x$value, .x$Algorithm, paired = T, p.adj = "bonf", alternative = type)$p.value, rowNames = T)}})
		df_perf_plot <- df %>% dplyr::filter(Measure == measure) %>% pivot_wider(names_from = "Algorithm", values_from = "value") %>% select(-Measure)
		num_comparison <- (ncol(df_perf_plot) - 1)*(ncol(df_perf_plot) - 2)/2
		legendBarDiag <- function(data, mapping, ...) {sample_size = data %>% pull(Sample) %>% n_distinct(); ggally_barDiag(data = data, mapping = mapping, ...) + annotate("text", label = paste0("hg38\nn=", sample_size), x = 1, y = Inf, vjust = 2.5) + annotate("text", label = paste0("T2T\nn=", sample_size), x = 2, y = Inf, vjust = 2.5) + annotate("text", label = paste0("YH\nn=", sample_size), x = 3, y = Inf, vjust = 2.5)}
		ablinePoints <- function(data, mapping, ...) {ggally_points(data = data, mapping = mapping, ...) + geom_abline(slope = 1, intercept = 0)}
		medianPairDiff <- function(x, y) {p1 <- wilcox.test(x, y, alternative = "greater", paired = T)$p.value;p2 <- wilcox.test(x, y, alternative = "less", paired = T)$p.value;return(paste(round(median(x - y), 4), case_when(p1 < 0.001 | p2 < 0.001 ~ "***", p1 < 0.01 | p2 < 0.01 ~ "**", p1 < 0.05 | p2 < 0.05 ~ "*", TRUE ~ "")))}
		lessBreaksHist <- function(data, mapping, ...) {ggally_facethist(data = data, mapping = mapping, ...) + scale_y_continuous(n.breaks = 4)}
		#medianPairDiff <- function(x, y) {p1 <- p.adjust(wilcox.test(x, y, alternative = "greater", paired = T)$p.value, method = "fdr", n = 10);p2 <- p.adjust(wilcox.test(x, y, alternative = "less", paired = T)$p.value, method = "fdr", n = 10);return(paste(round(median(x - y), 4), case_when(p1 < 0.001 | p2 < 0.001 ~ "***", p1 < 0.01 | p2 < 0.01 ~ "**", p1 < 0.05 | p2 < 0.05 ~ "*", TRUE ~ "")))}
		mat_plot <- ggpairs(df_perf_plot, mapping = aes(color = Genome), columns = 2:ncol(df_perf_plot), diag = list(discrete = wrap(legendBarDiag, alpha = 0.9), continuous = wrap("densityDiag", alpha = 0.8)), upper = list(continuous = wrap("statistic", text_fn = medianPairDiff, title = "p-Wil", display_grid = F), combo = wrap("box_no_facet", alpha = 0.9)), lower = list(continuous = wrap(ablinePoints, alpha = 0.9), combo = wrap(lessBreaksHist, alpha = 0.9, color = "black", bins = 20))) + theme_bw() + theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
		for(row in seq_len(mat_plot$nrow))
			for(col in seq_len(mat_plot$ncol))
				mat_plot[row, col] <- mat_plot[row, col] +
					scale_fill_brewer(palette = "Set2") +
					scale_color_brewer(palette = "Set2")
		#ggsave(paste0("matrix_plot/", sequencer, ".pdf"), mat_plot, h = 12, w = 16)
		assign(paste("p", measure, sequencer, sep = "_"), mat_plot, envir = globalenv())
		for(type in c("greater", "less")){
			sheetname = paste(sequencer, type, sep = "_")
			addWorksheet(my_wb, sheetName = sheetname)
			writeData(my_wb, sheetname, pairwise.wilcox.test(df$value, interaction(df$Genome, df$Algorithm), paired = T, p.adj = "bonf", alternative = type)$p.value, rowNames = T)
		}
		#saveWorkbook(my_wb2, file = paste0("stats_tests/compare_algorithm/", measure, "_", sequencer, ".xlsx"), overwrite = T)
		kw = kruskal.test(df$value, interaction(df$Genome, df$Algorithm))
		return(c(measure, sequencer, kw$p.value))
	})
	#saveWorkbook(my_wb, file = paste0("stats_tests/genome_cross_algorithm_", measure, ".xlsx"), overwrite = T)
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
		ggmatrix_gtable(p_MCC_ONT),
		ncol = 1,
		rel_heights = c(1.2, 1),
		labels = c("", "E")
	), 
	plot_grid(ggmatrix_gtable(p_MCC_ILMN), ggmatrix_gtable(p_MCC_MGI), ncol = 1, rel_heights = c(1, 1), labels = c("C", "D")), 
	ncol = 2, rel_widths = c(1, 1.2)
)

ggsave("Figure1.pdf", g1, width = 18, height = 18)
ggsave("Figure1.tiff", g1, dpi = 320, w = 400, h = 400, units = "mm", compression = "lzw")
ggsave("Figure1.png", g1, w = 400, h = 400, units = "mm")

gs_sens <- plot_grid(ggmatrix_gtable(p_sens_ILMN), ggmatrix_gtable(p_sens_MGI), ggmatrix_gtable(p_sens_ONT), ncol = 1, rel_heights = c(1.2, 1.2, 1), labels = c("I", "M", "O"))
ggsave("Sens.pdf", gs_sens, w = 9, h = 24)

gs_spec <- plot_grid(ggmatrix_gtable(p_spec_ILMN), ggmatrix_gtable(p_spec_MGI), ggmatrix_gtable(p_spec_ONT), ncol = 1, rel_heights = c(1.2, 1.2, 1), labels = c("I", "M", "O"))
ggsave("Spec.pdf", gs_spec, w = 9, h = 24)

rm(p_MCC_ILMN, p_MCC_MGI, p_MCC_ONT, p_sens_ILMN, p_sens_MGI, p_sens_ONT, p_spec_ILMN, p_spec_MGI, p_spec_ONT, g1, gs_sens, gs_spec)

# compare genomes
library(rstatix)
cmp_gnm_df <- bind_rows(
	ILMN %>% 
		select(c(Sequencer, Sample, Lab, ends_with(performance_measures))) %>% 
		pivot_longer(-c(Sequencer, Sample, Lab)) %>% 
		separate(name, c("Genome", "Algorithm", "Metric"), sep = "__"), 
	MGI %>% 
		select(c(Sequencer, Sample, Lab, ends_with(performance_measures))) %>% 
		pivot_longer(-c(Sequencer, Sample, Lab)) %>% 
		separate(name, c("Genome", "Algorithm", "Metric"), sep = "__"), 
	ONT %>% 
		select(c(Sequencer, Sample, Lab, ends_with(performance_measures))) %>% 
		pivot_longer(-c(Sequencer, Sample, Lab)) %>% 
		separate(name, c("Genome", "Algorithm", "Metric"), sep = "__")
	) %>% 
	dplyr::filter(Algorithm == "bwa_mem" | (Algorithm == "kraken2" & Sequencer == "ONT"))
cmp_gnm_res <- cmp_gnm_df %>% 
	group_by(Sequencer, Metric) %>% 
	group_modify(
		~.x %>% 
		pairwise_wilcox_test(
			value ~ Genome, 
			paired = T, 
			ref.group = "T2T", 
			alternative = "greater", 
			detailed = T)) %>% 
		select(
			Sequencer, 
			Metric, 
			estimate, 
			group1, 
			group2, 
			p, 
			p.adj, 
			p.adj.signif) %>% 
		mutate(Metric = factor(Metric, levels = c("sens", "spec", "MCC"))) %>% 
		arrange(Sequencer, Metric, group2)
write.table(cmp_gnm_res, "compare_genome_stats.tsv", sep = "\t", quote = F, row.names = F)

### QC
library(patchwork)

QC_ONT <- read.table("ONT_QC_summary.tsv", sep = "\t", h = T)
QC_SR <- read.table("SR_QC_summary.tsv", sep = "\t", h = T)
total_info <- bind_rows(
	QC_ONT %>% 
	select(Prefix, Sequencer, Run, Raw_Reads, Clean_Reads, Mean_Lengths, Q10, Mean_GC) %>%
	mutate(Q20 = NA, Dup = NA, Mean_GC = Mean_GC/100) %>% 
	rename(
		"Lab" = "Run", 
		"Mean length" = "Mean_Lengths", 
		"Mean GC content" = "Mean_GC"
	), 
	QC_SR %>% 
	select(Prefix, Sequencer, Lab, Raw_Reads, Clean_Reads, Length, Q20, GC, Dup) %>% 
	mutate(Q10 = NA) %>% rename("Mean length" = "Length", "Mean GC content" = "GC")) %>% mutate(
		Lab = gsub("clinical", "Clinical", 
				gsub("POOLING[12]", "Ref", 
					gsub("(Lab[0-9]+).*", "Ref_\\1", Lab)))
	) %>% 
	rename("Clean reads" = "Clean_Reads", "Raw reads" = "Raw_Reads") %>% 
	dplyr::filter(!grepl("NC|Undet|R[1-8]-1|RT[1-8]", Prefix) & `Clean reads` > 1000)
p_QC_ILMN <- ggplot(total_info %>% 
	dplyr::filter(Sequencer == "ILMN") %>% 
	select(-Prefix, -Q10) %>% 
	pivot_longer(-c(1:2), names_to = "key", values_to = "value") %>% 
	mutate(
		key = factor(key, 
			levels = c("Raw reads", "Clean reads", "Mean length", 
						"Q10", "Q20", "Mean GC content", "Dup"))
	), 
	aes(x = Lab, y = value, color = Lab)) + 
	geom_boxplot() + 
	labs(x = "", y = "", color = "Samples") + 
	facet_wrap(.~key, scales = "free") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 60, hjust = 1))
p_QC_MGI <- ggplot(total_info %>% 
	dplyr::filter(Sequencer == "MGI") %>% 
	select(-Prefix, -Q10) %>% 
	pivot_longer(-c(1:2), names_to = "key", values_to = "value") %>% 
	mutate(
		key = factor(key, 
			levels = c("Raw reads", "Clean reads", "Mean length", 
						"Q10", "Q20", "Mean GC content", "Dup"))
	), 
	aes(x = Lab, y = value, color = Lab)) + 
	geom_boxplot() + 
	labs(x = "", y = "", color = "Samples") + 
	facet_wrap(.~key, scales = "free") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 60, hjust = 1))
p_QC_ONT <- ggplot(total_info %>% 
	dplyr::filter(Sequencer == "ONT") %>% 
	select(-Prefix, -Q20, -Dup) %>% 
	pivot_longer(-c(1:2), names_to = "key", values_to = "value") %>% 
	mutate(
		key = factor(key, 
			levels = c("Raw reads", "Clean reads", "Mean length", 
						"Q10", "Q20", "Mean GC content", "Dup"))
	), 
	aes(x = Lab, y = value, color = Lab)) + 
	geom_boxplot() + 
	labs(x = "", y = "", color = "Samples") + 
	facet_wrap(.~key, scales = "free") + 
	theme_bw() + 
	theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave("Figure S1 Sequencing statistics.pdf", 
	p_QC_ILMN + p_QC_MGI + p_QC_ONT + 
	plot_annotation(tag_levels = "A") + 
	plot_layout(nrow = 1),
	h = 9, 
	w = 24))



### per read result compilation example
library(ggforce)
example_lab = "clinical"
example_prefix = "SD-003-8"
# example_prefix = "S4-1-5"
# example_lab = "HD"
host_taxid_hierarchy <- c(7711, 89593, 7742, 7776, 117570, 117571, 8287, 1338369, 32523, 32524, 40674, 32525, 9347, 1437010, 314146, 9443, 376913, 314293, 9526, 314295, 9604, 207598, 9605, 9606)
sample_lvls <- c("D0", "D1-1", "D1-2", "D1-2-1", "D1-2-2", "D1-2-3", "D1-3", "D1-2(10)", "D1-1(100)", "D2-1", "D2-2", "D2-2-1", "D2-2-2", "D2-2-3", "D2-3", "D2-2(10)", "D2-1(100)", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8")

example_stats <- read.table("ONT_performance_stats.tsv", sep = "\t", h = T) %>% 
	dplyr::filter(!grepl("NC", Sample)) %>% 
	mutate(
		Lab = factor(gsub("POOLING[12]_", "", Lab)), 
		Sample = factor(
		gsub("_D", "", Sample), 
		levels = sample_lvls)) %>% 
	dplyr::filter(Lab == example_lab, Prefix == example_prefix)
multi_label <- read.table(paste0("method_example/", example_prefix, "_multi_labels.tsv"), sep = "\t", h = T) %>% 
	select(-c("h", "nh", "kraken_h", "kraken_nh")) %>% 
	mutate(across(ends_with("kraken2"), ~if_else(. %in% host_taxid_hierarchy, 1, 0))) %>% 
	left_join(read.table(paste0("method_example/", example_prefix, "_resolved.tsv"), sep = "\t", h = T), by = "ReadID") %>% 
	mutate(Loosened_host_consensus = as.numeric(is.na(Resolved))) %>% 
	replace_na(list(Resolved = 1))
df_example <- multi_label %>% 
	mutate(o = NA) %>% 
	arrange(Resolved, Loosened_host_consensus) %>% 
	mutate(ReadID = 1:n()) %>% 
	rename("Ground_truth" = "Resolved", "Imperfect_consensus" = "Loosened_host_consensus") %>% 
	mutate(ReadID = factor(ReadID)) %>% 
	pivot_longer(-ReadID, names_to = "Method", values_to = "Label") %>% 
	mutate(Method = factor(gsub("__", "_", Method), levels = c("hg38_minimap2", "hg38_winnowmap", "T2T_minimap2", "T2T_winnowmap", "YH_minimap2", "YH_winnowmap", "hg38_STAT", "hg38_kraken2", "T2T_kraken2", "YH_kraken2", "o", "Imperfect_consensus", "Ground_truth")), Label = as.factor(gsub("1", "Host", gsub("0", "Non-host", Label))))
p_example <- ggplot(df_example) + 
	geom_bar(aes(y = ReadID, x = Method, fill = Label), stat = "identity", position = "stack") + 
	annotate("segment", x = 10.75, xend = 11.25, y = 1600000, yend = 1600000, arrow = arrow(length = unit(0.3, "cm"), type = "closed"), arrow.fill = "black", linewidth = 2.5) + 
	scale_fill_manual(values = c("Non-host" = "#A6CEE3", "Host" = "#FB9A99"), na.value = NA) + 
	scale_y_discrete(name = "Read", expand = expansion(add = c(0, 1000000))) + 
	scale_x_discrete(name = "", expand = expansion(add = c(0.6, 1)), labels = c("hg38_minimap2", "hg38_winnowmap", "T2T_minimap2", "T2T_winnowmap", "YH_minimap2", "YH_winnowmap", "hg38_STAT", "hg38_kraken2", "T2T_kraken2", "YH_kraken2", "", "Imperfect_consensus", "Ground_truth")) + 
	theme_classic() + 
	theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1), axis.ticks.x = element_blank(), axis.line = element_line(arrow = arrow(length = unit(0.2, "cm"))), plot.margin = margin(t = 0, r = 0, b = 0, l = 20))
library(rsvg)
rsvg_svg("method_example/T2T_logic_flow.svg", "method_example/T2T_logic_flow.cairo.svg")
library(grImport2)
p_logic_flow <- readPicture("method_example/T2T_logic_flow.cairo.svg")
ggsave("Figure S2 Example of collating host reads removal results of all methods.pdf", h = 7.5, w = 7.5)


### benchmarking
my_colors = c(TP = "#D8B365", FN = "#EBD9B2", TN = "#5AB4AC", FP = "#ACD9D5")

draw_benchmarking_plot = function(df, my_lab){
    ggplot(df %>% dplyr::filter(Lab == my_lab)) + 
		geom_bar(aes(x = fct_rev(Method), y = value, fill = Metrics, group = Metrics), stat = "identity", position = position_stack()) + 
		geom_hline(yintercept = 0) + 
		labs(x = "", y = "Reads") + 
		ggtitle(my_lab) + 
		labs(fill = "Category") + 
		scale_fill_manual(values = c(TP = "#D8B365", FN = "#EBD9B2", TN = "#5AB4AC", FP = "#ACD9D5")) + 
		scale_y_continuous(labels = function(x) abs(x), breaks = scales::extended_breaks(6)) + 
		theme(axis.text.x = element_text(angle = 60, hjust = 1), axis.ticks.y = element_blank(), panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA), legend.position = "none") + 
		coord_flip() + 
		facet_wrap(.~Sample, strip.position = "right", ncol = 1, dir = "v")
}

df_bm_ONT = ONT %>% 
	mutate(Lab = gsub("([n]?HD)", "Ref_\\1", Lab)) %>% 
	mutate(across(ends_with("TP"), function(x) x - Host_consensus), across(ends_with("TN"), function(x) x - Non_host_consensus)) %>% 
	select(c(Sample, Lab, ends_with(c("TP", "FP", "TN", "FN")))) %>% 
	pivot_longer(-c(Sample, Lab)) %>% 
	separate(name, c("Genome", "Algorithm", "Metrics"), sep = "__") %>% 
	unite("Method", c("Genome", "Algorithm")) %>% 
	mutate(value = if_else(grepl("TN|FP", Metrics), -value, value) , Metrics = factor(Metrics, levels = c("TN", "TP", "FN", "FP")))

g_ONT = wrap_plots(draw_benchmarking_plot(df_bm_ONT, "clinical"), draw_benchmarking_plot(df_bm_ONT, "Ref_HD"), draw_benchmarking_plot(df_bm_ONT, "Ref_nHD"), nrow = 1, guides = "collect")

df_bm_ILMN = ILMN %>% 
	mutate(Lab = gsub("(Lab.*)", "Ref_\\1", Lab)) %>% 
	mutate(across(ends_with("TP"), function(x) x - Host_consensus), across(ends_with("TN"), function(x) x - Non_host_consensus)) %>% 
	select(c(Sample, Lab, ends_with(c("TP", "FP", "TN", "FN")))) %>% 
	pivot_longer(-c(Sample, Lab)) %>% 
	separate(name, c("Genome", "Algorithm", "Metrics"), sep = "__") %>% 
	unite("Method", c("Genome", "Algorithm")) %>% 
	mutate(value = if_else(grepl("TN|FP", Metrics), -value, value) , Metrics = factor(Metrics, levels = c("TN", "TP", "FN", "FP")))

g_ILMN = wrap_plots(draw_benchmarking_plot(df_bm_ILMN, "clinical"), draw_benchmarking_plot(df_bm_ILMN, "Ref_Lab4"), draw_benchmarking_plot(df_bm_ILMN, "Ref_Lab34"), draw_benchmarking_plot(df_bm_ILMN, "Ref_Lab40"), draw_benchmarking_plot(df_bm_ILMN, "Ref_Lab49"), nrow = 1, guides = "collect")

df_bm_MGI = MGI %>% 
	mutate(Lab = gsub("(Lab.*)", "Ref_\\1", Lab)) %>% 
	mutate(across(ends_with("TP"), function(x) x - Host_consensus), across(ends_with("TN"), function(x) x - Non_host_consensus)) %>% 
	select(c(Sample, Lab, ends_with(c("TP", "FP", "TN", "FN")))) %>% 
	pivot_longer(-c(Sample, Lab)) %>% 
	separate(name, c("Genome", "Algorithm", "Metrics"), sep = "__") %>% 
	unite("Method", c("Genome", "Algorithm")) %>% 
	mutate(value = if_else(grepl("TN|FP", Metrics), -value, value) , Metrics = factor(Metrics, levels = c("TN", "TP", "FN", "FP")))

g_MGI = wrap_plots(draw_benchmarking_plot(df_bm_MGI, "clinical"), draw_benchmarking_plot(df_bm_MGI, "Ref_Lab8"), draw_benchmarking_plot(df_bm_MGI, "Ref_Lab22"), draw_benchmarking_plot(df_bm_MGI, "Ref_Lab28"), draw_benchmarking_plot(df_bm_MGI, "Ref_Lab43"), nrow = 1, guides = "collect")

legend <- get_legend(
  g_ILMN +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

ggsave("Figure S3 Benchmarking of the host reads removal methods.pdf", plot_grid(plot_grid(g_ILMN, g_MGI, g_ONT, labels = c("ILMN", "MGI", "ONT"), nrow = 1, rel_widths = c(5, 5, 3)), legend, ncol=1, rel_heights = c(1, .01)), height = 24, width = 36)


### QQ plot

p_qq <- wrap_plots(unlist(lapply(performance_measures, function(measure){
	unlist(lapply(sequencers, function(sequencer){
		p <- ggplot(eval(sym(sequencer)) %>%
				select(Lab, Sample, ends_with(measure)) %>% 
				unite("Sample", c(Sample, Lab)) %>% 
				mutate(Sample = factor(Sample, levels = unique(Sample))) %>% 
				mutate(baseline_col = .data[[paste0("T2T__kraken2__", measure)]]) %>% 
				mutate(across(ends_with(measure), ~ . - baseline_col)) %>% 
				pivot_longer(ends_with(measure), "Method") %>% 
				separate(Method, c("Genome", "Algorithm", "Measure"), sep = "__")) + 
			geom_qq(aes(sample = value)) + 
			stat_qq_line(aes(sample = value), color = "red") + 
			labs(
				y = paste(measure, "pairwise difference quantiles"), 
				x = "Normal distribution quantiles") + 
			facet_wrap(.~Genome+Algorithm)
		return(list(p))
		}), recursive = F)
}), recursive = F), nrow = 3) + plot_annotation(tag_levels = "A")
ggsave("Figure S8 QQ plots of MCC.pdf", p_qq, h = 16, w = 16)


### Unique reads

# call compile.R

#ideogram
library(RIdeogram)

chr_df <- read.table("ideogram/chr_size.txt", sep = "\t", h = T) %>% 
	mutate(Chr = gsub("Chr", "", Chr)) %>% 
	mutate(Chr = gsub("Chr", "", Chr))
uniq_df <- read.table("ideogram/uniq_region.txt", sep = "\t", h = T) %>% 
	mutate(Chr = gsub("Chr", "", Chr)) %>% 
	mutate(Chr = gsub("Chr", "", Chr))
#sat
sat_df <- read.xlsx("ideogram/Human_chr_stats.xlsx", sheet = "Satellites") %>% 
	mutate(Chr = gsub("Chr", "", Chr))
#orange-red c("#fdb57d", "#f07e40")
#cyan-blue c("#85C4C9", "#3B738F")
ideogram(chr_df, overlaid = sat_df, colorset1 = c("#85C4C9", "#3B738F"))
convertSVG("chromosome.svg", device = "pdf")
#
sequencers <- c("ILMN", "MGI", "ONT")
window_size <- 10000
cutoff_ratio <- 0.05

df <- read.table("Uniq_reads_GC_loc.tsv", sep = "\t", h = T) %>% 
	mutate(Chr = gsub("Chr", "", as.character(Chr)), Pos = as.numeric(Pos)) %>% 
	dplyr::filter(!is.na(Pos) & Chr != "") %>% 
	mutate(Segment = Pos %/% window_size)
lapply(sequencers, function(sequencer){
	seg_hg38_df <- df %>% 
		dplyr::filter(Sequencer == sequencer & hg38 == 1)  %>% 
		group_by(Chr) %>% count(Segment)
	#write.table(seg_hg38_df, paste0(sequencer, "_hg38_", window_size, ".txt"), sep = "\t", quote = F, row.names = F)
	seg_hg38_df <- seg_hg38_df %>% 
		dplyr::filter(n >= max(n)*cutoff_ratio) %>% 
		transmute(
			Chr = Chr, 
			Start = Segment * window_size + 1, 
			End = (Segment + 1) * window_size + 1, 
			Value = n)

	seg_YH_df <- df %>% 
		dplyr::filter(Sequencer == sequencer & YH == 1) %>% 
		group_by(Chr) %>% count(Segment) 
	#write.table(seg_YH_df, paste0(sequencer, "_YH_", window_size, ".txt"), sep = "\t", quote = F, row.names = F)	
	seg_YH_df <- seg_YH_df %>% 
		dplyr::filter(n >= max(n)*cutoff_ratio) %>% 
		transmute(
			Chr = Chr, 
			Start = Segment * window_size + 1, 
			End = (Segment + 1) * window_size + 1, 
			Value = n)

	ideogram(chr_df, overlaid = seg_hg38_df, label = uniq_df, label_type = "marker", colorset1 = c("#95D0A6", "#33806C"))
	convertSVG("chromosome.svg", device = "pdf")
	system(paste0("mv chromosome.pdf ", sequencer, "_hg38.pdf"))
	ideogram(chr_df, overlaid = seg_hg38_df, colorset1 = c("#95D0A6", "#33806C"))
	convertSVG("chromosome.svg", device = "pdf")
	system(paste0("mv chromosome.pdf ", sequencer, "_hg38_nolabel.pdf"))

	ideogram(chr_df, overlaid = seg_YH_df, label = uniq_df, label_type = "marker", colorset1 = c("#FDDDAE", "#EC5D55"))
	convertSVG("chromosome.svg", device = "pdf")
	system(paste0("mv chromosome.pdf ", sequencer, "_YH.pdf"))
	ideogram(chr_df, overlaid = seg_YH_df, colorset1 = c("#FDDDAE", "#EC5D55"))
	convertSVG("chromosome.svg", device = "pdf")
	system(paste0("mv chromosome.pdf ", sequencer, "_YH_nolabel.pdf"))
})	


# GC-content (Length?) : violin
GC_df <- read.table("Uniq_reads_GC.tsv", sep = "\t", h = T) %>% 
	pivot_longer(c(T2T, hg38, YH), names_to = "Genome", values_to = "Bool") %>% 
	dplyr::filter(Bool == 1) %>% 
	select(-Bool)
p_GC <- ggplot(GC_df, aes(x = Genome, y = GC, fill = Genome)) + 
	geom_violin(bw = 4) + 
	geom_boxplot(width = 0.1, coef = Inf) + 
	labs(y = "GC (%)") + 
	ggpubr::stat_compare_means(
		aes(label = ..p.signif..), 
		method = "t.test", 
		comparisons = list(
			c("T2T", "hg38"), 
			c("T2T", "YH"), 
			c("hg38", "YH"))) + 
	ggpubr::stat_compare_means(
		method = "anova", 
		label.y = 140) + 
	scale_fill_brewer(palette = "Set2") + 
	theme_classic() + 
	theme(
		axis.title.x = element_blank(), 
		axis.text.y = element_text(angle = 90, hjust = 0.5), 
		legend.position = "none") + 
	facet_wrap(.~Sequencer)


# TaxID
host_taxid_hierarchy <- c(7711, 89593, 7742, 7776, 117570, 117571, 8287, 1338369, 32523, 32524, 40674, 32525, 9347, 1437010, 314146, 9443, 376913, 314293, 9526, 314295, 9604, 207598, 9605, 9606)
library(tidyverse)
sample_df <- read.table("/nas/Analyse/wanglei/T2T/analysis/QC/samples.tsv", sep = "\t", h = T)
prefix2sample <- sample_df %>% pull(Sample)
names(prefix2sample) <- sample_df %>% pull(Prefix)
sequencers = c("ILMN", "MGI", "ONT")
files <- unlist(lapply(sequencers, function(sequencer){
	files <- Sys.glob(paste0(sequencer, "/*/*"))
	if(sequencer == "ONT") {
		files <- files[!grepl("NC", files)]
	} else if(sequencer == "ILMN") {
		files <- files[!grepl("NC|S3|R[1-8]-1|Undet", files)]
	} else {
		files <- files[!grepl("NC|RT", files)]
	}
	return(files)
}))
df <- purrr::reduce(
	lapply(files, function(x)
		{
			sequencer = unlist(str_split(x, "/"))[1]
			lab = unlist(str_split(x, "/"))[2]
			prefix = str_replace(unlist(str_split(x, "/"))[3], "_genomic_stats.tsv", "")
			read.table(x, sep = "\t", h = T) %>% dplyr::filter(T2T == 1 | hg38 == 1 | YH == 1) %>% select(starts_with(c("T2T", "hg38", "YH"))) %>% mutate(Sequencer = sequencer, Lab = lab, Sample = prefix2sample[prefix])}
		), bind_rows
) %>% relocate(Sequencer, Lab, Sample)
write.table(df, "TaxID/Per_genome_FN_TaxID.tsv", sep = "\t", quote = F, row.names = F)
# Sequencer       Lab     Sample  T2T     T2T_kraken2_TaxID       hg38    hg38_kraken2_TaxID      YH      YH_kraken2_TaxID
tax_df <- bind_rows(
	df %>% dplyr::filter(T2T == 1) %>% group_by(Sequencer) %>% count(T2T_kraken2_TaxID) %>% transmute(Sequencer = Sequencer, Genome = "T2T", TaxID = T2T_kraken2_TaxID, n = n),
	df %>% dplyr::filter(hg38 == 1) %>% group_by(Sequencer) %>% count(hg38_kraken2_TaxID) %>% transmute(Sequencer = Sequencer, Genome = "hg38", TaxID = hg38_kraken2_TaxID, n = n),
	df %>% dplyr::filter(YH == 1) %>% group_by(Sequencer) %>% count(YH_kraken2_TaxID) %>% transmute(Sequencer = Sequencer, Genome = "YH", TaxID = YH_kraken2_TaxID, n = n)
)
write.table(tax_df, "TaxID_summary.tsv", sep = "\t", quote = F, row.names = F)
tax_df <- read.table("TaxID_summary.tsv", sep = "\t", h = T) %>% mutate(TaxID = as.character(TaxID))
tax_top10_df <- tax_df %>% group_by(Sequencer, Genome) %>% arrange(desc(n)) %>% mutate(rank = 1:n()) %>% mutate(TaxID = if_else(rank <= 10, TaxID, "other")) %>% group_by(Sequencer, Genome, TaxID) %>% summarize_all(sum) %>% mutate(rank = if_else(rank <= 10, rank, as.integer(11))) %>% unite("Value", c(TaxID, n)) %>% pivot_wider(names_from = Genome, values_from = Value) %>% separate("hg38", c("hg38_TaxID", "hg38_count")) %>% separate("T2T", c("T2T_TaxID", "T2T_count")) %>% separate("YH_TaxID", c("YH", "YH_count")) %>% group_by(Sequencer) %>% arrange(rank, .by_group = T)
write.table(tax_top10_df, "TaxID_top10.tsv", sep = "\t", quote = F, row.names = F)

df <- read.table("Uniq_reads_TaxID.tsv", sep = "\t", h = T)
tax_df <- df %>% 
	pivot_longer(
		ends_with("TaxID"), 
		names_to = "Genome", 
		values_to = "TaxID") %>% 
	dplyr::filter(!TaxID %in% host_taxid_hierarchy) %>% 
	mutate(
		TaxID = as.character(TaxID), 
		Genome = gsub("_.*", "", Genome)) %>% 
	group_by(Sequencer, Genome) %>% 
	count(TaxID) %>% 
	arrange(desc(n)) %>% 
	mutate(rank = 1:n()) %>% 
	mutate(TaxID = if_else(rank <= 5, TaxID, "other")) %>% 
	group_by(Sequencer, Genome, TaxID) %>% 
	summarize_all(sum)
write.table(tax_df, "TaxID_summ.tsv", sep = "\t", quote = F, row.names = F)

tax_df <- read.table("TaxID_summ.tsv", sep = "\t", h = T) %>% 
	mutate(TaxID = if_else(rank <= 3, TaxID, "other")) %>% 
	group_by(Sequencer, Genome, TaxID) %>% 
	summarize_all(sum) %>% 
	mutate(TaxID = factor(TaxID)) 
tax_cols <- c("0" = "#A1216B", "131567" = "#618685", "2759" = "#FFEF96", "316" = "#E78AC3", "446" = "#66C2A5", "508771" = "#FFE6E6", "other" = "#87BDD8")
#tax_cols <- c("0" = "#A1216B", "1" = "#E78AC3", "131567" = "#618685", "2759" = "#FFEF96", "33154" = "#002DB3", "316" = "#A31AFF", "446" = "#66C2A5", "1280" = "#FFD92F", "1641" = "#FFB366", "33033" = "#FFE6E6", "40214" = "#E78AC3", "508771" = "#87BDD8", "1969841" = "#A6D854", "other" = "grey30")
arrow_df <- tax_df %>% 
	dplyr::filter(TaxID == "other") %>% 
	ungroup() %>% 
	mutate(y = n, yend = tax_df %>% 
	group_by(Sequencer, Genome) %>% 
	summarize(total = sum(n), .groups = "drop") %>% 
	pull(total)) %>% 
	mutate(perc = round(100*(yend - y)/yend, 2)) %>% 
	mutate(Genome = factor(Genome, levels = c("YH", "T2T", "hg38")))

p_TaxID <- ggplot() + 
	geom_bar(
		data = tax_df,
		mapping =aes(
			x = fct_rev(Genome), 
			y = n, 
			fill = fct_relevel(
				fct_rev(fct_reorder(TaxID, n, .fun = sum)), 
				"other", 
				after = Inf)),
		width = 0.2, 
		stat = "identity", 
		position = "stack") + 
	geom_segment(
		data = arrow_df, 
		mapping = aes(
			x = as.numeric(Genome) - 0.2, 
			xend = as.numeric(Genome) - 0.2, 
			y = y, yend = yend), 
		arrow = arrow(length = unit(0.2, "cm"))) + 
	geom_segment(
		data = arrow_df, 
		mapping = aes(
			x = as.numeric(Genome) - 0.2, 
			xend = as.numeric(Genome) - 0.2, 
			y = yend, 
			yend = y), 
		arrow = arrow(length = unit(0.2, "cm"))) + 
	geom_text(
		data = arrow_df, 
		mapping = aes(
			x = as.numeric(Genome) - 0.4, 
			y = (y + yend)/2, 
			label = paste0(perc, "%")), 
			angle = 90) + 
	labs(y = "Reads", fill = "Taxon") + 
	facet_wrap(.~Sequencer, scales = "free", nrow = 1) + 
	scale_fill_manual(
		values = tax_cols, 
		labels = c(
			"Unclassified", 
			"Cellular organisms", 
			"Eukaryota", 
			expression(italic("S. stutzeri")), 
			expression(italic("L. pneumophila")), 
			expression(paste(italic("T. gondii"), " ME49")), 
			"Other")) + 
	guides(fill = guide_legend(nrow = 3, byrow = TRUE)) + 
	theme_classic() + 
	theme(
		axis.title.x = element_blank(), 
		axis.text.y = element_text(angle = 90, hjust = 0.5), 
		legend.position = "bottom")
ggsave("Uniq_reads_TaxID.pdf", p_TaxID, h = 6, w = 12)


# kmer plot
libary(ggrepel)

for(k in 2:5){
	assign(paste0("df_", k, "mer"), bind_rows(bind_rows(read.table(paste0("kmers/ONT/", k, "mer/hg38.tab")) %>% mutate(Genome = "hg38"), read.table(paste0("kmers/ONT/", k, "mer/T2T.tab")) %>% mutate(Genome = "T2T"), read.table(paste0("kmers/ONT/", k, "mer/YH.tab")) %>% mutate(Genome = "YH")) %>% mutate(Sequencer = "ONT"), bind_rows(read.table(paste0("kmers/ILMN/", k, "/hg38.tab")) %>% mutate(Genome = "hg38"), read.table(paste0("kmers/ILMN/", k, "/T2T.tab")) %>% mutate(Genome = "T2T"), read.table(paste0("kmers/ILMN/", k, "/YH.tab")) %>% mutate(Genome = "YH")) %>% mutate(Sequencer = "ILMN"), bind_rows(read.table(paste0("kmers/MGI/", k, "/hg38.tab")) %>% mutate(Genome = "hg38"), read.table(paste0("kmers/MGI/", k, "/T2T.tab")) %>% mutate(Genome = "T2T"), read.table(paste0("kmers/MGI/", k, "/YH.tab")) %>% mutate(Genome = "YH")) %>% mutate(Sequencer = "MGI")) %>% rename(kmer = V1, Count = V2) %>% group_by(kmer, Genome) %>% summarize(Frequency = sum(Count), .groups = "drop"))
}


p_kmer_top <- (ggplot(df_2mer) + geom_line(aes(x = kmer, y = Frequency, color = Genome, group = Genome), size = 1.5) + labs(x = "2mers") + scale_color_brewer(palette = "Set2") + theme_classic() + theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5), axis.ticks = element_blank(), panel.grid.major.x = element_line(color = "grey90", size = 0.2, linetype = 2)) | 
ggplot(df_3mer) + geom_line(aes(x = kmer, y = Frequency, color = Genome, group = Genome), size = 1.5) + labs(x = "3mers") + scale_color_brewer(palette = "Set2") + theme_classic() + theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5), axis.ticks = element_blank(), panel.grid.major.x = element_line(color = "grey90", size = 0.2, linetype = 2)) | 
ggplot(df_4mer) + geom_line(aes(x = kmer, y = Frequency, color = Genome, group = Genome), size = 1.5) + labs(x = "4mers") + scale_color_brewer(palette = "Set2") + theme_classic() + theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5), axis.ticks = element_blank(), panel.grid.major.x = element_line(color = "grey90", size = 0.2, linetype = 2))) + plot_layout(widths = c(10, 32, 136))

kmer_arrow_df <- data.frame(x = c(58, 207, 217, 270, 509), y = rep(10000000, 5), label = c("AATGG", "ATGGA", "ATTCC", "CATTC", "TGGAA"))
p_kmer_bottom <- ggplot(df_5mer %>% mutate(Genome = factor(Genome, levels = c("YH", "hg38", "T2T")))) + geom_line(aes(x = kmer, y = Frequency, color = Genome, group = Genome), size = 1.5) + geom_text_repel(data = kmer_arrow_df, aes(x = x, y = y, label = label), size = 5, box.padding = 1.5, min.segment.length = 0, arrow = arrow(length = unit(0.3, "cm"))) + labs(x = "Different 5mers") + scale_color_manual(values = c("hg38" = "#66C2A5", "T2T" = "#FC8D62", "YH" = "#8DA0CB")) + theme_classic() + theme(axis.text = element_blank(), axis.ticks = element_blank())

p_kmer <- p_kmer_top / p_kmer_bottom + plot_layout(guides = "collect") & theme(legend.position = "bottom")
#ggsave("Figure S6 Spectrum of kmers in false negative reads for each genome.pdf", p_kmer_top / p_kmer_bottom + plot_layout(guides = "collect") & theme(legend.position = "bottom"), h = 6, w = 20)


# Ideogram
"Adobe AI 
x coordinates: 
	right column (YH): 50*chr_no
	left column (hg38): -12.5
	left annotation (sat): -20.5
y coordinates:
	-7.5 chr_no
	498 chr_no y
	300.64 1y
	305.07 2y
	337.47 3y
	343.3 4y
	352.305 5y
	360.025 6y
	368.985 7y
	380.155 8y
	376.76 9y
	389.105 10y
	388.78 11y
	390.185 12y
	405.595 13y
	415.215 14y
	416.345 15y
	419.01 16y
	428.35 17y
	431.295 18y
	445.95 19y
	442.41 20y
	458.85 21y
	453.995 22y
	373.89 xy
	445.33 yy
不含label高 397.57 中点y坐标 300.675
不含label不含编号高 386.59 中点y坐标 295.185 sat中点x坐标 604.5

总图：
	总体中点x坐标 616.0669 
	编号x 617.815
	与chr1柱子dy 202.815
	y1 217.2678
	编号y1^ 420.0828
	y2 637.6953
	编号y2^ 840.5103
	y3 1057.6953
	编号y3^ 1260.5103
"

#library(grImport2)
#library(rsvg)

#rsvg_svg("ideogram/final/Ideogram.svg", "ideogram/final/Ideo_Cairo.svg")
#ideo_graph <- readPicture("ideogram/final/Ideo_Cairo.svg")

#g2 <- plot_grid(plot_grid(pictureGrob(ideo_graph), plot_grid(p_GC, p_TaxID, ncol = 1, rel_heights = c(1, 1.5), labels = c("B", "C")), nrow = 1, rel_widths = c(2.75, 1), labels = c("A", "")), p_kmer, ncol = 1, rel_heights = c(5, 1), labels = c("", "D"))
#ggsave("test.pdf", plot_grid(plot_grid(ggplot() + geom_blank() + theme_void(), plot_grid(ggplot() + geom_blank() + theme_void(), p_TaxID, ncol = 1, rel_heights = c(1, 1), labels = c("B", "C")), nrow = 1, rel_widths = c(2.75, 1), labels = c("A", "")), p_kmer, ncol = 1, rel_heights = c(3, 1), labels = c("", "D")), w = 20, h = 18)
g2 <- plot_grid(plot_grid(ggplot() + geom_blank() + theme_void(), plot_grid(p_GC, p_TaxID, ncol = 1, rel_heights = c(1, 1), labels = c("B", "C")), nrow = 1, rel_widths = c(2.75, 1), labels = c("A", "")), p_kmer, ncol = 1, rel_heights = c(3, 1), labels = c("", "D"))
ggsave("Figure2_raw.pdf", g2, w = 20, h = 18)

