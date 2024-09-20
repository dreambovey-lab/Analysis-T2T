library(RColorBrewer)
library(patchwork)
options(warn = -1)
	
performance_measures = c("sens", "spec", "MCC")
sequencers = c("ILMN", "MGI", "ONT")

ONT = read.table("ONT_performance_stats.tsv", sep = "\t", h = T) %>% mutate(Lab = gsub("POOLING[12]_", "", Lab), Sample = factor(Sample, levels = c("NC_C", "NC_HD", "NC_nHD", "S1_D", "S2_D", "S3_D", "S4_D", "S5_D", "S6_D", "S7_D", "S8_D", "D0", "D1-1", "D1-2", "D1-3", "D1-2(10)", "D1-1(100)", "D2-1", "D2-2", "D2-3", "D2-2(10)", "D2-1(100)"))) %>% arrange(Lab, Sample) %>% mutate(Serial = 1:n(), .after = "Sequencer")
#write.xlsx(ONT, "ONT.xlsx") 
#ONT = read.xlsx("ONT.xlsx", sheet = 1)
ILMN = read.xlsx("ILMN.xlsx", sheet = 1) %>% mutate(Lab = gsub("_.*", "", Lab))
MGI = read.xlsx("MGI.xlsx", sheet = 1) %>% mutate(Lab = gsub("_.*", "", Lab))

p_values = lapply(performance_measures, function(measure) {
	unlist(lapply(sequencers, function(sequencer){
		df = eval(sym(sequencer)) %>% select(Lab, Sample, ends_with(measure)) %>% unite("Sample", c(Sample, Lab)) %>% mutate(Sample = factor(Sample, levels = unique(Sample))) %>% pivot_longer(ends_with(measure), "Method") %>% mutate(Method = gsub("bowtie2_v", "bowtie2@v", gsub("bwa_m", "bwa@m", Method))) %>% separate(Method, c("Genome", "Algorithm", "Measure"), sep = "_") %>% mutate(Algorithm = gsub("bowtie2@v", "bowtie2_v", gsub("bwa@m", "bwa_m", Algorithm))) %>% mutate(Genome = factor(Genome, levels = c("hg38", "T2T", "YH")), Algorithm = factor(Algorithm, levels = c("bowtie2_vf", "bowtie2_vs", "bwa_mem", "minimap2", "winnowmap", "STAT", "kraken2")))
		#ddf = eval(sym(sequencer)) %>% select(Lab, Sample, ends_with(measure)) %>% mutate(across(ends_with(measure), function(x) x - eval(sym(paste0("T2T_kraken2_", measure))))) %>% unite("Sample", c(Sample, Lab)) %>% mutate(Sample = factor(Sample, levels = unique(Sample))) %>% pivot_longer(ends_with(measure), "Method") %>% mutate(Method = gsub("bowtie2_v", "bowtie2@v", gsub("bwa_m", "bwa@m", Method))) %>% separate(Method, c("Genome", "Algorithm", "Measure"), sep = "_") %>% mutate(Algorithm = gsub("bowtie2@v", "bowtie2_v", gsub("bwa@m", "bwa_m", Algorithm))) %>% mutate(Genome = factor(Genome, levels = c("hg38", "T2T", "YH")), Algorithm = factor(Algorithm, levels = c("bowtie2_vf", "bowtie2_vs", "bwa_mem", "minimap2", "winnowmap", "STAT", "kraken2")))		
		#pheatmap(ddf %>% group_by(Genome, Algorithm) %>% summarize(measure = median(value), .groups = "drop") %>% pivot_wider(names_from = "Algorithm", values_from = "measure") %>% column_to_rownames("Genome"), main = paste0(measure, "_", sequencer), legend = F, cluster_rows = F, cluster_cols = F, display_numbers = T, number_format = "%.4f", color = colorRampPalette(brewer.pal(n = 3, name ="Greens"))(100), filename = paste0("all_samples/", measure, "_", sequencer, ".png"), width = 2, height = 2, fontsize = 8)
		kw = kruskal.test(df$value, interaction(df$Genome, df$Algorithm))
		pw_wil_g = pairwise.wilcox.test(df$value, interaction(df$Genome, df$Algorithm), paired = T, p.adj = "bonf", alternative = "greater")
		write.table(pw_wil_g$p.value, paste0("stats_tests/genome_cross_algorithm/", measure, "_", sequencer, "_greater.txt"), sep = "\t", quote = F, row.names = T)
		pw_wil_l = pairwise.wilcox.test(df$value, interaction(df$Genome, df$Algorithm), paired = T, p.adj = "bonf", alternative = "less")
		write.table(pw_wil_l$p.value, paste0("stats_tests/genome_cross_algorithm/", measure, "_", sequencer, "_less.txt"), sep = "\t", quote = F, row.names = T)
		list(measure, sequencer, kw$p.value)
	}))
})

#ONT Kruskal-Wallis p-value = 0.0002099
 "sens"                         "ONT"                         "2.76589811870283e-21"         
 "sens"                         "ILMN"                         "1.55434197825468e-154"       
 "sens"                         "MGI"                         "8.19428202467519e-119"        
 "spec"                        "ONT"                          "1.21845234915736e-25"        
"spec"                         "ILMN"                        "1.45858473312019e-52"         
"spec"                         "MGI"                          "7.23533882751251e-56"        
"PPV"                          "ONT"                         "3.29889088374883e-16"         
"PPV"                          "ILMN"                         "8.02195133520427e-75"        
"PPV"                          "MGI"                         "2.42613393900818e-79"         
"NPV"                         "ONT"                          "0.00000000000130501162841061"
"NPV"                          "ILMN"                        "6.08921915134577e-164"        
"NPV"                         "MGI"                          "9.24660143780548e-139"       
"MCC"                          "ONT"                         "0.00223645918689595"          
"MCC"                         "ILMN"                         "1.09840396472437e-165"       
"MCC"                          "MGI"                         "1.1515908350305e-134"  

#per-sample heatmap
lapply(performance_measures, function(measure) {
	lapply(sequencers, function(sequencer){
		for(genome in c("hg38", "T2T", "YH")){
			if(genome == "hg38") {my_color = "BuPu"}
			else if(genome == "T2T") {my_color = "Greens"}
			else {my_color = "Oranges"}
			df = eval(sym(sequencer)) %>% select(Lab, Sample, Serial, ends_with(measure)) %>% dplyr::filter(!grepl("Undetermined", Sample)) %>% unite("Sample", c(Sample, Lab)) %>% mutate(Sample = factor(Sample, levels = unique(Sample))) %>% pivot_longer(ends_with(measure), "Method") %>% mutate(Method = gsub("bowtie2_v", "bowtie2@v", gsub("bwa_m", "bwa@m", Method))) %>% separate(Method, c("Genome", "Algorithm", "Measure"), sep = "_") %>% mutate(Algorithm = gsub("bowtie2@v", "bowtie2_v", gsub("bwa@m", "bwa_m", Algorithm))) %>% mutate(Genome = factor(Genome, levels = c("hg38", "T2T", "YH")), Algorithm = factor(Algorithm, levels = c("bowtie2_vf", "bowtie2_vs", "bwa_mem", "minimap2", "winnowmap", "STAT", "kraken2")))
			pheatmap(df %>% dplyr::filter(Genome == genome) %>% select(-c(Genome, Measure)) %>% pivot_wider(names_from = "Algorithm", values_from = "value") %>% dplyr::filter(!grepl("NA", Sample)) %>% arrange(Serial) %>% select(-Serial) %>% column_to_rownames("Sample"), cluster_row = F, cluster_cols = F, display_numbers = T, number_format = "%.4f", color = colorRampPalette(brewer.pal(n = 3, name = my_color))(100), filename = paste0("per_sample/", measure, "/", measure, "_", sequencer, "_", genome, "_per_sample.png"), width = 6, height = 12, fontsize = 10)
		}
	})
})


#line graph
lapply(performance_measures, function(measure) {
	lapply(sequencers, function(sequencer){
		df = eval(sym(sequencer)) %>% select(Lab, Sample, Serial, ends_with(measure)) %>% dplyr::filter(!grepl("Undetermined", Sample)) %>% mutate(SampleType = case_when(grepl("clinical", Lab) ~ "clinical", grepl("Lab", Lab) ~ "D_series", grepl("nHD", Lab) ~ "D_series_nHD", grepl("_HD", Lab) ~ "D_series_HD", TRUE ~ "NA")) %>% mutate(Sample = factor(Sample, levels = unique(Sample))) %>% pivot_longer(ends_with(measure), "Method") %>% mutate(Method = gsub("bowtie2_v", "bowtie2@v", gsub("bwa_m", "bwa@m", Method))) %>% separate(Method, c("Genome", "Algorithm", "Measure"), sep = "_") %>% mutate(Algorithm = gsub("bowtie2@v", "bowtie2_v", gsub("bwa@m", "bwa_m", Algorithm))) %>% mutate(Genome = factor(Genome, levels = c("hg38", "T2T", "YH")), Algorithm = factor(Algorithm, levels = c("bowtie2_vf", "bowtie2_vs", "bwa_mem", "minimap2", "winnowmap", "STAT", "kraken2")))
		#df %>% group_by(SampleType) %>% group_walk(~ ggsave(paste0("line_graph/compare_genomes/", measure, "_", sequencer, "_", .y$SampleType, ".png"), ggplot(.x %>% group_by(Sample, Genome, Algorithm) %>% summarize(value = mean(value), .groups = "drop")) + geom_line(aes(x = Sample, y = value, group = Genome, color = Genome), alpha = 0.8) + labs(y = measure) + theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + facet_wrap(.~Algorithm, scales = "free_y"), height = 10, width = 10))
		df %>% group_by(SampleType) %>% group_walk(~ ggsave(paste0("line_graph/compare_algorithms/", measure, "_", sequencer, "_", .y$SampleType, ".png"), ggplot(.x %>% group_by(Sample, Genome, Algorithm) %>% summarize(value = mean(value), .groups = "drop")) + geom_line(aes(x = Sample, y = value, group = Algorithm, color = Algorithm), alpha = 0.8) + labs(y = measure) + theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + facet_wrap(.~Genome, scales = "free_y"), height = 10, width = 10))
		return(NULL)
	})
})

#line graph 2
lapply(sequencers, function(sequencer){
	plots = unlist(lapply(performance_measures, function(measure) {
		df = eval(sym(sequencer)) %>% select(Lab, Sample, Serial, ends_with(measure)) %>% dplyr::filter(!grepl("(Undetermined|NC|_R)", Sample)) %>% mutate(SampleType = case_when(grepl("clinical", Lab) ~ "clinical", grepl("Lab", Lab) ~ "D_series", grepl("nHD", Lab) ~ "D_series_nHD", grepl("[^n]HD", Lab) ~ "D_series_HD", TRUE ~ "NA")) %>% mutate(Sample = factor(Sample, levels = unique(Sample))) %>% pivot_longer(ends_with(measure), "Method") %>% mutate(Method = gsub("bowtie2_v", "bowtie2@v", gsub("bwa_m", "bwa@m", Method))) %>% separate(Method, c("Genome", "Algorithm", "Measure"), sep = "_") %>% mutate(Algorithm = gsub("bowtie2@v", "bowtie2_v", gsub("bwa@m", "bwa_m", Algorithm))) %>% mutate(Genome = factor(Genome, levels = c("hg38", "T2T", "YH")), Algorithm = factor(Algorithm, levels = c("bowtie2_vf", "bowtie2_vs", "bwa_mem", "minimap2", "winnowmap", "STAT", "kraken2")))
		#df %>% group_by(SampleType) %>% group_walk(~ ggsave(paste0("line_graph/compare_genomes/", measure, "_", sequencer, "_", .y$SampleType, ".png"), ggplot(.x %>% group_by(Sample, Genome, Algorithm) %>% summarize(value = mean(value), .groups = "drop")) + geom_line(aes(x = Sample, y = value, group = Genome, color = Genome), alpha = 0.8) + labs(y = measure) + theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + facet_wrap(.~Algorithm, scales = "free_y"), height = 10, width = 10))
		return(df %>% group_by(SampleType) %>% group_map(~ggplot(.x %>% group_by(Sample, Genome, Algorithm) %>% summarize(value = mean(value), .groups = "drop")) + geom_line(aes(x = Sample, y = value, group = Algorithm, color = Algorithm), alpha = 0.8, size = 2) + labs(y = measure, x = .y$SampleType) + scale_color_manual(values = brewer.pal(n = length(unique(.x$Algorithm)), name = "Dark2")) + theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + facet_wrap(.~Genome, scales = "free_y")))
	}), recursive = F)
	if(sequencer == "ONT") {my_width = 25}
	else {my_width = 20}
	ggsave(paste0("line_graph/wrapped/compare_algorithms_", sequencer, ".png"), wrap_plots(plots, nrow = 3, guides = "collect", tag_level = "keep"), height = 15, width = my_width)
})
lapply(sequencers, function(sequencer){
	plots = unlist(lapply(performance_measures, function(measure) {
		df = eval(sym(sequencer)) %>% select(Lab, Sample, Serial, ends_with(measure)) %>% dplyr::filter(!grepl("(Undetermined|NC|_R)", Sample)) %>% mutate(SampleType = case_when(grepl("clinical", Lab) ~ "clinical", grepl("Lab", Lab) ~ "D_series", grepl("nHD", Lab) ~ "D_series_nHD", grepl("[^n]HD", Lab) ~ "D_series_HD", TRUE ~ "NA")) %>% mutate(Sample = factor(Sample, levels = unique(Sample))) %>% pivot_longer(ends_with(measure), "Method") %>% mutate(Method = gsub("bowtie2_v", "bowtie2@v", gsub("bwa_m", "bwa@m", Method))) %>% separate(Method, c("Genome", "Algorithm", "Measure"), sep = "_") %>% mutate(Algorithm = gsub("bowtie2@v", "bowtie2_v", gsub("bwa@m", "bwa_m", Algorithm))) %>% mutate(Genome = factor(Genome, levels = c("hg38", "T2T", "YH")), Algorithm = factor(Algorithm, levels = c("bowtie2_vf", "bowtie2_vs", "bwa_mem", "minimap2", "winnowmap", "STAT", "kraken2")))
		#df %>% group_by(SampleType) %>% group_walk(~ ggsave(paste0("line_graph/compare_genomes/", measure, "_", sequencer, "_", .y$SampleType, ".png"), ggplot(.x %>% group_by(Sample, Genome, Algorithm) %>% summarize(value = mean(value), .groups = "drop")) + geom_line(aes(x = Sample, y = value, group = Genome, color = Genome), alpha = 0.8) + labs(y = measure) + theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + facet_wrap(.~Algorithm, scales = "free_y"), height = 10, width = 10))
		return(df %>% group_by(SampleType) %>% group_map(~ggplot(.x %>% group_by(Sample, Genome, Algorithm) %>% summarize(value = mean(value), .groups = "drop")) + geom_line(aes(x = Sample, y = value, group = Genome, color = Genome), alpha = 0.8, size = 2) + labs(y = measure, x = .y$SampleType) + scale_color_manual(values = brewer.pal(n = length(unique(.x$Genome)), name = "Dark2")) + theme_bw() + theme(axis.text.x = element_text(angle = 60, hjust = 1)) + facet_wrap(.~Algorithm, scales = "free_y")))
	}), recursive = F)
	if(sequencer == "ONT") {my_width = 25}
	else {my_width = 20}
	ggsave(paste0("line_graph/wrapped/compare_genomes_", sequencer, ".png"), wrap_plots(plots, nrow = 3, guides = "collect", tag_level = "keep"), height = 15, width = my_width)
})