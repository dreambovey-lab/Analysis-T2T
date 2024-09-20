ONT = read.xlsx("ONT.xlsx", sheet = 1)
ILMN = read.xlsx("ILMN.xlsx", sheet = 1) %>% mutate(Lab = gsub("_.*", "", Lab))
MGI = read.xlsx("MGI.xlsx", sheet = 1) %>% mutate(Lab = gsub("_.*", "", Lab))

my_colors = c(TP = "#D8B365", FN = "#EBD9B2", TN = "#5AB4AC", FP = "#ACD9D5")

draw_benchmarking_plot = function(df, my_lab){
    ggplot(df %>% dplyr::filter(Lab == my_lab)) + geom_bar(aes(x = fct_rev(Method), y = value, fill = Metrics, group = Metrics), stat = "identity", position = position_stack()) + geom_hline(yintercept = 0) + labs(x = "Method", y = "Number of reads") + ggtitle(my_lab) + scale_fill_manual(values = c(TP = "#D8B365", FN = "#EBD9B2", TN = "#5AB4AC", FP = "#ACD9D5")) + scale_y_continuous(labels = function(x) abs(x), breaks = scales::extended_breaks(6)) + theme(axis.ticks.y = element_blank(), panel.background = element_blank(), panel.border = element_rect(color = "black", fill = NA)) + coord_flip() + facet_wrap(.~Sample, strip.position = "right", ncol = 1, dir = "v")
}

df_ONT = ONT %>% dplyr::filter(!grepl("NC", Sample)) %>% mutate(Lab = gsub("POOLING[12]_", "", Lab)) %>% mutate(across(ends_with("TP"), function(x) x - Host_consensus), across(ends_with("TN"), function(x) x - Non_host_consensus)) %>% select(c(Sample, Lab, ends_with(c("TP", "FP", "TN", "FN")))) %>% pivot_longer(-c(Sample, Lab)) %>% separate(name, c("Method", "Metrics"), sep = "__") %>% mutate(value = if_else(grepl("TN|FP", Metrics), -value, value), Method = factor(Method, levels = c("T2T_minimap2", "T2T_winnowmap", "T2T_kraken2", "hg38_minimap2", "hg38_winnowmap", "hg38_kraken2", "hg38_STAT", "YH_minimap2", "YH_winnowmap", "YH_kraken2")), Metrics = factor(Metrics, levels = c("TN", "TP", "FN", "FP")))

g_ONT = wrap_plots(draw_benchmarking_plot(df_ONT, "clinical"), draw_benchmarking_plot(df_ONT, "HD"), draw_benchmarking_plot(df_ONT, "nHD"), nrow = 1, guides = "collect")

df_ILMN = ILMN %>% dplyr::filter(!grepl("NC", Sample)) %>% mutate(across(ends_with("TP"), function(x) x - Host_consensus), across(ends_with("TN"), function(x) x - Non_host_consensus)) %>% select(c(Sample, Lab, ends_with(c("TP", "FP", "TN", "FN")))) %>% pivot_longer(-c(Sample, Lab)) %>% separate(name, c("Method", "Metrics"), sep = "__") %>% mutate(value = if_else(grepl("TN|FP", Metrics), -value, value), Method = factor(Method, levels = c("T2T_bowtie2_vf", "T2T_bowtie2_vs", "T2T_bwa_mem", "T2T_kraken2", "hg38_bowtie2_vf", "hg38_bowtie2_vs",  "hg38_bwa_mem", "hg38_kraken2", "hg38_STAT", "YH_bowtie2_vf", "YH_bowtie2_vs", "YH_bwa_mem", "YH_kraken2")), Metrics = factor(Metrics, levels = c("TN", "TP", "FN", "FP")))

g_ILMN = wrap_plots(draw_benchmarking_plot(df_ILMN, "clinical"), draw_benchmarking_plot(df_ILMN, "Lab4"), draw_benchmarking_plot(df_ILMN, "Lab34"), draw_benchmarking_plot(df_ILMN, "Lab40"), draw_benchmarking_plot(df_ILMN, "Lab49"), nrow = 1, guides = "collect")

df_MGI = MGI %>% dplyr::filter(!grepl("NC", Sample)) %>% mutate(across(ends_with("TP"), function(x) x - Host_consensus), across(ends_with("TN"), function(x) x - Non_host_consensus)) %>% select(c(Sample, Lab, ends_with(c("TP", "FP", "TN", "FN")))) %>% pivot_longer(-c(Sample, Lab)) %>% separate(name, c("Method", "Metrics"), sep = "__") %>% mutate(value = if_else(grepl("TN|FP", Metrics), -value, value), Method = factor(Method, levels = c("T2T_bowtie2_vf", "T2T_bowtie2_vs", "T2T_bwa_mem", "T2T_kraken2", "hg38_bowtie2_vf", "hg38_bowtie2_vs",  "hg38_bwa_mem", "hg38_kraken2", "hg38_STAT", "YH_bowtie2_vf", "YH_bowtie2_vs", "YH_bwa_mem", "YH_kraken2")), Metrics = factor(Metrics, levels = c("TN", "TP", "FN", "FP")))

g_MGI = wrap_plots(draw_benchmarking_plot(df_MGI, "clinical"), draw_benchmarking_plot(df_MGI, "Lab8"), draw_benchmarking_plot(df_MGI, "Lab22"), draw_benchmarking_plot(df_MGI, "Lab28"), draw_benchmarking_plot(df_MGI, "Lab43"), nrow = 1, guides = "collect")

ggsave("benchmarking.pdf", g_ONT + g_ILMN + g_MGI + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A"), height = 16, width = 16)
