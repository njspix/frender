library(tidyverse)

reads_by_type = function(file, top_n_perc= 0.05){
  
  expt_name = str_remove(file, "^.*mismatches_") %>% str_remove(".csv")
  data = read_csv(file, col_types = "ccccffdl")
  
  all_demux_ok = data %>% 
    pull(demux_ok) %>% all()
  
  demux_message = if_else(all_demux_ok, "", "\nWARNING: Some files appear to be incorrectly demuxed!")
  
  temp = if_else(all_demux_ok, "", paste0("\nThe actual distribution of barcodes among demux files is shown in red; the correct distribution is shown in green."))
  
  caption_message = paste0("All reads found in supplied files or directory are shown grouped by type.\n", 
                           "Barcodes comprising more than ", top_n_perc*100, "% of reads in their category are shown separately.", temp)
  
  data %>% 
    filter(read_type != "demuxable") %>% 
    arrange(desc(reads)) %>% 
    mutate(keep = if_else(reads > sum(reads)*top_n_perc, T, F)) %>% 
    arrange(desc(reads)) -> non_demuxable
  
  non_demuxable %>% 
    filter(!keep) %>% 
    group_by(read_type) %>% 
    summarize(read_type = read_type, reads = sum(reads), demux_ok = all(demux_ok), keep = F) %>% 
    distinct() -> non_demux_summary
  
  data %>% 
    filter(read_type == 'demuxable') %>% 
    group_by(sample_name) %>% 
    summarize(matched_idx1 = matched_idx1, matched_idx2 = matched_idx2, read_type = read_type, sample_name = sample_name, reads = sum(reads), demux_ok = all(demux_ok)) %>% 
    distinct() -> demux_summary
  
  data %>% 
    filter(read_type == 'demuxable' & demux_ok) %>% 
    group_by(sample_name) %>% 
    summarize(matched_idx1 = matched_idx1, matched_idx2 = matched_idx2, read_type = read_type, sample_name = sample_name, reads = sum(reads), demux_ok = all(demux_ok)) %>% 
    distinct() %>% 
    mutate(demux_ok= as.logical(all_demux_ok^demux_ok)) %>% 
    bind_rows(non_demux_summary) %>% 
    bind_rows(non_demuxable %>% filter(keep)) %>% 
    ungroup()%>% 
    mutate(reads = reads/1E6,
           label = case_when(
             read_type == 'demuxable' ~ "",
             read_type == 'undetermined' ~ if_else(is.na(idx1), "all others", paste0(idx1, "+\n", idx2)),
             read_type == 'index_hop' ~ if_else(is.na(matched_idx1), "all others", paste0(matched_idx1, "+\n", matched_idx2)),
             read_type == 'ambiguous' ~ if_else(is.na(matched_idx1), "all others", paste0(matched_idx1, "+\n", matched_idx2))
           ),
           sample_name = case_when(
             label == 'all others' ~ label,
             is.na(sample_name) ~ as.character(row_number()),
             TRUE ~ as.character(sample_name)
           ),
           sample_name = fct_reorder(sample_name, reads) %>% fct_relevel("all others", after = 0)) %>% 
    arrange(read_type, desc(sample_name)) %>% 
    group_by(read_type) %>% 
    mutate(cumulative_sum = cumsum(reads), label_height = cumulative_sum-(0.5*reads),
           total = sum(reads)) %>% 
    ungroup() %>% 
    mutate(read_type = fct_reorder(read_type, total, .desc = T)) %>% 
    arrange(read_type)-> demux_ok_summary
  
  non_demuxable %>% 
    filter(keep) %>% 
    bind_rows(non_demux_summary) %>% 
    bind_rows(demux_summary) %>% 
    ungroup()%>% 
    mutate(reads = reads/1E6,
           label = case_when(
             read_type == 'demuxable' ~ "",
             read_type == 'undetermined' ~ if_else(is.na(idx1), "all others", paste0(idx1, "+\n", idx2)),
             read_type == 'index_hop' ~ if_else(is.na(matched_idx1), "all others", paste0(matched_idx1, "+\n", matched_idx2)),
             read_type == 'ambiguous' ~ if_else(is.na(matched_idx1), "all others", paste0(matched_idx1, "+\n", matched_idx2))
           ),
           sample_name = case_when(
             label == 'all others' ~ label,
             is.na(sample_name) ~ as.character(row_number()),
             TRUE ~ as.character(sample_name)
           ),
           sample_name = fct_reorder(sample_name, reads) %>% fct_relevel("all others", after = 0)) %>% 
    arrange(read_type, desc(sample_name)) %>% 
    group_by(read_type) %>% 
    mutate(cumulative_sum = cumsum(reads), label_height = cumulative_sum-(0.5*reads),
           total = sum(reads)) %>% 
    ungroup() %>% 
    mutate(read_type = fct_reorder(read_type, total, .desc = T)) %>% 
    arrange(read_type) -> plot_data
  
  if (all_demux_ok){
    test = F
  } else { test = expr(guide_legend())}
  
  ggplot(plot_data)+
    geom_col(aes(read_type, reads, fill = read_type))+
    scale_fill_manual(values = c('demuxable' = "#7FC97F", 'undetermined'= '#BEAED4', 'index_hop' = '#FDC086', 'ambiguous' = '#FFFF99'))+
    geom_col(data = demux_ok_summary, aes(read_type, reads, color = demux_ok), fill = NA)+
    scale_color_manual(values = c('FALSE' = 'red', 'TRUE' = 'white'))+
    geom_text(aes(read_type, label_height, label = label), size = 3, vjust = "center")+
    labs(fill = "Read type", x = "", y = "Reads (million)", color = "Correctly\ndemuxed?", caption = caption_message)+
    guides(col= eval(test))+
    ggtitle(paste0(expt_name, demux_message))+
    theme(plot.caption = element_text(hjust = 0))
}

barcodes_by_prevalence = function(file, cutoff = 0.99){
  
  expt_name = str_remove(file, "^.*mismatches_") %>% str_remove(".csv")
  data = read_csv(file, col_types = "ccccffdl")
  
  total_reads = data %>% 
    pull(reads) %>% sum()
  
  all_demux_ok = data %>% 
    pull(demux_ok) %>% all()
  
  demux_message = if_else(all_demux_ok, "All files appear to be correctly demuxed", "WARNING! Some files appear to be incorrectly demuxed")
  
  data %>% 
    mutate(matched = paste0(matched_idx1, "+", matched_idx2)) %>% 
    group_by(matched)%>% 
    summarize(read_type = read_type, sample_name = sample_name, reads = sum(reads), demux_ok = all(demux_ok)) %>% 
    distinct() %>% 
    arrange(desc(reads)) %>% 
    mutate(name = case_when(
      matched == "NA+NA" ~ "undetermined",
      !is.na(sample_name) ~ as.character(sample_name),
      TRUE ~ matched
    )) %>% 
    ungroup() %>% 
    mutate(test = reads/sum(reads),
           test2 = cumsum(test),
           test3 = test2 < cutoff) %>% 
    filter(as.logical(test3)) %>% 
    select(name, read_type, reads, demux_ok) -> filtered
  
  filtered %>% 
    mutate(name = fct_reorder(as_factor(name), reads)) %>% 
    ggplot()+
    geom_col(aes(name, reads/1E6, fill = read_type))+
    geom_col(aes(name, reads/1E6, fill = NA, color = demux_ok))+
    coord_flip()+
    ggtitle(expt_name, subtitle = paste0("Files shown account for ", round((sum(filtered$reads))*100/total_reads, 2),"% of all reads\n", demux_message))+
    labs(x = "", y = "Reads (million)", fill = "Read type", color = "Correctly\ndemuxed?")+
    scale_fill_manual(values = c('demuxable' = "#7FC97F", 'undetermined'= '#BEAED4', 'index_hop' = '#FDC086', 'ambiguous' = '#FFFF99'))+
    scale_color_manual(values = c('TRUE' = 'white', 'FALSE' = 'red'))+
    guides(color = guide_none())
}

pdf("plots.pdf")

for (file in list.files(pattern = "frender-scan-results.*.csv", full.names = T)) {
  print(reads_by_type(file))
  print(barcodes_by_prevalence(file))
}
dev.off()
