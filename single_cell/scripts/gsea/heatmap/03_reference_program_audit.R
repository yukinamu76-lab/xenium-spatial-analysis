suppressPackageStartupMessages({library(dplyr);library(tidyr);library(readr);library(stringr);library(purrr);library(tibble)})
source("config.R")
gsea<-read_csv(file.path(output_dir,"01_all_clusters_GSEA_long.csv"),show_col_types=FALSE)|>mutate(cluster=as.character(cluster),padj=pmax(as.numeric(padj),1e-300),NES=as.numeric(NES))
classification<-read_csv(file.path(output_dir,"02_classification_used_v3.csv"),show_col_types=FALSE)
program_hits<-imap_dfr(reference_program_patterns,function(pattern,program){classification|>filter(str_detect(str_to_upper(pathway),pattern)|str_detect(str_to_upper(pathway_label),pattern))|>transmute(program,pathway,pathway_label,functional_block,eligible_for_v3_selection,exclusion_reason)})|>distinct()
audit_long<-program_hits|>
  left_join(
    gsea |> select(pathway,cluster,NES,padj) |> distinct(pathway, cluster, .keep_all = TRUE),
    by="pathway",
    relationship="many-to-many"
  ) |>
  distinct(program, pathway, cluster, .keep_all = TRUE) |>
  mutate(significant=padj<fdr_cutoff&abs(NES)>=nes_cutoff)|>arrange(program,pathway,as.integer(cluster))
write_csv(audit_long,file.path(output_dir,"08_reference_program_audit_long_v3.csv"))
audit_wide<-audit_long|>select(program,functional_block,pathway,pathway_label,eligible_for_v3_selection,exclusion_reason,cluster,NES,padj)|>pivot_wider(names_from=cluster,values_from=c(NES,padj),names_glue="{.value}_C{cluster}")
write_csv(audit_wide,file.path(output_dir,"09_reference_program_audit_wide_v3.csv"))
program_summary<-audit_long|>group_by(program,cluster)|>summarise(n_pathways=n_distinct(pathway),n_significant=sum(significant,na.rm=TRUE),median_NES=median(NES,na.rm=TRUE),max_abs_NES=max(abs(NES),na.rm=TRUE),best_pathway=pathway[which.max(abs(NES))][1],best_NES=NES[which.max(abs(NES))][1],best_FDR=padj[which.max(abs(NES))][1],.groups="drop")|>arrange(program,as.integer(cluster))
write_csv(program_summary,file.path(output_dir,"10_reference_program_cluster_summary_v3.csv"))
message("v3 reference-program audit completed.")
