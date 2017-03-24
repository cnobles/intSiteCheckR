#' Serial processing of integration site data will use one core by definition
#' and may take longer to analyze a requested dataset.
process_intsite_data <- function(uniq_sites_data, specimen_data, args, 
                                 config, gene_references){
  uniq_sites_data <- uniq_sites_data[, !duplicated(colnames(uniq_sites_data))]
  uniq_sites_data$gtsp <- sapply(
    strsplit(as.character(uniq_sites_data$sampleName), "-"), "[[", 1)
  uniq_sites_data <- left_join(
    uniq_sites_data,
    specimen_data[,c("patient", "celltype", "timepoint", "specimenaccnum")],
    by = c("gtsp" = "specimenaccnum"))

  process_data <- split(uniq_sites_data, uniq_sites_data$sampleName)
  
  #' Refine breakpoints for each replicate independently
  process_data <- unname(unlist(GRangesList(lapply(
    process_data, 
    function(sites, max_gap){
      sites <- db_to_granges(sites)
      mcols(sites) <- mcols(sites)[,c("samplename", "specimen", "patient", 
                                      "celltype", "timepoint", "count")]
      sites <- suppressMessages(
        refine_breakpoints(sites, sata.gap = max_gap, counts = "count"))
      sites$called.bp <- sites$adj.bp <- NULL
      unique_granges(sites, sum.counts = TRUE, counts.col = "count")
      }, max_gap = config$breakpoint_window
  ))))
  
  #' Standardize integration site positions across each patient
  process_data <- split(process_data, process_data$patient)
  process_data <- unname(unlist(GRangesList(lapply(
    process_data, 
    function(sites, max_gap){
      sites <- suppressMessages(
        standardize_intsites(sites, sata.gap = max_gap))
      sites$clusID <- sites$called.pos <- sites$adj.pos <- NULL
      unique_granges(sites, sum.counts = TRUE, counts.col = "count")
    }, max_gap = config$standardize_window
  ))))
  
  #' Condense the data to integration site locations and determine each clone's 
  #' abundance by the specified method.
  if(args$abundance_method != "readCount"){
    if(args$abundance_method == "uniqFragLength"){
      method <- "fragLen"
    }else if(args$abundance_method == "sonicAbundance"){
      method <- "estAbund"
    }
    
    cond_sites <- suppressMessages(suppressWarnings(try(condense_intsites(
        process_data, grouping = "specimen", return.abundance = TRUE, 
        method = method, replicates = "samplename"), 
      silent = TRUE)))
    
    if(class(cond_sites) == "try-error"){
      cond_sites <- suppressMessages(suppressWarnings(try(condense_intsites(
        process_data, grouping = "specimen", return.abundance = TRUE, 
        method = "fragLen", replicates = "samplename"), 
        silent = TRUE)))
      if(config$print_deviations){
        pandoc.emphasis(descript$sonicAbundFail)}
    }
    
    cond_sites$samplename <- cond_sites$count <- NULL
  }else{
    cond_sites <- flank(process_data, -1, start = TRUE)
    cond_sites$samplename <- NULL
    cond_sites <- unique_granges(
      cond_sites, sum.counts = TRUE, counts.col = "count")
    cond_sites$estAbund <- as.integer(cond_sites$count)
    cond_sites$relAbund <- cond_sites$estAbund/sum(cond_sites$estAbund)
    cond_sites$relRank <- rank(-1*cond_sites$relAbund, ties.method = "max")
    cond_sites$posID <- generate_posid(cond_sites)
    cond_sites$count <- NULL
  }
  
  if(args$position_ids_exclusively){
    cond_sites$relAbund <- cond_sites$relRank <- NULL
  }
  
  #' If Reference Genes are available, give gene names and markers to 
  #' integration sites.
  if(class(gene_references$ref_genes) != "try-error"){
    #' Annotate Sites with Gene names for within and nearest genes
    cond_sites <- getSitesInFeature(
      cond_sites,
      gene_references$ref_genes,
      colnam = "in_gene",
      feature.colnam = "name2"
    )
    cond_sites <- getNearestFeature(
      cond_sites,
      gene_references$ref_genes,
      colnam = "nearest_gene",
      feature.colnam = "name2"
    )
    
    #' Add marks ("*" for in_gene, "~" for onco gene, and "!" for bad_actors)
    cond_sites$gene_id_wo_mark <- ifelse(
      cond_sites$in_gene == "FALSE",
      cond_sites$nearest_gene,
      cond_sites$in_gene
    )
    
    cond_sites$gene_mark <- ifelse(cond_sites$in_gene == "FALSE", "", "*")
    
    cond_sites$gene_mark <- ifelse(
      cond_sites$gene_id_wo_mark %in% gene_references$onco_genes,
      paste0(cond_sites$gene_mark, "~"),
      cond_sites$gene_mark
    )
    
    cond_sites$gene_mark <- ifelse(
      cond_sites$gene_id_wo_mark %in% gene_references$bad_actors,
      paste0(cond_sites$gene_mark, "!"),
      cond_sites$gene_mark
    )
    
    cond_sites$gene_id <- paste0(
      sapply(strsplit(cond_sites$gene_id_wo_mark, ","), "[[", 1), 
      " ", 
      cond_sites$gene_mark)
  }
  
  list(
    "all_sites" = process_data,
    "sites_final" = cond_sites[order(cond_sites$estAbund, decreasing = TRUE)]
  )
}
