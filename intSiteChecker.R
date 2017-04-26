#' Check specific integration site data with a commandline script. The script
#' connects to an INSPIIRED database, acquires all integration sites for the 
#' single specimen (GTSP????), refines breakpoints, standardizes integration
#' site positions, annotates against current onco-related gene and bad actor 
#' lists, and then returns a small table with the number of most abundant clones
#' and respecive information. Parameters and databases can be adjusted in the
#' config.yml file. For efficient processing, only unique integration sites and 
#' not multihits are analyzed.
options(stringsAsFactors = FALSE)
suppressMessages(library("argparse"))
suppressMessages(library("yaml"))
suppressMessages(library("pander"))

code_dir <- dirname(
  sub("--file=", "", grep("--file=", commandArgs(trailingOnly=FALSE), value=T)))

#' Set up and gather command line arguments
config <- suppressWarnings(
  yaml.load_file(file.path(code_dir, "config.yml")))
descript <- suppressWarnings(
  yaml.load_file(file.path(code_dir, "descriptions.yml")))

parser <- ArgumentParser(description = descript$program_short_description)
parser$add_argument(
  "specimen", type = "character", nargs = "+", 
  help = descript$specimen)
parser$add_argument(
  "-p", "--position_ids", type = "character", nargs = "+",
  help = descript$position_ids)
parser$add_argument(
  "-e", "--position_ids_exclusively", action = "store_true", 
  help = descript$position_ids_exclusively)
parser$add_argument(
  "-a", "--abundance_method", type = "character", 
  default = config$abundance_method, help = descript$abundance_method)
parser$add_argument(
  "-n", "--num_of_clones", type = "integer", default = config$num_of_top_clones,
  help = descript$num_of_clones)
parser$add_argument(
  "-f", "--fuzzy_window", type = "integer", default = config$fuzzy_window,
  help = descript$fuzzy_window)
parser$add_argument(
  "-r", "--ref_genome", type = "character", default = config$ref_genome,
  help = descript$ref_genome)
parser$add_argument(
  "-c", "--cores", type = "integer", default = config$cores,
  help = descript$cores)
parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", help = descript$output)

#' Dev command
#cmdline <- c("GTSP0518", "GTSP0853", "GTSP1322", "-p", "chr12-4270498", "-e")
#args <- parser$parse_args(cmdline)

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))
args$specimen <- unique(unlist(strsplit(args$specimen, " ")))
args$specimen_database <- ifelse(
  config$specimen_db_type == "mysql",
  config$specimen_db_group, config$specimen_db_name)
args$intsite_database <- ifelse(
  config$intsite_db_type == "mysql",
  config$intsite_db_group, config$intsite_db_name)

if(length(args$position_ids) > 0){
  args$position_ids <- unique(unlist(strsplit(args$position_ids, " ")))}
args <- args[c("specimen", "position_ids", "position_ids_exclusively", 
               "num_of_clones", "ref_genome", "abundance_method", 
               "fuzzy_window", "cores", "intsite_database", 
               "specimen_database","output")]

if(args$cores > 0){
  suppressMessages(library("parallel"))
  max_cores <- detectCores()
  if(args$cores > max_cores){ args$cores <- max_cores }
  source(file.path(
    code_dir, "support_scripts", "parallel_intsite_processing.R"))
}else{
  source(file.path(
    code_dir, "support_scripts", "serial_intsite_processing.R"))
}

#' Print inputs to the console
if(config$print_params){
  input_table <- data.frame(
    "Variables" = paste0(names(args), " :"), 
    "Values" = sapply(1:length(args), function(i){
      paste(args[[i]], collapse = ", ")
    }))
}else{
  input_table <- data.frame(
    "Variables" = paste0(names(args), " :")[1:6], 
    "Values" = sapply(1:6, function(i){
      paste(args[[i]], collapse = ", ")
    }))
}

if(config$print_inputs){
  if(length(args$position_ids) > 0){
    pandoc.title("Analytical Inputs")
    pandoc.table(input_table, justify = c("left", "center"))
  }else{
    exclude <- match(
      c("position_ids :", "position_ids_exclusively :"), 
      input_table$Variables)
    input_table <- input_table[-exclude,]
    row.names(input_table) <- NULL
    pandoc.title("Analytical Inputs")
    pandoc.table(input_table, justify = c("left", "center"))
    rm(exclude)
  }
}

#' Load additional dependencies
add_dependencies <- c("GenomicRanges", "dplyr", "hiAnnotator", 
                      "sonicLength", "gintools", "DBI")

null <- suppressMessages(
  sapply(add_dependencies, library, character.only = TRUE))

rm(null)

#' Download specimen information from specimen database
if(config$specimen_db_type == "mysql"){
  pack <- suppressMessages(require("RMySQL"))
  if(!pack) stop("Could not load MySQL package.")
  junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
  dbConn <- dbConnect(MySQL(), group = args$specimen_database)
  stopifnot(dbGetQuery(dbConn, "SELECT 1") == 1)
  rm(pack)
}else if(config$specimen_db_type == "sqlite"){
  pack <- suppressMessages(require("RSQLite"))
  if(!pack) stop("Could not load SQLite package.")
  #junk <- sapply(dbListConnections(SQLite()), dbDisconnect)
  dbConn <- dbConnect(SQLite(), dbname = args$specimen_database)
  stopifnot(dbGetQuery(dbConn, "SELECT 1") == 1)
  rm(pack)
}else{
  stop("Could not establish a specimen database connection. 
       Check config for database type ('mysql' or 'sqlite').")
}

meta_cols <- c("Trial", "Patient", "CellType",
               "Timepoint", "VCN", "SpecimenAccNum")
query_selection <- paste(
  "SELECT",
  paste(meta_cols, collapse = ", "),
  "FROM gtsp",
  sep = " ")
specimen_cond <- paste(args$specimen, collapse = "', '")
query_condition <- paste0("WHERE SpecimenAccNum IN ('", specimen_cond, "')")
query <- paste(query_selection, query_condition, sep = " ")
message(query)
specimen_data <- dbGetQuery(dbConn, query)
names(specimen_data) <- tolower(names(specimen_data))

if(config$specimen_db_type == "mysql"){
  junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
  detach("package:RMySQL")
}else if(config$specimen_db_type == "sqlite"){
  dbDisconnect(dbConn)
  detach("package:RSQLite")
}
rm(query, query_condition, query_selection, meta_cols, dbConn)

#' Print specimen information to the console
if(config$print_summary){
  pandoc.title("Specimen Summary")
  pandoc.table(
    t(specimen_data), justify = c("left", rep("center", length(args$specimen))))
}

#' Download integration sites from INSPIIRED database
if(config$intsite_db_type == "mysql"){
  pack <- suppressMessages(require("RMySQL"))
  if(!pack) stop("Could not load MySQL package.")
  dbConn <- dbConnect(MySQL(), group = args$intsite_database)
  stopifnot(dbGetQuery(dbConn, "SELECT 1") == 1)
  rm(pack)
}else if(config$insite_db_type == "sqlite"){
  pack <- suppressMessages(require("RSQLite"))
  if(!pack) stop("Could not load SQLite package.")
  dbConn <- dbConnect(SQLite(), dbname = args$intsite_database)
  stopifnot(dbGetQuery(dbConn, "SELECT 1") == 1)
  rm(pack)
}else{
  stop("Could not establish a intsite database connection. 
       Check config for database type ('mysql' or 'sqlite').")
}

query_uniq <- sprintf("SELECT * FROM samples
      JOIN sites ON samples.sampleID = sites.sampleID
      JOIN pcrbreakpoints ON pcrbreakpoints.siteID = sites.siteID
      WHERE LEFT(samples.sampleName, INSTR(samples.sampleName, '-')-1)
      IN ('%1$s')",
      specimen_cond)

message(query_uniq)
uniq_sites_data <- dbGetQuery(dbConn, query_uniq)

if(config$intsite_db_type == "mysql"){
  junk <- sapply(dbListConnections(MySQL()), dbDisconnect)
  detach("package:RMySQL")
}else if(config$intsite_db_type == "sqlite"){
  dbDisconnect(dbConn)
  detach("package:RSQLite")
}
rm(query, query_condition, query_selection, meta_cols, dbConn)

#' Gather reference materials: onco-related gene list, bad actors, ref_genes
ref_genes <- suppressWarnings(try(
  readRDS(config$ref_genes_path), silent = TRUE))
if(class(ref_genes) != "GRanges" & config$print_deviations){
  pandoc.emphasis(descript$refGenesNotLoaded)}

onco_genes_data <- read.delim(
  config$oncogene_list,
  header = TRUE,
  sep = "\t")
onco_genes <- unique(onco_genes_data[,"symbol"])

bad_actors_data <- read.delim(
  config$bad_actor_list,
  header = TRUE,
  sep = "\t")
bad_actors <- unique(bad_actors_data[,"symbol"])

gene_references <- list(
  "ref_genes" = ref_genes, "onco_genes" = onco_genes, "bad_actors" = bad_actors)

#' If exclusively processing integration sites for specific position ids, then
#' the uniq_sites_data needs to be filtered for sites near the input ids. 
if(length(args$position_ids) > 0){
  chr_pos <- strsplit(args$position_ids, "[+-]", perl = TRUE)
  strands <- ifelse(grepl("+", args$position_ids, fixed = TRUE), "+", "-")
  posids <- data.frame(
    "chr" = sapply(chr_pos, "[[", 1),
    "strand" = strands,
    "position" = as.numeric(sapply(chr_pos, "[[", 2)))
  posids$min_position <- posids$position - args$fuzzy_window
  posids$max_position <- posids$position + args$fuzzy_window
  if(args$position_ids_exclusively){
    select_sites <- do.call(
      rbind, 
      lapply(1:nrow(posids), function(i, data, posids){
        data[
          data$chr == posids[i, "chr"] & 
          data$strand == posids[i, "strand"] & 
          data$position >= posids[i, "min_position"] & 
          data$position <= posids[i, "max_position"],]
    }, data = uniq_sites_data, posids = posids))
    
    uniq_sites_data <- unique(select_sites)
  }
}

#' Process integration sites for analysis. This step includes refining 
#' breakpoints, standardizing integration site positions, annotating against
#' the ref_genes dataset (if available), and developing a gene id for the site.
sites <- process_intsite_data(
  uniq_sites_data, specimen_data, args, config, gene_references)

#' Print summary table(s) of the top clones or position ids
summary_table <- GenomicRanges::as.data.frame(sites$sites_final) %>%
  select(-seqnames, -start, -end, -width, -strand, -in_gene, -in_geneOrt, 
         -nearest_geneDist, -nearest_gene, -nearest_geneOrt, -gene_id_wo_mark,
         -gene_mark) %>%
  group_by(specimen) %>%
  top_n(n = ifelse(args$num_of_clones == 0, nrow(.), args$num_of_clones), 
        estAbund) %>%
  ungroup()
  
if(!args$position_ids_exclusively){
  summary_table <- mutate(
    summary_table, relAbund = round(relAbund, digits = 3))
}

if(class(ref_genes) == "try-error"){
  summary_table <- mutate(
    summary_table, gene_id = rep(NA, n()))
}

if(length(args$position_ids) > 0){
  posid_table <- do.call(rbind, lapply(
    1:nrow(posids), function(i, data, posids){
      data_chr_pos <- strsplit(data$posID, "[+-]", perl = TRUE)
      data_strands <- ifelse(grepl("+", data$posID, fixed = TRUE), "+", "-")
      data_posids <- data.frame(
        "chr" = sapply(data_chr_pos, "[[", 1),
        "strand" = data_strands,
        "position" = sapply(data_chr_pos, "[[", 2)
      )
      ids_to_keep <- which(data_posids$chr == posids[i, "chr"] &
                           data_posids$strand == posids[i, "strand"] &
                           data_posids$position >= posids[i, "min_position"] &
                           data_posids$position <= posids[i, "max_position"])
      data[ids_to_keep,]
    }, data = summary_table, posids = posids))
  
  summary_table <- distinct(bind_rows(posid_table, summary_table))
}
  
summary_table <- split(summary_table, summary_table$specimen)
  
if(!args$position_ids_exclusively){
  tables <- lapply(summary_table, function(sum_tab){
    title_cols <- c("patient", "timepoint", "celltype", "specimen")
    table_cols <- c("gene_id", "posID", "estAbund", "relAbund")
    col_nums <- match(title_cols, names(sum_tab))
    pandoc.title(paste0(paste(sum_tab[1, title_cols], collapse = " - "), "\n", 
                        "Most Abundant Clones within Specimen ", 
                        unique(sum_tab$specimen)))
    pandoc.table(arrange(sum_tab[, table_cols], desc(estAbund)))
  })
}else{
  tables <- lapply(summary_table, function(sum_tab){
    title_cols <- c("patient", "timepoint", "celltype", "specimen")
    table_cols <- c("gene_id", "posID", "estAbund")
    pandoc.title(paste0(paste(sum_tab[1, title_cols], collapse = " - "), "\n", 
                        "Most Abundant Clones within Specimen ", 
                        unique(sum_tab$specimen)))
    pandoc.table(arrange(sum_tab[, table_cols], desc(estAbund)))
  })
}

if(length(args$output) > 0){
  sites$summary_table <- bind_rows(summary_table)
  saveRDS(sites, file = paste0(args$output, ".rds"))
}
q()
