#########################################
# Blood Plasma Marker Selection Pipeline
# BioMed Innovation Lab
# University of Guelph, Canada
# July, 2021
#########################################

################################
# System Requirement
# OS: Linux or BSD (e.g., Mac) with zcat, tail, awk and zgrep commands available.
# Access to Internet is required
################################


# Load the required packages
bpmsp.libs = c('ggrepel', 'data.table', 'stringr', 'dplyr', 'cgdsr')
for (lib in bpmsp.libs){
  if( !is.element(lib, .packages(all.available = TRUE)) ) {
    install.packages(lib)
  }
  library(lib,character.only = TRUE)
}


# Function to load blood plasma marker data 
bpmsp.load_blood_plasma_markers <- function() {
  # Source: https://www.proteinatlas.org/humanproteome/blood/proteins+detected+by+immunoassay
  proteins_detected_by_immunoassay = read.delim('human_protein_atlas/proteins_detected_by_immunoassay.tsv', sep='\t', header = TRUE, stringsAsFactors = FALSE)
  
  # Source: https://www.proteinatlas.org/humanproteome/blood/proteins+detected+in+ms
  proteins_detected_by_ms = read.delim('human_protein_atlas/proteins_detected_in_ms.tsv', sep='\t', header = TRUE, stringsAsFactors = FALSE)
  
  # Source: https://www.proteinatlas.org/humanproteome/blood/proteins+detected+by+pea
  proteins_detected_by_pea = read.delim('human_protein_atlas/proteins_detected_by_pea.tsv', sep='\t', header = TRUE, stringsAsFactors = FALSE)
  
  # Source: https://www.proteinatlas.org/search/protein_class:Plasma+proteins, sourced from Plasma Protein Database as indicated 
  # in https://www.proteinatlas.org/humanproteome/proteinclasses
  proteins_from_plasma_protein_db = read.delim('human_protein_atlas/protein_class_Plasma.tsv', sep='\t', header = TRUE, stringsAsFactors = FALSE)
  
  plasma_protein_sets = list(
    Immunoassay = proteins_detected_by_immunoassay$Gene,
    MS = proteins_detected_by_ms$Gene,
    PEA = proteins_detected_by_pea$Gene,
    Plasma_Protein_DB = proteins_from_plasma_protein_db$Gene
  )
  
  # Save the data in a global variable
  bpmsp.all_plasma_protein_genes <<- Reduce(union, plasma_protein_sets)
  # write(paste(bpmsp.all_plasma_protein_genes, collapse="\n"), file="plasma_genes.txt")
  message("TCGA/GTEx data saved in global variable bpmsp.all_plasma_protein_genes")
}


# Function to load the re-aligned TCGA and GTEX data from UCSC Xena (recomputed by Toil)
bpmsp.load_recomputed_tcga_gtex_data <- function() {
  
  # First, download the Ensembl ID to gene names mapping from GTEX, if not already
  if (!file.exists('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz')) {
    message("Downloading GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz from GTEx, this will take a while...")
    download.file('https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz', 
                  'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz')
  }

  # Then, run a shell script to extract the Ensembl IDs and the gene names, if not already
  if (!file.exists('gtex_ensembl_id_gene_names.txt')) {
    message("Extracting Ensembl IDs...")
    system("zcat GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz | awk '{print $1,$2}' > gtex_ensembl_id_gene_names.txt")
  }
  gtex_ensembl_id_gene_names = fread("gtex_ensembl_id_gene_names.txt", header=FALSE, sep = " ")
  
  # Find the unique ones
  unique_gtex_genes = unique(gtex_ensembl_id_gene_names[,2])$V2 # Get the gene names
  
  # Short list the markers to include those blood plasma markers only
  bpmsp.gtex_plasma_genes <<- intersect(unique_gtex_genes, bpmsp.all_plasma_protein_genes) 
  gtex_plasma_genes_with_ensembl_ids = gtex_ensembl_id_gene_names[gtex_ensembl_id_gene_names$V2 %in% bpmsp.gtex_plasma_genes, ]
  
  # Remove irrelevant records
  gtex_plasma_genes_with_ensembl_ids = gtex_plasma_genes_with_ensembl_ids[!grepl("(.*_PAR_Y)|ENSG00000228741.2", gtex_plasma_genes_with_ensembl_ids$V1), ]
  gtex_plasma_genes_with_ensembl_ids = gtex_plasma_genes_with_ensembl_ids[match(bpmsp.gtex_plasma_genes, gtex_plasma_genes_with_ensembl_ids$V2), ] #Keeping the same order
  
  # Save to a file for later use
  gtex_plasma_genes = gtex_plasma_genes_with_ensembl_ids$V2
  names(gtex_plasma_genes) = gtex_plasma_genes_with_ensembl_ids$V1
  bpmsp.gtex_plasma_genes <<- gtex_plasma_genes # Expose as a global variable
  write(paste(names(bpmsp.gtex_plasma_genes), collapse="\n"), file=str_glue("gtex_plasma_gene_ensembl_ids.txt"))
  
  # Now, download the Toil-recomputed dataset, if not already
  if (!file.exists('TcgaTargetGtex_rsem_gene_tpm.gz')) {
    message("Downloading 'TcgaTargetGtex_rsem_gene_tpm.gz' from USSC Xena, this will take a while...")
    download.file('https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/TcgaTargetGtex_rsem_gene_tpm.gz', 'TcgaTargetGtex_rsem_gene_tpm.gz')
  }
  
  # Then, run shell scripts to extract all the data matching the in scope markers, if not already
  if (!file.exists('in_scope_tcga_gtex_rsem_gene_tpm.tsv')) {
    message("Extracting in scope TCGA/GTEx data, this will take a while...")
    gc()
    system("zcat TcgaTargetGtex_rsem_gene_tpm.gz 2>/dev/null | head -n1 > in_scope_tcga_gtex_rsem_gene_tpm.tsv")
    system("zgrep -F -f gtex_plasma_gene_ensembl_ids.txt TcgaTargetGtex_rsem_gene_tpm.gz >> in_scope_tcga_gtex_rsem_gene_tpm.tsv")
  }
  
  # Save it in a global variable
  bpmsp.short_listed_tcga_gtex_rsem_gene_tpm <<- fread(str_glue("in_scope_tcga_gtex_rsem_gene_tpm.tsv"), header=TRUE, sep = "\t")
  message("TCGA/GTEx data saved in global variable bpmsp.short_listed_tcga_gtex_rsem_gene_tpm")
}

# Function to get cancer case study TCGA barcodes by the cancer name, and optionally a cancer_subtype_code 
# (i.e., case study prefix such as luad for adenocarcinoma)
bpmsp.get_cancer_study_barcodes <- function(cancer_name, cancer_subtype_code) {
  mycgds = CGDS("http://www.cbioportal.org/")
  cgds_cancer_studies = getCancerStudies(mycgds) # cgdsr function
  
  # We need to find the "TCGA" study. For now, use the study with "tcga_pan_can_atlas_2018" in its name
  my_cancer_studies = cgds_cancer_studies[grep(str_glue('^{cancer_name}'), cgds_cancer_studies$name, ignore.case = T), ]
  case_study_name_pattern = '.*_tcga_pan_can_atlas_2018'
  if (!missing(cancer_subtype_code)) {
    case_study_name_pattern = str_glue('{cancer_subtype_code}.*_tcga_pan_can_atlas_2018')
  }
  
  cancer_study_id = my_cancer_studies[grep(case_study_name_pattern, my_cancer_studies$cancer_study_id, ignore.case = T), 1]
  if(length(cancer_study_id) != 1) {
    stop("Please explore and choose the right case study. To explore more cases, run View(my_cancer_studies)")
  }
  
  case_lists = getCaseLists(mycgds, cancer_study_id)
  case_name_pattern = '.*_tcga_pan_can_atlas_2018_rna_seq_v2_mrna'
  if (!missing(cancer_subtype_code)) {
    case_name_pattern = str_glue('{cancer_subtype_code}.*_tcga_pan_can_atlas_2018_rna_seq_v2_mrna')
  }
  
  
  case_id=case_lists[grep(case_name_pattern, case_lists$case_list_id, ignore.case=T), 1]
  if(length(case_id) != 1) {
    stop("Please explore and choose the right case. To explore more cases, run View(case_lists)")
  }
  
  genetic_profiles = getGeneticProfiles(mycgds, cancer_study_id)
  # Looking for the one with rna_seq_v2_mran_median_Zscore as it has less data to download compared to the seq_v2_mrna one
  genetic_profile_name_pattern = '.*_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_Zscores'
  if (!missing(cancer_subtype_code)) {
    genetic_profile_name_pattern = str_glue('{cancer_subtype_code}.*_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_Zscores')
  }
  
  genetic_profile_id = genetic_profiles[grep(genetic_profile_name_pattern, genetic_profiles$genetic_profile_id, ignore.case=T), 1]
  
  # Assuming the first 100 gene records should be enough to uncover the TCGA barcodes
  tcga_barcodes = rownames(getProfileData(mycgds, bpmsp.gtex_plasma_genes[c(1:100)], genetic_profile_id, case_id))
  return(tcga_barcodes)
}

# Function to select the TCGA/GTEx data the barcodes
bpmsp.get_tcga_gtex_data_by_barcodes <- function(tcga_barcodes) {
  tcga_barcodes_reformatted = gsub('.', '-', tcga_barcodes, fixed=T) # Somehow the format is slightly different
  tcga_gtex_tpm = bpmsp.short_listed_tcga_gtex_rsem_gene_tpm  %>% select(matches("GTEX-.*|sample") | any_of(tcga_barcodes_reformatted))
  tcga_barcodes_ensembl_id = tcga_gtex_tpm$sample
  tcga_gtex_tpm = tcga_gtex_tpm[, -c(1)]
  rownames(tcga_gtex_tpm) = bpmsp.gtex_plasma_genes[tcga_barcodes_ensembl_id]
  return(tcga_gtex_tpm)
}

# Function to get the top 5% markers with the greatest expression level difference between tumor and normal samples
bpmsp.short_list_markers_by_expression_level_difference <- function(tcga_gtex_tpm) {
  tcga_gtex_tpm_normal_tpm_log2 = rowMeans(tcga_gtex_tpm %>% select(matches("GTEX-.*")), na.rm = T)
  tcga_gtex_tpm_tumo_tpm_log2 = rowMeans(tcga_gtex_tpm %>% select(matches("TCGA-.*")), na.rm = T)
  
  tpm_data = data.frame(x = tcga_gtex_tpm_tumo_tpm_log2, y = tcga_gtex_tpm_normal_tpm_log2, z = tcga_gtex_tpm_tumo_tpm_log2 - tcga_gtex_tpm_normal_tpm_log2)
  rownames(tpm_data) = rownames(tcga_gtex_tpm)
  
  tpm_diff_shortlisted=rownames(tpm_data[tpm_data$z > quantile(tpm_data$z, 0.95), ])
  return(tpm_diff_shortlisted)
}

# Function to download IHC data
bpmsp.load_ihc_data <-function() {
  # Download normal and pathology data from HPA, if not already
  if (!file.exists('normal_tissue.tsv')) {
    message("Downloading 'normal_tissue.tsv.zip' from HPA")
    download.file('https://www.proteinatlas.org/download/normal_tissue.tsv.zip', 'normal_tissue.tsv.zip')
    unzip('normal_tissue.tsv.zip')
  }
  
  if (!file.exists('pathology.tsv')) {
    message("Downloading 'pathology.tsv.zip' from HPA")
    download.file('https://www.proteinatlas.org/download/pathology.tsv.zip', 'pathology.tsv.zip')
    unzip('pathology.tsv.zip')
  }
  
  bpmsp.normal_tissue_ihc <<- read.delim('normal_tissue.tsv', sep='\t', header = TRUE, stringsAsFactors = FALSE)
  bpmsp.tumor_tissue_ihc <<- read.delim('pathology.tsv', sep='\t', header = TRUE, stringsAsFactors = FALSE)
  message("IHC data are loaded into the global variables bpmsp.normal_tissue_ihc and bpmsp.tumor_tissue_ihc")
}

bpmsp.short_list_markers_by_IHC_score_z_scores <- function(cancer_name) {
  # Select pathology data by cancer name
  cancer_tissue_ihc = bpmsp.tumor_tissue_ihc[grep(str_glue('^{cancer_name}'), bpmsp.tumor_tissue_ihc$Cancer), ]
  
  # Assign scores to High, Medium, Low, Not Detected
  scores = diag(c(3, 2, 1, 0))
  cancer_tissue_ihc$meanScores = rowMeans(as.matrix(cancer_tissue_ihc[, c("High", "Medium", "Low", "Not.detected")]) %*% scores)
  
  # Remove NA and duplicates
  sanitized_cancer_tissue_ihc = aggregate(meanScores ~ Gene.name, cancer_tissue_ihc, mean)
  
  # Limit the scope to those in blood plasma marker list
  genes_selected_for_analysis = intersect(bpmsp.all_plasma_protein_genes, sanitized_cancer_tissue_ihc$Gene.name)
  
  # Further limit the scope to those in normal IHC data
  genes_in_normal_tissue_ihc = unique(bpmsp.normal_tissue_ihc$Gene.name)
  genes_selected_for_analysis = intersect(genes_selected_for_analysis, genes_in_normal_tissue_ihc)
  
  # Now compute the Z-score of tumor tissue IHC scores
  tumor_tissue_ihc_z_scores = as.data.frame(scale(sanitized_cancer_tissue_ihc$meanScores))$V1
  names(tumor_tissue_ihc_z_scores) = sanitized_cancer_tissue_ihc$Gene.name
  selected_tumor_ihc_z_scores = tumor_tissue_ihc_z_scores[genes_selected_for_analysis]
  
  
  # Normal tissue IHC data format is a bit different, so need to count the occurence to compute the scores
  selected_normal_tissue_ihc_records = bpmsp.normal_tissue_ihc[bpmsp.normal_tissue_ihc$Gene.name %in% genes_selected_for_analysis,]
  selected_normal_tissue_ihc = table(selected_normal_tissue_ihc_records$Gene.name, selected_normal_tissue_ihc_records$Level)
  selected_normal_tissue_ihc.meanScores = rowMeans(as.matrix(selected_normal_tissue_ihc[, c("High", "Medium", "Low", "Not detected")]) %*% scores)
  
  # Now compute the Z-score of normal tissue IHC scores
  normal_tissue_ihc_z_scores =  as.data.frame(scale(selected_normal_tissue_ihc.meanScores))$V1
  names(normal_tissue_ihc_z_scores) = names(selected_normal_tissue_ihc.meanScores)
  selected_normal_tissue_ihc_z_scores = normal_tissue_ihc_z_scores[genes_selected_for_analysis]
  
  ihc_data = data.frame(x = selected_tumor_ihc_z_scores, y = selected_normal_tissue_ihc_z_scores, z = selected_tumor_ihc_z_scores - selected_normal_tissue_ihc_z_scores)
  
  ihc_z_score_diff_shortlisted=rownames(ihc_data[ihc_data$z > quantile(ihc_data$z, 0.95), ])
  return(ihc_z_score_diff_shortlisted)
}
