# Helper: Parse VCF anchor-base notation to HGVS
vcf_to_hgvs <- function(chrom, pos, ref, alt) {
  # Aggressively strip all hidden whitespace from vectors
  chrom <- gsub("\\s+", "", as.character(chrom))
  pos   <- as.numeric(gsub("\\s+", "", as.character(pos)))
  ref   <- gsub("\\s+", "", as.character(ref))
  alt   <- gsub("\\s+", "", as.character(alt))
  
  mapply(function(c, p, r, a) {
    if (is.na(p) || is.na(r) || is.na(a)) return(NA_character_)
    
    # Format chromosome (chrM -> m, chr1 -> 1)
    c_fmt <- sub("^chr", "", c, ignore.case = TRUE)
    c_fmt <- ifelse(toupper(c_fmt) == "M", "m", c_fmt)
    prefix <- paste0(c_fmt, ".")
    
    len_r <- nchar(r)
    len_a <- nchar(a)
    
    if (len_r == 1 && len_a == 1) {
      # SNV
      return(paste0(prefix, p, r, ">", a))
    } else if (len_r > len_a) {
      # Deletion
      del_seq <- substr(r, len_a + 1, len_r)
      start_pos <- p + len_a
      end_pos <- p + len_r - 1
      if (start_pos == end_pos) {
        return(paste0(prefix, start_pos, "del", del_seq))
      } else {
        return(paste0(prefix, start_pos, "_", end_pos, "del", del_seq))
      }
    } else if (len_r < len_a) {
      # Insertion
      ins_seq <- substr(a, len_r + 1, len_a)
      return(paste0(prefix, p, "_", p + 1, "ins", ins_seq))
    } else {
      # MNV or complex substitution
      return(paste0(prefix, p, "_", p + len_r - 1, "delins", a))
    }
  }, chrom, pos, ref, alt, USE.NAMES = FALSE)
}

# Helper: Parse Annotation (Annovar) notation to HGVS
anno_to_hgvs <- function(chrom, start, end, ref, alt) {
  # Aggressively strip all hidden whitespace from vectors
  chrom <- gsub("\\s+", "", as.character(chrom))
  start <- as.numeric(gsub("\\s+", "", as.character(start)))
  end   <- as.numeric(gsub("\\s+", "", as.character(end)))
  ref   <- gsub("\\s+", "", as.character(ref))
  alt   <- gsub("\\s+", "", as.character(alt))
  
  mapply(function(c, s, e, r, a) {
    if (is.na(s)) return(NA_character_)
    
    c_fmt <- sub("^chr", "", c, ignore.case = TRUE)
    c_fmt <- ifelse(toupper(c_fmt) == "M", "m", c_fmt)
    prefix <- paste0(c_fmt, ".")
    
    if (r == "-" && a != "-") {
      # Insertion
      return(paste0(prefix, s, "_", s + 1, "ins", a))
    } else if (a == "-" && r != "-") {
      # Deletion
      if (s == e) {
        return(paste0(prefix, s, "del", r))
      } else {
        return(paste0(prefix, s, "_", e, "del", r))
      }
    } else if (nchar(r) == 1 && nchar(a) == 1) {
      # SNV
      return(paste0(prefix, s, r, ">", a))
    } else {
      # Substitution
      if (s == e) {
        return(paste0(prefix, s, "delins", a))
      } else {
        return(paste0(prefix, s, "_", e, "delins", a))
      }
    }
  }, chrom, start, end, ref, alt, USE.NAMES = FALSE)
}

# Main Pipeline
run_etl_pipeline <- function(samples_dir, sample_map_dir, output_dir){
  suppressPackageStartupMessages({
    library(tidyverse)
    library(vcfR)
    library(arrow)
  })
  
  options(stringsAsFactors = FALSE, scipen = 999)
  
  tryCatch({
    if (file.exists(sample_map_dir)) {
      sample_map <- readxl::read_excel(sample_map_dir)
      
      sample_map <- sample_map %>%
        mutate(GROUP = sapply(str_split(NAME, "-"), `[`, 1))
    }
  }, error = function(e){
    stop("Sample map file not found.")
  }) 
  
  all_variants <- list()
  all_annotations <- list()
  
  samples <- list.dirs(samples_dir)[-1]
  
  for (i in seq_along(samples)) {
    
    folder <- samples[i]
    sample_id <- basename(folder)
    message(sprintf("[%d/%d] Processing %s", i, length(samples), sample_id))
    
    vcf_file <- list.files(folder, pattern = "\\.(vcf|vcf\\.gz)?$", full.names = TRUE)
    if (!length(vcf_file)) {
      warning("No .vcf or .vcf.gz files found in: ", sample_id, call. = FALSE)
    }
    
    anno_file <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
    if (!length(anno_file)) {
      warning("No annotation file found in: ", sample_id, call. = FALSE)
    }
    
    tryCatch({
      vcf <- read.vcfR(vcf_file, verbose = FALSE)
      anno_data <- read_csv(anno_file, col_types = cols(.default = col_character()))
    }, error = function(e){
      warning("No file available for:", sample_id, call. = F)
    })
    
    # 2. Convert VCF to HGVS
    vcf_fix <- as.data.frame(vcf@fix) %>%
      mutate(VariantID = vcf_to_hgvs(CHROM, POS, REF, ALT))
    
    gt_data <- extract_gt_tidy(vcf, gt_column_prepend = "", 
                               VariantKey = vcf_fix$VariantID, verbose = F)
    
    variant_info <- vcf_fix %>%
      left_join(gt_data, by = c("VariantID" = "VariantKey"))
    
    all_variants[[sample_id]] <- variant_info
    
    # 3. Convert Annotation to HGVS (Safely handles missing 'end' columns)
    anno_data <- anno_data %>% 
      mutate(
        VariantID = anno_to_hgvs(
          chrom = chrom, 
          start = start, 
          end = if("end" %in% names(.)) end else start, 
          ref = REF, 
          alt = ALT
        ),
        Sample_ID = sample_id
      )
    
    all_annotations[[sample_id]] <- anno_data
  }
  
  # 4. Combine and deduplicate Variants table
  final_variants <- bind_rows(all_variants) %>% distinct(VariantID, .keep_all = TRUE)
  final_annotations <- bind_rows(all_annotations)
  
  # 5. Persistent backup creation
  write_dataset(final_variants, path = file.path(output_dir, "variants"), format = "parquet")
  write_dataset(final_annotations, path = file.path(output_dir, "annotations"), format = "parquet")
  write_dataset(sample_map, path = file.path(output_dir, "metadata"), format = "parquet")
  
  return(TRUE)
}