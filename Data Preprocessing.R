### -------- DATA PRE-PROCESSING -------------- ###

# Library load
load_libraries <- function() {
  library(tidyverse) # for data manipulation
  library(missForest) # for data imputation for UKB data
  library(caret) # for correlation filteration
}

# Set working directory
set_working_directory <- function(directory_path) {
  setwd(dir = directory_path)
}

# Function to read data for Cohort 1 (UK Biobank)
read_cohort1_data <- function(features_file, diagnosis_file) {
  BBdata_metadata <- read.csv(features_file, header = TRUE, row.names = NULL)
  BBdiagnosis_metadata <- read.csv(diagnosis_file, header = TRUE, row.names = NULL)
  
  return(list(BBdata_metadata = BBdata_metadata, BBdiagnosis_metadata = BBdiagnosis_metadata))
}

# Function to preprocess Cohort 1 (UK Biobank) data
preprocess_cohort1_data <- function(BBdata_metadata, BBdiagnosis_metadata) {
  
  # Cleaning and Preprocessing
  # ---------------------------------
  
  # Subsetting the metabolite section
  my_string <- "Metabolite"
  my_cols <- grep(my_string, colnames(BBdata_metadata), value = TRUE)
  BBmetabolomics <- BBdata_metadata[, my_cols]
  rownames(BBmetabolomics) <- BBdata_metadata$f.eid
  
  # Formatting the names of the metabolites
  my_string <- ".0.0.Metabolites"
  new_colnames <- gsub(my_string, "", colnames(BBmetabolomics))
  colnames(BBmetabolomics) <- new_colnames
  
  # Cleaning unwanted objects
  rm(list = c("my_string", "my_cols", "new_colnames"))
  
  # Subset to keep only metabolite samples for UK Biobank
  max_null_percentage <- 50  # The metabolite samples are over 90% filled for all metabolites
  null_percentages <- rowMeans(is.na(BBmetabolomics)) * 100
  BBmetabolomics <- BBmetabolomics[null_percentages < max_null_percentage, ]
  
  # Formatting the Metadata
  BBmetadata <- BBdata_metadata[, c((ncol(BBdata_metadata) - 5):ncol(BBdata_metadata))]
  BBmetadata <- BBmetadata[BBmetadata$f.eid %in% rownames(BBmetabolomics), ]
  BBdiagnosis_metadata <- BBdiagnosis_metadata[BBdiagnosis_metadata$f.eid %in% rownames(BBmetabolomics), ]
  
  head(BBdiagnosis_metadata)
  count(distinct(BBdiagnosis_metadata, f.eid))  # How many are distinct?
  
  # Function to convert date format to a unified one of YYYY-MM-DD
  convert_date_format <- function(date) {
    if (grepl("^\\d{4}\\.\\d$", date)) {
      date <- str_replace(date, "\\.", "-0")
      date <- paste(date, "15", sep = "-")
    } else if (grepl("^\\d{4}\\.\\d{2}$", date)) {
      date <- str_replace(date, "\\.", "-")
      date <- paste(date, "15", sep = "-")
    } else if (grepl("^\\d{4}$", date)) {
      date <- paste(date, "06-15", sep = "-")
    }
    return(date)
  }
  
  # Applying the conversion function to the 'Date_value' column
  BBdiagnosis_metadata$Date_value <- sapply(BBdiagnosis_metadata$Date_value, convert_date_format)
  
  # Subset the dataframe to keep only the sample duplicates with the latest dates
  BBdiagnosis_metadata$Date_value <- as.Date(BBdiagnosis_metadata$Date_value)  # Convert to Date
  
  # Filtering to keep only the rows with the maximum Date_value within each f.eid group
  BBdiagnosis_metadata <- BBdiagnosis_metadata %>%
    group_by(f.eid) %>%
    arrange(desc(Date_value)) %>%
    slice(if (max(Date_value) %in% Date_value) 1 else 0)
  
  # Combining the diagnosis data with the metadata to form the new BBmetadata
  BBmetadata <- merge(BBmetadata, BBdiagnosis_metadata, by = 'f.eid', all = TRUE)
  rownames(BBmetadata) <- BBmetadata$f.eid
  BBmetadata <- select(BBmetadata, -f.eid)
  
  # Formatting BBmetadata columns
  BBmetadata$Ethnicity <- as.character(BBmetadata$Ethnicity)
  BBmetadata$Sex <- as.character(BBmetadata$Sex)
  BBmetadata$casecontrol <- as.character(BBmetadata$casecontrol)
  # change name of "casecontrol" to "Label" for MetaboAnalyst format
  names(BBmetadata)[names(BBmetadata) == "casecontrol"] <- "Label"
  
  # Formatting the BBdata_metadata
  rownames(BBdata_metadata) <- BBdata_metadata$f.eid
  BBdata_metadata <- BBdata_metadata[rownames(BBdata_metadata) %in% rownames(BBmetabolomics), ]
  
  # -----
  
  # Renaming the BBmetabolite names
  # -------------------------------
  
  name_mapping <- c(
    "Total.Cholesterol" = "Total-C",
    "Total.Cholesterol.Minus.HDL.C" = "non-HDL-C",
    "Remnant.Cholesterol..Non.HDL..Non.LDL..Cholesterol." = "Remnant-C",
    "VLDL.Cholesterol" = "VLDL-C",
    "Clinical.LDL.Cholesterol" = "Clinical LDL-C",
    "LDL.Cholesterol" = "LDL-C",
    "HDL.Cholesterol" = "HDL-C",
    "Total.Triglycerides" = "Total-TG",
    "Triglycerides.in.VLDL" = "VLDL-TG",
    "Triglycerides.in.LDL" = "LDL-TG",
    "Triglycerides.in.HDL" = "HDL-TG",
    "Total.Phospholipids.in.Lipoprotein.Particles" = "Total-PL",
    "Phospholipids.in.VLDL" = "VLDL-PL",
    "Phospholipids.in.LDL" = "LDL-PL",
    "Phospholipids.in.HDL" = "HDL-PL",
    "Total.Esterified.Cholesterol" = "Total-CE",
    "Cholesteryl.Esters.in.VLDL" = "VLDL-CE",
    "Cholesteryl.Esters.in.LDL" = "LDL-CE",
    "Cholesteryl.Esters.in.HDL" = "HDL-CE",
    "Total.Free.Cholesterol" = "Total-FC",
    "Free.Cholesterol.in.VLDL" = "VLDL-FC",
    "Free.Cholesterol.in.LDL" = "LDL-FC",
    "Free.Cholesterol.in.HDL" = "HDL-FC",
    "Total.Lipids.in.Lipoprotein.Particles" = "Total-L",
    "Total.Lipids.in.VLDL" = "VLDL-L",
    "Total.Lipids.in.LDL" = "LDL-L",
    "Total.Lipids.in.HDL" = "HDL-L",
    "Total.Concentration.of.Lipoprotein.Particles" = "Total-P",
    "Concentration.of.VLDL.Particles" = "VLDL-P",
    "Concentration.of.LDL.Particles" = "LDL-P",
    "Concentration.of.HDL.Particles" = "HDL-P",
    "Average.Diameter.for.VLDL.Particles" = "VLDL size",
    "Average.Diameter.for.LDL.Particles" = "LDL size",
    "Average.Diameter.for.HDL.Particles" = "HDL size",
    "Phosphoglycerides" = "Phosphoglyc",
    "Total.Cholines" = "Cholines",
    "Phosphatidylcholines" = "Phosphatidylc",
    "Sphingomyelins" = "Sphingomyelins",
    "Apolipoprotein.B" = "ApoB",
    "Apolipoprotein.A1" = "ApoA1",
    "Total.Fatty.Acids" = "Total-FA",
    "Degree.of.Unsaturation" = "Unsaturation",
    "Omega.3.Fatty.Acids" = "Omega-3",
    "Omega.6.Fatty.Acids" = "Omega-6",
    "Polyunsaturated.Fatty.Acids" = "PUFA",
    "Monounsaturated.Fatty.Acids" = "MUFA",
    "Saturated.Fatty.Acids" = "SFA",
    "Linoleic.Acid" = "LA",
    "Docosahexaenoic.Acid" = "DHA",
    "Alanine" = "Alanine",
    "Glutamine" = "Glutamine",
    "Glycine" = "Glycine",
    "Histidine" = "Histidine",
    "Total.Concentration.of.Branched.Chain.Amino.Acids..Leucine...Isoleucine...Valine." = "Total BCAA",
    "Isoleucine" = "Isoleucine",
    "Leucine" = "Leucine",
    "Valine" = "Valine",
    "Phenylalanine" = "Phenylalanine",
    "Tyrosine" = "Tyrosine",
    "Glucose" = "Glucose",
    "Lactate" = "Lactate",
    "Pyruvate" = "Pyruvate",
    "Citrate" = "Citrate",
    "X3.Hydroxybutyrate" = "bOHbutyrate",
    "Acetate" = "Acetate",
    "Acetoacetate" = "Acetoacetate",
    "Acetone" = "Acetone",
    "Creatinine" = "Creatinine",
    "Albumin" = "Albumin",
    "Glycoprotein.Acetyls" = "GlycA",
    "Concentration.of.Chylomicrons.and.Extremely.Large.VLDL.Particles" = "XXL-VLDL-P",
    "Total.Lipids.in.Chylomicrons.and.Extremely.Large.VLDL" = "XXL-VLDL-L",
    "Phospholipids.in.Chylomicrons.and.Extremely.Large.VLDL" = "XXL-VLDL-PL",
    "Cholesterol.in.Chylomicrons.and.Extremely.Large.VLDL" = "XXL-VLDL-C",
    "Cholesteryl.Esters.in.Chylomicrons.and.Extremely.Large.VLDL" = "XXL-VLDL-CE",
    "Free.Cholesterol.in.Chylomicrons.and.Extremely.Large.VLDL" = "XXL-VLDL-FC",
    "Triglycerides.in.Chylomicrons.and.Extremely.Large.VLDL" = "XXL-VLDL-TG",
    "Concentration.of.Very.Large.VLDL.Particles" = "XL-VLDL-P",
    "Total.Lipids.in.Very.Large.VLDL" = "XL-VLDL-L",
    "Phospholipids.in.Very.Large.VLDL" = "XL-VLDL-PL",
    "Cholesterol.in.Very.Large.VLDL" = "XL-VLDL-C",
    "Cholesteryl.Esters.in.Very.Large.VLDL" = "XL-VLDL-CE",
    "Free.Cholesterol.in.Very.Large.VLDL" = "XL-VLDL-FC",
    "Triglycerides.in.Very.Large.VLDL" = "XL-VLDL-TG",
    "Concentration.of.Large.VLDL.Particles" = "L-VLDL-P",
    "Total.Lipids.in.Large.VLDL" = "L-VLDL-L",
    "Phospholipids.in.Large.VLDL" = "L-VLDL-PL",
    "Cholesterol.in.Large.VLDL" = "L-VLDL-C",
    "Cholesteryl.Esters.in.Large.VLDL" = "L-VLDL-CE",
    "Free.Cholesterol.in.Large.VLDL" = "L-VLDL-FC",
    "Triglycerides.in.Large.VLDL" = "L-VLDL-TG",
    "Concentration.of.Medium.VLDL.Particles" = "M-VLDL-P",
    "Total.Lipids.in.Medium.VLDL" = "M-VLDL-L",
    "Phospholipids.in.Medium.VLDL" = "M-VLDL-PL",
    "Cholesterol.in.Medium.VLDL" = "M-VLDL-C",
    "Cholesteryl.Esters.in.Medium.VLDL" = "M-VLDL-CE",
    "Free.Cholesterol.in.Medium.VLDL" = "M-VLDL-FC",
    "Triglycerides.in.Medium.VLDL" = "M-VLDL-TG",
    "Concentration.of.Small.VLDL.Particles" = "S-VLDL-P",
    "Total.Lipids.in.Small.VLDL" = "S-VLDL-L",
    "Phospholipids.in.Small.VLDL" = "S-VLDL-PL",
    "Cholesterol.in.Small.VLDL" = "S-VLDL-C",
    "Cholesteryl.Esters.in.Small.VLDL" = "S-VLDL-CE",
    "Free.Cholesterol.in.Small.VLDL" = "S-VLDL-FC",
    "Triglycerides.in.Small.VLDL" = "S-VLDL-TG",
    "Concentration.of.Very.Small.VLDL.Particles" = "XS-VLDL-P",
    "Total.Lipids.in.Very.Small.VLDL" = "XS-VLDL-L",
    "Phospholipids.in.Very.Small.VLDL" = "XS-VLDL-PL",
    "Cholesterol.in.Very.Small.VLDL" = "XS-VLDL-C",
    "Cholesteryl.Esters.in.Very.Small.VLDL" = "XS-VLDL-CE",
    "Free.Cholesterol.in.Very.Small.VLDL" = "XS-VLDL-FC",
    "Triglycerides.in.Very.Small.VLDL" = "XS-VLDL-TG",
    "Concentration.of.IDL.Particles" = "IDL-P",
    "Total.Lipids.in.IDL" = "IDL-L",
    "Phospholipids.in.IDL" = "IDL-PL",
    "Cholesterol.in.IDL" = "IDL-C",
    "Cholesteryl.Esters.in.IDL" = "IDL-CE",
    "Free.Cholesterol.in.IDL" = "IDL-FC",
    "Triglycerides.in.IDL" = "IDL-TG",
    "Concentration.of.Large.LDL.Particles" = "L-LDL-P",
    "Total.Lipids.in.Large.LDL" = "L-LDL-L",
    "Phospholipids.in.Large.LDL" = "L-LDL-PL",
    "Cholesterol.in.Large.LDL" = "L-LDL-C",
    "Cholesteryl.Esters.in.Large.LDL" = "L-LDL-CE",
    "Free.Cholesterol.in.Large.LDL" = "L-LDL-FC",
    "Triglycerides.in.Large.LDL" = "L-LDL-TG",
    "Concentration.of.Medium.LDL.Particles" = "M-LDL-P",
    "Total.Lipids.in.Medium.LDL" = "M-LDL-L",
    "Phospholipids.in.Medium.LDL" = "M-LDL-PL",
    "Cholesterol.in.Medium.LDL" = "M-LDL-C",
    "Cholesteryl.Esters.in.Medium.LDL" = "M-LDL-CE",
    "Free.Cholesterol.in.Medium.LDL" = "M-LDL-FC",
    "Triglycerides.in.Medium.LDL" = "M-LDL-TG",
    "Concentration.of.Small.LDL.Particles" = "S-LDL-P",
    "Total.Lipids.in.Small.LDL" = "S-LDL-L",
    "Phospholipids.in.Small.LDL" = "S-LDL-PL",
    "Cholesterol.in.Small.LDL" = "S-LDL-C",
    "Cholesteryl.Esters.in.Small.LDL" = "S-LDL-CE",
    "Free.Cholesterol.in.Small.LDL" = "S-LDL-FC",
    "Triglycerides.in.Small.LDL" = "S-LDL-TG",
    "Concentration.of.Very.Large.HDL.Particles" = "XL-HDL-P",
    "Total.Lipids.in.Very.Large.HDL" = "XL-HDL-L",
    "Phospholipids.in.Very.Large.HDL" = "XL-HDL-PL",
    "Cholesterol.in.Very.Large.HDL" = "XL-HDL-C",
    "Cholesteryl.Esters.in.Very.Large.HDL" = "XL-HDL-CE",
    "Free.Cholesterol.in.Very.Large.HDL" = "XL-HDL-FC",
    "Triglycerides.in.Very.Large.HDL" = "XL-HDL-TG",
    "Concentration.of.Large.HDL.Particles" = "L-HDL-P",
    "Total.Lipids.in.Large.HDL" = "L-HDL-L",
    "Phospholipids.in.Large.HDL" = "L-HDL-PL",
    "Cholesterol.in.Large.HDL" = "L-HDL-C",
    "Cholesteryl.Esters.in.Large.HDL" = "L-HDL-CE",
    "Free.Cholesterol.in.Large.HDL" = "L-HDL-FC",
    "Triglycerides.in.Large.HDL" = "L-HDL-TG",
    "Concentration.of.Medium.HDL.Particles" = "M-HDL-P",
    "Total.Lipids.in.Medium.HDL" = "M-HDL-L",
    "Phospholipids.in.Medium.HDL" = "M-HDL-PL",
    "Cholesterol.in.Medium.HDL" = "M-HDL-C",
    "Cholesteryl.Esters.in.Medium.HDL" = "M-HDL-CE",
    "Free.Cholesterol.in.Medium.HDL" = "M-HDL-FC",
    "Triglycerides.in.Medium.HDL" = "M-HDL-TG",
    "Concentration.of.Small.HDL.Particles" = "S-HDL-P",
    "Total.Lipids.in.Small.HDL" = "S-HDL-L",
    "Phospholipids.in.Small.HDL" = "S-HDL-PL",
    "Cholesterol.in.Small.HDL" = "S-HDL-C",
    "Cholesteryl.Esters.in.Small.HDL" = "S-HDL-CE",
    "Free.Cholesterol.in.Small.HDL" = "S-HDL-FC",
    "Triglycerides.in.Small.HDL" = "S-HDL-TG"
  )
  
  colnames(BBmetabolomics) <- unname(sapply(colnames(BBmetabolomics), function(x) name_mapping[x]))
  # -----
  
  # Subset to the Complete Cases of Diet Samples
  # --------------------------------------------
  
  # Dealing with the diet data
  my_string <- "Diet"
  my_cols <- grep(my_string, colnames(BBdata_metadata), value = TRUE)
  BBdiet <- BBdata_metadata[, my_cols]
  
  # Formatting the names of the metabolites
  my_string <- ".Diet...0.0" 
  new_colnames <- gsub(my_string, "", colnames(BBdiet))
  colnames(BBdiet) <- new_colnames
  
  # filter out columns not to be used (have numerous NAs {type}, hard to quantify, not a focus point {salt})
  columns_to_remove <- grepl("type|Hot|Variation|Major|questionnaires|Never|Salt", colnames(BBdiet))
  BBdiet.feat <- BBdiet[, !columns_to_remove]
  
  # filter out incomplete samples across the remaining columns
  BBdiet.feat.samp <- BBdiet.feat[complete.cases(BBdiet.feat), ]
  rows_to_remove <- apply(BBdiet.feat.samp == -3 | BBdiet.feat.samp == -1, 1, any)
  BBdiet.feat.samp <- BBdiet.feat.samp[!rows_to_remove, ]
  
  # subset the metadata and metabolite to the common samples
  BBmetabolomics <- BBmetabolomics[rownames(BBmetabolomics) %in% rownames(BBdiet.feat.samp), ]
  BBmetadata <- BBmetadata[rownames(BBmetadata) %in% rownames(BBdiet.feat.samp), ]
  # -----

  # Removing non-baseline IBD individuals
  # -------------
  # Adding a column of difference between age at diagnosis and age at measurement
  BBmetadata <- BBmetadata %>%
    mutate(difference = ifelse(!is.na(Age_value) & !is.na(Age.at.measurement), Age.at.measurement - Age_value, NA))
  
  # Storing the nonIBD sample names
  nonIBD_sample_names <- rownames(BBmetadata[which(BBmetadata$Label == "0"), ])
  
  # Keeping IBD samples that had a history of IBD at baseline
  IBD_sample_names_to_keep <- rownames(BBmetadata[which(!is.na(BBmetadata$difference) & BBmetadata$difference >= 0), ])
  
  # Combining the two sets of sample names
  all_sample_names_to_keep <- union(nonIBD_sample_names, IBD_sample_names_to_keep)
  
  # Subsetting BBmetabolomics and BBmetadata and BBdiet based on the combined sample names
  BBmetabolomics <- BBmetabolomics[rownames(BBmetabolomics) %in% all_sample_names_to_keep, ]
  BBmetadata <- BBmetadata[rownames(BBmetadata) %in% all_sample_names_to_keep, ]
  BBdiet.feat.samp <- BBdiet.feat.samp[rownames(BBdiet.feat.samp) %in% all_sample_names_to_keep, ]
  # -----
  
  # Continue editting the Diet Data
  # -------------------------------
  
  BBdiet_numerical_features_vector <- c(
    'Cooked_vegetable_intake', 
    'Salad_._raw_vegetable_intake',  
    'Fresh_fruit_intake', 
    'Dried_fruit_intake', 
    'Bread_intake', 
    'Cereal_intake', 
    'Tea_intake', 
    'Coffee_intake', 
    'Water_intake' 
  )
  BBdiet_categorical_features_vector <- colnames(BBdiet.feat.samp)[!(colnames(BBdiet.feat.samp) %in% 
                                                                       BBdiet_numerical_features_vector)]
  BBdiet_numerical <- BBdiet.feat.samp[, BBdiet_numerical_features_vector]
  BBdiet_categorical <- BBdiet.feat.samp[, BBdiet_categorical_features_vector]
  # replacing -10 with a numerical average equivalence of "less than one"
  BBdiet_numerical <- replace(BBdiet_numerical, (BBdiet_numerical == "-10"), 0.5)
  
  # Robust Scaling function
  robust_scale <- function(x) {
    (x - median(x)) / IQR(x)
  }
  
  BBdiet_numerical_robust_scaled <- as.data.frame(apply(BBdiet_numerical, 2, robust_scale))
  
  # to replace the category forms with numerically equivalent values
  value_mapping <- list(
    `0` = 0, 
    `1` = 0.07, 
    `2` = 0.14, 
    `3` = 0.43, 
    `4` = 0.79, 
    `5` = 1
  )
  
  BBdiet_final <- cbind(BBdiet_numerical_robust_scaled, BBdiet_categorical)
  # -----
  
  # Renaming BBdiet names
  # ---------------------
  
  name_mapping <- c(
    "Cooked_vegetable_intake" = "Cooked Vegetables",
    "Salad_._raw_vegetable_intake" = "Raw Vegetables",
    "Fresh_fruit_intake" = "Fresh Fruits",
    "Dried_fruit_intake" = "Dried Fruits",
    "Bread_intake" = "Bread",
    "Cereal_intake" = "Cereal",
    "Tea_intake" = "Tea",
    "Coffee_intake" = "Coffee",
    "Water_intake" = "Water",
    "Beef_intake" = "Beef",
    "Cheese_intake" = "Cheese",
    "Lamb.mutton_intake" = "Lamb-Mutton",
    "Non_oily_fish_intake" = "Non-Oily Fish",
    "Oily_fish_intake" = "Oily Fish",
    "Pork_intake" = "Pork",
    "Poultry_intake" = "Poultry",
    "Processed_meat_intake" = "Processed Meat"
  )
  
  colnames(BBdiet_final) <- unname(sapply(colnames(BBdiet_final), function(x) name_mapping[x]))
  # -----
  
  return(list(BBmetabolomics = BBmetabolomics, BBmetadata = BBmetadata, BBdiet = BBdiet_final))
  
}

# Function to impute, standardize, and transform Cohort 1 (UK Biobank) data
impute_standardize_transform_cohort1 <- function(BBmetabolomics) {
  # Calculate the percentage of missing and 'zero'values for each column in BBmetabolomics
  BBmetabolomics_missing_perc <- mean(is.na(BBmetabolomics)) * 100
  BBmetabolomics_zero_perc <- mean(BBmetabolomics == 0, na.rm = TRUE) * 100 
  
  # Impute missing values using missForest
  BBmetabolomics.imp <- missForest(BBmetabolomics, maxiter = 10, ntree = 100)
  BBmetabolomics.imp <- BBmetabolomics.imp$ximp

  # Perform Pareto scaling (column-wise)
  BBmetabolomics.imp.sc <- scale(BBmetabolomics.imp, center = FALSE, scale = sqrt(apply(BBmetabolomics.imp, 2, var)))
  
  # Apply log transformation after shifting by 1
  BBmetabolomics.imp.sc.trans <- as.data.frame(log(BBmetabolomics.imp.sc + abs(min(BBmetabolomics.imp.sc) + 1)))
  
  # Feature Selection
  # Remove columns that are a sum of other columns to avoid collinearity and 
  # to help the correlation remover make better decisions
  
  # Find column indices containing the string "Total" (don't consider C, if C = A + B + etc.)
  # and don't consider particle size
  columns_to_remove <- grep("Total|size", colnames(BBmetabolomics.imp.sc.trans))
  
  # Remove columns with "Total" in their names
  BBmetabolomics.imp.sc.trans.filtered <- BBmetabolomics.imp.sc.trans[, -columns_to_remove]
  
  # Remove highly correlated features
  # use spearman because of the possibly non linear nature of the data
  cor_mat <- cor(BBmetabolomics.imp.sc.trans.filtered, method = "spearman")
  
  highly_corr_feat <- findCorrelation(cor_mat, cutoff=0.9)
  highly_corr_feat = sort(highly_corr_feat)
  
  removed_highly_cor_feat_data <- BBmetabolomics.imp.sc.trans.filtered[, -c(highly_corr_feat)]
  return(list(BBprocessed = BBmetabolomics.imp.sc.trans,
              BBprocessed_filtered = removed_highly_cor_feat_data,
              BBhigh_corr_features = highly_corr_feat,
              BBcorr_mat = cor_mat,
              BBmissingno = BBmetabolomics_missing_perc,
              BBzerono = BBmetabolomics_zero_perc))
}

# Function to read Cohort 2 (HMP2 data)
read_cohort2_data <- function(hmp2_metabolomics_file, hmp2_metadata_file) {
  HMPmetabolomics <- read.csv(hmp2_metabolomics_file, header = TRUE, row.names = NULL)
  HMPmetadata <- read.csv(hmp2_metadata_file, header = TRUE, row.names = NULL)
  
  return(list(HMPmetabolomics = HMPmetabolomics, 
              HMPmetadata = HMPmetadata))
}

# Function to preprocess Cohort 2 (HMP2 data)
preprocess_cohort2_data <- function(HMPmetabolomics, HMPmetadata) {
  
  # Preprocess the Metadata
  # ------------------------
  
  # Change name of "diagnosis" to "Label" for MetaboAnalyst format
  names(HMPmetadata)[names(HMPmetadata) == "diagnosis"] <- "Label"
  
  # Get the complete cases of the HMPdiet
  HMPdiet <- HMPmetadata[, c(which(colnames(HMPmetadata) == "External.ID"), 72:81, 83:92, 101)]
  HMPdiet <- HMPdiet[, -which(colnames(HMPdiet) == "Probiotic")]
  HMPdiet[HMPdiet == ""] <- NA
  HMPdiet <- HMPdiet[complete.cases(HMPdiet), ]
  HMPmetadata <- HMPmetadata[(HMPmetadata$External.ID %in% HMPdiet$External.ID), ]
  # -----
  
  # Remove non-baseline IBD individuals
  # -----------------------------------
  
  HMPmetadata <- HMPmetadata %>% 
    mutate(difference = ifelse(!is.na(consent_age) & !is.na(Age.at.diagnosis), consent_age - Age.at.diagnosis, NA))
  
  # Storing the nonIBD sample names
  nonIBD_sample_names <- HMPmetadata[which(HMPmetadata$Label == "nonIBD"), "External.ID"]
  
  # Keeping IBD samples that had a history of IBD at baseline
  IBD_sample_names_to_keep <- HMPmetadata[which(!is.na(HMPmetadata$difference) & HMPmetadata$difference >= 0), "External.ID"]
  
  # Combining the two sets of sample names
  all_sample_names_to_keep <- union(nonIBD_sample_names, IBD_sample_names_to_keep)
  
  # Subset HMPmetadata based on the combined sample names
  HMPmetadata <- HMPmetadata[HMPmetadata$External.ID %in% all_sample_names_to_keep, ]
  
  # Sub HMP metadata to the ones with the metabolites
  HMPmetadata <- HMPmetadata %>% filter(data_type %in% "metabolomics")
  rownames(HMPmetadata) <- HMPmetadata$External.ID
  HMPmetadata$External.ID <- NULL
 
  # Edit the metadata to fill in the NA sections
  HMPmetadata <- replace(HMPmetadata, is.na(HMPmetadata), "NA")
  HMPmetadata <- replace(HMPmetadata, HMPmetadata == "", "NA")
  # ------
  
  # Pre-process the Metabolomics Data
  # ---------------------------------
  
  # Select one metabolomic method (HILIC-pos) in HMPmetabolomics
  HMPmetabolomics <- subset(HMPmetabolomics, Method == "HILIC-pos")
  
  # Subset HMPmetabolomics to only the metabolites that exist
  HMPmetabolomics <- subset(HMPmetabolomics, Metabolite != "")
  
  # Keep only the samples and the metabolite columns in HMPmetabolomics
  HMPmetabolomics <- HMPmetabolomics[, -c(1:5,7), drop = FALSE]
  
  # Remove duplicated rows (metabolites), keeping the one with fewer NA values
  na_counts <- rowSums(is.na(HMPmetabolomics))
  duplicated_rows <- duplicated(HMPmetabolomics$Metabolite) | duplicated(HMPmetabolomics$Metabolite, fromLast = TRUE)
  keep_rows <- !duplicated_rows | (duplicated_rows & na_counts == min(na_counts[duplicated_rows]))
  HMPmetabolomics <- HMPmetabolomics[keep_rows, ]
  
  # Set row names and drop the metabolite column
  row.names(HMPmetabolomics) <- HMPmetabolomics$Metabolite
  HMPmetabolomics <- HMPmetabolomics %>% select(-Metabolite)
  
  # Transpose HMPmetabolomics to fit the form of BBmetabolomics (samples as rows, metabolites as columns)
  HMPmetabolomics <- as.data.frame(t(HMPmetabolomics))
  
  # --> INSERT FILTERING STEP HERE (E.G. SAMPLES THAT HAVE >50% MISSING)
  # This metabolomic dataset has about <5% missing values, so no need for filtering
  
  # Sub HMPmetabolomics to available samples in HMPmetadata
  HMPmetabolomics <- HMPmetabolomics[rownames(HMPmetabolomics) %in% rownames(HMPmetadata), ]
  
  # In case there was some filtering step in metabolomics that took out more samples
  HMPmetadata <- HMPmetadata[rownames(HMPmetadata) %in% rownames(HMPmetabolomics), ]
  # ------
  
  # Renaming the BBmetabolite names
  # -------------------------------
  
  # Mapping of old names to new names
  name_mapping <- c(
    "trimethylbenzene" = "Trimethylbenzene",
    "bilirubin" = "Bilirubin",
    "cytidine" = "Cytidine",
    "hippurate" = "Hippurate",
    "hypoxanthine" = "Hypoxan",
    "inosine" = "Inosine",
    "pantothenate" = "Pantothene",
    "glycine" = "Glycine",
    "alanine" = "Alanine",
    "serine" = "Serine",
    "methionine" = "Methionine",
    "aspartate" = "Aspartate",
    "glutamate" = "Glutamate",
    "glutamine" = "Glutamine",
    "histidine" = "Histidine",
    "arginine" = "Arginine",
    "lysine" = "Lysine",
    "valine" = "Valine",
    "leucine" = "Leucine",
    "isoleucine" = "Isoleucine",
    "phenylalanine" = "Phenylalanine",
    "tyrosine" = "Tyrosine",
    "tryptophan" = "Tryptophan",
    "proline" = "Proline",
    "hydroxyproline" = "Hydroxyproline",
    "ornithine" = "Ornithine",
    "citrulline" = "Citrulline",
    "taurine" = "Taurine",
    "5-hydroxytryptophan" = "5-HTP",
    "serotonin" = "Serotonin",
    "dimethylglycine" = "DMG",
    "ADMA/SDMA" = "ADMA/SDMA",
    "NMMA" = "NMMA",
    "aminoisobutyric acid/GABA" = "AIB/GABA",
    "kynurenic acid" = "KA",
    "1-methylhistamine" = "1-MH",
    "histamine" = "Histamine",
    "N-carbamoyl-beta-alanine" = "N-CBA",
    "niacinamide" = "NAM",
    "betaine" = "Betaine",
    "choline" = "Choline",
    "alpha-glycerophosphocholine" = "alpha-GPC",
    "acetylcholine" = "Acetylcholine",
    "spermidine" = "Spermidine",
    "creatine" = "Creatine",
    "creatinine" = "Creatinine",
    "trimethylamine-N-oxide" = "TMAO",
    "adenosine" = "Adenosine",
    "cytosine" = "Cytosine",
    "2-deoxyadenosine" = "2-dAdo",
    "2-deoxycytidine" = "2-dC",
    "cotinine" = "Cotinine",
    "pipecolic acid" = "Pipecolic acid",
    "5-aminolevulinic acid" = "5-ALA",
    "1-methylnicotinamide" = "1-MNA",
    "butyrobetaine" = "Butyrobetaine",
    "putrescine" = "Putrescine",
    "methionine sulfoxide" = "MetO",
    "beta-alanine" = "Beta-ala",
    "anserine" = "Anserine",
    "carnitine" = "Carn",
    "C2 carnitine" = "C2 Carnitine",
    "C3 carnitine" = "C3 Carnitine",
    "C3-DC-CH3 carnitine" = "C3-DC-CH3 Carn",
    "C4 carnitine" = "C4 Carnitine",
    "C5 carnitine" = "C5 Carnitine",
    "C8 carnitine" = "C8 Carnitine",
    "C9 carnitine" = "C9 Carnitine",
    "C10 carnitine" = "C10 Carnitine",
    "C10:2 carnitine" = "C10:2 Carnitine",
    "C12 carnitine" = "C12 Carnitine",
    "C12:1 carnitine" = "C12:1 Carnitine",
    "C14 carnitine" = "C14 Carnitine",
    "C14:1 carnitine" = "C14:1 Carnitine",
    "C14:2 carnitine" = "C14:2 Carnitine",
    "C16 carnitine" = "C16 Carnitine",
    "C16-OH carnitine" = "C16-OH Carnitine",
    "C18 carnitine" = "C18 Carnitine",
    "C18:1 carnitine" = "C18:1 Carnitine",
    "C18:1-OH carnitine" = "C18:1-OH Carn",
    "C18:2 carnitine" = "C18:2 Carnitine",
    "C20 carnitine" = "C20 Carnitine",
    "C20:4 carnitine" = "C20:4 Carnitine",
    "1-methylguanine" = "1-MeG",
    "1-methylguanosine" = "1-MeGuo",
    "21-deoxycortisol" = "21-Deoxycortisol",
    "2-aminoheptanoic acid" = "2-AHA",
    "2'-O-methyladenosine" = "2'-O-MeAdo",
    "3-methylhistidine" = "3-MeHis",
    "3-methylxanthine" = "3-MeX",
    "3'-O-methyladenosine" = "3'-O-MeAdo",
    "4-aminophenol" = "4-Aminophenol",
    "4-guanidinobutanoic acid" = "4-GuanidinoBA",
    "4-hydroxy-3-methylacetophenone" = "4H-3MeAcetophenone",
    "5-acetylamino-6-amino-3-methyluracil" = "5-Ac-6-Am-3-MeU",
    "7-methylguanine" = "7-MeG",
    "7-methylxanthine" = "7-MeX",
    "8-hydroxy-deoxyguanosine" = "8-OH-dGuo",
    "acetaminophen" = "Acetaminophen",
    "acetyl-galactosamine" = "Acetyl-GalN",
    "agmatine" = "Agmatine",
    "alloisoleucine" = "Alloisoleucine",
    "alpha-hydroxymetoprolol" = "alpha-OH-Metoprolol",
    "atenolol" = "Atenolol",
    "beta-guanidinopropionic acid" = "beta-GPA",
    "biliverdin" = "Biliverdin",
    "biotin" = "Biotin",
    "cadaverine" = "Cadaverine",
    "caffeine" = "Caffeine",
    "cortisol" = "Cortisol",
    "alanylalanine" = "AlaAla",
    "diacetylspermine" = "DAS",
    "ectoine" = "Ectoine",
    "gabapentin" = "Gabapentin",
    "guanidoacetic acid" = "GAA",
    "histidinol" = "Histidinol",
    "homoarginine" = "HomoArg",
    "homocitrulline" = "Homocit",
    "hydroxycotinine" = "Hydroxycotinine",
    "imidazole propionate" = "ImP",
    "imidazoleacetic acid" = "IAA",
    "linoleoyl ethanolamide" = "LEA",
    "metformin" = "Metformin",
    "methylimidazole acetic acid" = "MIAPA",
    "metronidazole" = "Metronidazole",
    "N1,N12-diacetylspermine" = "N1,N12-DAS",
    "N1-acetylspermidine" = "N1-ASPD",
    "N1-acetylspermine" = "N1-ASP",
    "N1-methyl-2-pyridone-5-carboxamide" = "NMPC",
    "N2,N2-dimethylguanosine" = "N2,N2-DiMeGuo",
    "N6-acetyllysine" = "N6-AcLys",
    "N-acetyalanine" = "N-AcAla",
    "N-acetylaspartic acid" = "N-AcAsp",
    "N-acetylglutamine" = "N-AcGln",
    "N-acetylhistamine" = "N-AcHSM",
    "N-acetylhistidine" = "N-AcHis",
    "N-acetylornithine" = "N-AcOrn",
    "N-acetylputrescine" = "N-AcPut",
    "N-alpha-acetylarginine" = "N-alpha-AcArg",
    "N6,N6-dimethyllysine" = "N6 N6-DiMeLys",
    "N6,N6,N6-trimethyllysine" = "N6 N6 N6-TriMeLys",
    "nicotinuric acid" = "NUA",
    "N-methylproline" = "N-MePro",
    "phytosphingosine" = "Phytosphing",
    "piperine" = "Piperine",
    "proline betaine" = "ProBet",
    "pterin" = "Pterin",
    "pyridoxine" = "Pyridoxine",
    "quinine" = "Quinine",
    "S-adenosylmethionine" = "SAM",
    "S-methylcysteine-S-oxide" = "S-MeCys-S-Oxide",
    "sphinganine" = "Sphinganine",
    "sulfapyridine" = "Sulfapyridine",
    "threosphingosine" = "Threo-Sph",
    "trigonelline" = "Trigonelline",
    "urobilin" = "Urobilin",
    "C14:0 LPC" = "C14:0 LPC",
    "C16:1 LPC" = "C16:1 LPC",
    "C16:0 LPC" = "C16:0 LPC",
    "C18:3 LPC" = "C18:3 LPC",
    "C18:2 LPC" = "C18:2 LPC",
    "C18:1 LPC" = "C18:1 LPC",
    "C18:0 LPC" = "C18:0 LPC",
    "C20:1 LPC" = "C20:1 LPC",
    "C16:1 LPC plasmalogen" = "C16:1 LPC-P",
    "C18:1 LPC plasmalogen" = "C18:1 LPC-P",
    "C16:0 LPE" = "C16:0 LPE",
    "C18:0 LPE-A" = "C18:0 LPE-A",
    "C18:0 LPE-B" = "C18:0 LPE-B",
    "C22:6 LPE" = "C22:6 LPE",
    "C36:2 PC" = "C36:2 PC",
    "sphingosine-isomer1" = "Sphingo-iso1",
    "sphingosine-isomer2" = "Sphingo-iso2",
    "sphingosine-isomer3" = "Sphingo-iso3",
    "C14:0 SM" = "C14:0 SM",
    "C16:0 SM" = "C16:0 SM"
  )
  
  # Replace column names
  colnames(HMPmetabolomics) <- unname(sapply(colnames(HMPmetabolomics), function(x) name_mapping[x]))
  # ------
  
  # Continue editing the Diet Data
  # ------------------------------
  
  # remove duplicated rows
  HMPdiet <- HMPdiet[!duplicated(HMPdiet), ]
  
  # Sub the diet data to the final samples in the metadata
  HMPdiet <- HMPdiet[HMPdiet$External.ID %in% rownames(HMPmetadata), ]
  rownames(HMPdiet) <- HMPdiet$External.ID
  HMPdiet$External.ID <- NULL
  
  value_mapping <- list(
    `No, I did not consume these products in the last 7 days` = 0,
    `Within the past 4 to 7 days` = 0.2,
    `Within the past 2 to 3 days` = 0.58,
    `Yesterday, 1 to 2 times` = 0.9,
    `Yesterday, 3 or more times` = 1
  )
  
  HMPdiet_replaced <- HMPdiet %>%
    mutate_at(vars(everything()), ~ as.numeric(value_mapping[as.character(.)]))
  # ------
  
  # Edit the Diet names
  # -------------------
  
  name_mapping <- c(
    "Soft.drinks..tea.or.coffee.with.sugar..corn.syrup..maple.syrup..cane.sugar..etc." = "Soft Drinks",
    "Diet.soft.drinks..tea.or.coffee.with.sugar..Stevia..Equal..Splenda.etc." = "Diet Soft Drinks",
    "Fruit.juice..orange..apple..cranberry..prune.etc.." = "Fruit Juice",
    "Water" = "Water",
    "Alcohol..beer..brandy..spirits..hard.liquor..wine..aperitif..etc.." = "Alcohol",
    "Yogurt.or.other.foods.containing.active.bacterial.cultures..kefir..sauerkraut." = "Foods with Active Bacteria",
    "Dairy..milk..cream..ice.cream..cheese..cream.cheese." = "Dairy",
    "Fruits..no.juice...Apples..raisins..bananas..oranges..strawberries..blueberries" = "Fruits",
    "Vegetables..salad..tomatoes..onions..greens..carrots..peppers..green.beans..etc." = "Vegetables",
    "Beans..tofu..soy..soy.burgers..lentils..Mexican.beans..lima.beans.etc." = "Beans",
    "Whole.grains..wheat..oats..brown.rice..rye..quinoa..wheat.bread..wheat.pasta." = "Whole Grains",
    "Starch..white.rice..bread..pizza..potatoes..yams..cereals..pancakes..etc.." = "Starch",
    "Eggs" = "Eggs",
    "Processed.meat..other.red.or.white.meat.such.as.lunch.meat..ham..salami..bologna" = "Processed Meat",
    "Red.meat..beef..hamburger..pork..lamb." = "Red Meat",
    "White.meat..chicken..turkey..etc.." = "White Meat",
    "Shellfish..shrimp..lobster..scallops..etc.." = "Shellfish",
    "Fish..fish.nuggets..breaded.fish..fish.cakes..salmon..tuna..etc.." = "Fish",
    "Sweets..pies..jam..chocolate..cake..cookies..etc.." = "Sweets",
    "Tea.or.coffee.no.sugar.and.no.sugar.replacement" = "Beverage without Sugar"
  )
  
  colnames(HMPdiet_replaced) <- unname(sapply(colnames(HMPdiet_replaced), 
                                              function(x) name_mapping[x]))
  
  # Not interested in Alcohol
  HMPdiet_replaced <- HMPdiet_replaced[, !(colnames(HMPdiet_replaced) == "Alcohol")]
  # ------
  
  # # Trim down the metadata
  # # -----------------------
  # 
  # # Define the threshold for removing columns with one unique value
  # unique_threshold <- 1  # All values are the same
  # # Define the threshold for removing columns with more than 50% "NA" values
  # NA_threshold <- 0.5  # More than 50% "NA"
  # # Combine all removal operations in a single pipeline
  # HMPmetadata_selected_names <- HMPmetadata %>%
  #   select(-where(~ n_distinct(.) <= unique_threshold)) %>%
  #   select(-where(~ mean(. == "NA", na.rm = TRUE) > NA_threshold))
  # 
  # # DIET IDs to remove
  # dietIDs_to_remove <- c(17:26, 28:37, 46)
  # 
  # # Unnecessary columns to remove
  # unnecessary_columns_to_remove <- c(1:10, 55, 57:63)
  # 
  # HMPmetadata_selected_names <- HMPmetadata_selected_names %>%
  #   select(-all_of(union(dietIDs_to_remove, unnecessary_columns_to_remove)))
  # 
  # HMPmetadata_names <- names(HMPmetadata_selected_names)
  # ------
  
  return(list(HMPdiet = HMPdiet_replaced,
              HMPmetabolomics = HMPmetabolomics, 
              HMPmetadata = HMPmetadata))
}

# Function to impute, standardize, and transform Cohort 2 (HMP2 data) data
impute_standardize_transform_cohort2 <- function(HMPmetabolomics) {
  # Calculate missing and zero percentages
  HMPmetabolomics_missing_perc <- mean(is.na(HMPmetabolomics)) * 100
  HMPmetabolomics_zero_perc <- mean(HMPmetabolomics == "0", na.rm = TRUE) * 100
  
  # Change zero values with half of the minimum positive value
  replacezero <- function(x) "[<-"(x, !x | is.na(x), min(x[x > 0], na.rm = TRUE) / 2)
  HMPmetabolomics.imp <- as.data.frame(apply(HMPmetabolomics, 2, replacezero))
  
  # Perform scaling and transform
  HMPmetabolomics.imp.sc <- scale(HMPmetabolomics.imp, center = FALSE, scale = sqrt(apply(HMPmetabolomics.imp, 2, var)))
  HMPmetabolomics.imp.sc.trans <- as.data.frame(log(HMPmetabolomics.imp.sc + abs(min(HMPmetabolomics.imp.sc)) + 1))
  
  # Remove highly collinear features
  cor_mat <- cor(HMPmetabolomics.imp.sc.trans, method = "spearman")
  
  highly_corr_feat <- findCorrelation(cor_mat, cutoff=0.9)
  highly_corr_feat = sort(highly_corr_feat)
  
  removed_highly_cor_feat_data <- HMPmetabolomics.imp.sc.trans[, -c(highly_corr_feat)]
  return(list(HMPprocessed = HMPmetabolomics.imp.sc.trans,
              HMPprocessed_filtered = removed_highly_cor_feat_data,
              HMPhigh_corr_features = highly_corr_feat,
              HMPcorr_mat = cor_mat,
              HMPmissingno = HMPmetabolomics_missing_perc,
              HMPzerono = HMPmetabolomics_zero_perc))
}

# Function to find which features were highly correlated to the ones that stayed
which_correlated <- function(metabolomics_data, high_corr_feat_vector, corr_matrix) {
  # Create an empty list to store results
  result_list <- list()
  
  # Loop through each metabolite in BBmetabolites
  for (i in colnames(metabolomics_data)) {
    # Initialize an empty character vector to store highly correlated features
    high_corr_features <- character(0)
    
    # Loop through each high correlation feature in BBhighcorrfeatures
    for (j in high_corr_feat_vector) {
      
      # Skip the current metabolite
      if (i == j) {
        next
      }
      
      # Check if correlation is >= 0.9
      if (corr_matrix[i, j] >= 0.9) {
        # Append highly correlated feature to the vector
        high_corr_features <- c(high_corr_features, j)
      }
    }
    
    # If there are highly correlated features, create a string with comma-separated values
    if (length(high_corr_features) > 0) {
      result_list[[i]] <- paste(high_corr_features, collapse = ", ")
    }
  }
  
  # Create a data frame from the result list
  result_df <- data.frame(
    metabolite = names(result_list),
    highly_correlated_features = unlist(result_list)
  )
  
  return(result_df)
  
}

# Function to format labels for both Cohort 1 and Cohort 2 data
format_labels <- function(cohort_metadata) {
  # Define a function for find and replace
  find_replace <- function(x) {
    x <- gsub("CD", "IBD", x)
    x <- gsub("UC", "IBD", x)
    x <- gsub("1", "IBD", x)
    x <- gsub("0", "nonIBD", x)
    return(x)
  }
  
  # Change the IBD and nonIBD labels to a unified form
  cohort_metadata$Label <- find_replace(cohort_metadata$Label)
  
  # Retrieve the label column from the Metadata
  label_data <- cohort_metadata %>% select(Label)
  
  return(list(label_data, cohort_metadata))
}

# Function to merge the data with the labels
merge_data_with_label <- function(BBmetabolomics_label, 
                                  BBmetabolomics_imp_sc_trans_filtered, BBmetabolomics_imp_sc_trans,
                                  HMPmetabolomics_label, 
                                  HMPmetabolomics_imp_sc_trans_filtered, HMPmetabolomics_imp_sc_trans) {
  
  # For Machine Learning and Exploration
  BBmetabolomics.imp.sc.trans.filt_with_label <- merge(BBmetabolomics_label, BBmetabolomics_imp_sc_trans_filtered, by = "row.names")
  HMPmetabolomics.imp.sc.trans.filt_with_label <- merge(HMPmetabolomics_label, HMPmetabolomics_imp_sc_trans_filtered, by = "row.names")
 
  # For keep-sake
  BBmetabolomics.imp.sc.trans_with_label <- merge(BBmetabolomics_label, BBmetabolomics_imp_sc_trans, by = "row.names")
  HMPmetabolomics.imp.sc.trans_with_label <- merge(HMPmetabolomics_label, HMPmetabolomics_imp_sc_trans, by = "row.names")
  
  # Return the merged datasets with labels
  return(list(BBmetabolomics_imp_sc_trans_with_label = BBmetabolomics.imp.sc.trans_with_label,
              HMPmetabolomics_imp_sc_trans_with_label = HMPmetabolomics.imp.sc.trans_with_label,
              BBmetabolomics_imp_sc_trans_filt_with_label = BBmetabolomics.imp.sc.trans.filt_with_label,
              HMPmetabolomics_imp_sc_trans_filt_with_label = HMPmetabolomics.imp.sc.trans.filt_with_label))
}

# Function to export processed data
format_and_export_data <- function(data_list, 
                                   BBmetadata, HMPmetadata, BBdiet, HMPdiet, 
                                   BBcor, HMPcor, BBCorrFeats, HMPCorrFeats, 
                                   output_directory) {
  
  BB_with_label <- data_list$BBmetabolomics_imp_sc_trans_with_label
  HMP_with_label <- data_list$HMPmetabolomics_imp_sc_trans_with_label
  BB.filt_with_label <- data_list$BBmetabolomics_imp_sc_trans_filt_with_label
  HMP.filt_with_label <- data_list$HMPmetabolomics_imp_sc_trans_filt_with_label
 
  # Set the rownames as the Sample IDs
  rownames(BB_with_label) <- BB_with_label$Row.names
  BB_with_label$Row.names <- NULL
  rownames(HMP_with_label) <- HMP_with_label$Row.names
  HMP_with_label$Row.names <- NULL
  rownames(BB.filt_with_label) <- BB.filt_with_label$Row.names
  BB.filt_with_label$Row.names <- NULL
  rownames(HMP.filt_with_label) <- HMP.filt_with_label$Row.names
  HMP.filt_with_label$Row.names <- NULL

  # Export the normal data
  write.csv(BB_with_label, file = paste0(output_directory, "BB_imp_sc_trans_label.csv"), row.names = TRUE)
  write.csv(HMP_with_label, file = paste0(output_directory, "HMP_imp_sc_trans_label.csv"), row.names = TRUE)
  write.csv(BB.filt_with_label, file = paste0(output_directory, "BB_imp_sc_trans_filt_label.csv"), row.names = TRUE)
  write.csv(HMP.filt_with_label, file = paste0(output_directory, "HMP_imp_sc_trans_filt_label.csv"), row.names = TRUE)
  write.csv(BBcor, file = paste0(output_directory, "BB_cor_mat.csv"), row.names = TRUE)
  write.csv(HMPcor, file = paste0(output_directory, "HMP_cor_mat.csv"), row.names = TRUE)
  write.csv(BBCorrFeats, file = paste0(output_directory, "BB_high_cor_features.csv"), row.names = TRUE)
  write.csv(HMPCorrFeats, file = paste0(output_directory, "HMP_high_cor_features.csv"), row.names = TRUE)
 
  # Export the metadata
  write.csv(BBmetadata[, c(1,2,4,5,10)], file = paste0(output_directory, "BBmetadata.csv"), row.names = TRUE) # only sections of interest
  write.csv(HMPmetadata, file = paste0(output_directory, "HMPmetadata.csv"), row.names = TRUE)
  write.csv(HMPmetadata[, c("Participant.ID", "week_num", "Age.at.diagnosis", "consent_age", "Label")], 
            file = paste0(output_directory, "HMPmetadata_timepoint.csv"), row.names = TRUE)

  # Export the Diet Data
  write.csv(BBdiet, file = paste0(output_directory, "BBdiet.csv"), row.names = TRUE)
  write.csv(HMPdiet, file = paste0(output_directory, "HMPdiet.csv"), row.names = TRUE)
}


#### ------------ BEGIN RUNNING -----------------
#### --------------------------------------------

# Uncomment to LOAD ALL PROCESSED DATA
# load (file = "Data Preprocessing.RData")

load_libraries()
# >> input working directory path <<
set_working_directory("C:/Users/lenovo/Documents/GitHub/XAIMetabolomeDiet/R Processed Data")

# Step 1: Read data for Cohort 1 (UK Biobank)
cohort1_data <- read_cohort1_data("../Cohort Data/UKBiobank Data/FeaturesData.csv", 
                                  "../Cohort Data/UKBiobank Data/DiagnosisData.csv")
BBdata_metadata <- cohort1_data$BBdata_metadata
BBdiagnosis_metadata <- cohort1_data$BBdiagnosis_metadata

# Step 2: Preprocess data for Cohort 1 (UK Biobank)
cohort1_processed_data <- preprocess_cohort1_data(BBdata_metadata, BBdiagnosis_metadata)
BBmetabolomics <- cohort1_processed_data$BBmetabolomics
BBdiet <- cohort1_processed_data$BBdiet
BBmetadata <- cohort1_processed_data$BBmetadata

# Step 3: Impute, standardize, and transform data for Cohort 1 (UK Biobank)
cohort1_transformed_data <- impute_standardize_transform_cohort1(BBmetabolomics)
BBmetabolomics_missing_perc <- cohort1_transformed_data$BBmissingno
BBmetabolomics_zero_perc <- cohort1_transformed_data$BBzerono
BBmetabolomics_imp_sc_trans <- cohort1_transformed_data$BBprocessed
BBmetabolomics_imp_sc_trans_filtered <- cohort1_transformed_data$BBprocessed_filtered
BBcor <- cohort1_transformed_data$BBcorr_mat
BBhighlycorrfeatures <- colnames(BBcor)[cohort1_transformed_data$BBhigh_corr_features]

# Step 4: Read in data for Cohort 2 (HMP2 data)
cohort2_data <- read_cohort2_data("../Cohort Data/HMP2DB Data/iHMP_metabolomics.csv", 
                                  "../Cohort Data/HMP2DB Data/hmp2_metadata.csv")
HMPmetabolomics_pre <- cohort2_data$HMPmetabolomics
HMPmetadata_pre <- cohort2_data$HMPmetadata

# Step 5: Preprocess data for Cohort 2 (HMP2 data)
cohort2_processed_data <- preprocess_cohort2_data(HMPmetabolomics_pre, HMPmetadata_pre)
HMPmetabolomics <- cohort2_processed_data$HMPmetabolomics
HMPdiet <- cohort2_processed_data$HMPdiet
HMPmetadata <- cohort2_processed_data$HMPmetadata 

# Step 6: Impute, standardize, and transform metabolomics data for Cohort 2 (HMP2 data)
cohort2_transformed_data <- impute_standardize_transform_cohort2(HMPmetabolomics)
HMPmetabolomics_missing_perc <- cohort2_transformed_data$HMPmissingno
HMPmetabolomics_zero_perc <- cohort2_transformed_data$HMPzerono
HMPmetabolomics_imp_sc_trans <- cohort2_transformed_data$HMPprocessed
HMPmetabolomics_imp_sc_trans_filtered <- cohort2_transformed_data$HMPprocessed_filtered
HMPcor <- cohort2_transformed_data$HMPcorr_mat
HMPhighlycorrfeatures <- colnames(HMPcor)[cohort2_transformed_data$HMPhigh_corr_features]

# Step 7: Find out which is highly correlated and was filtered out
BBCorrFeats <- which_correlated(BBmetabolomics_imp_sc_trans_filtered, BBhighlycorrfeatures, BBcor)
HMPCorrFeats <- which_correlated(HMPmetabolomics_imp_sc_trans_filtered, HMPhighlycorrfeatures, HMPcor)

# Step 8: Format labels for both Cohort 1 and Cohort 2 data
BBmetabolomics_label <- format_labels(BBmetadata)[[1]]
BBmetadata <- format_labels(BBmetadata)[[2]]
HMPmetabolomics_label <- format_labels(HMPmetadata)[[1]]
HMPmetadata <- format_labels(HMPmetadata)[[2]]

# Step 9: Merge processed data with labels
data_with_labels <- merge_data_with_label(BBmetabolomics_label, 
                                          BBmetabolomics_imp_sc_trans_filtered, BBmetabolomics_imp_sc_trans,
                                          HMPmetabolomics_label, 
                                          HMPmetabolomics_imp_sc_trans_filtered, HMPmetabolomics_imp_sc_trans)

# Step 10: Format and export processed data
format_and_export_data(data_with_labels, 
                       BBmetadata, HMPmetadata, BBdiet, HMPdiet, 
                       BBcor, HMPcor, BBCorrFeats, HMPCorrFeats, "")

save.image(file = "Data Preprocessing Jan22 neat backup.RData")
