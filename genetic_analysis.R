# Set working directory
setwd("~/Documents/Estatistica_mestrado")

# Load packages
require(openxlsx)
require(dplyr)
require(HardyWeinberg)
require(stringr)
require(SNPassoc)
require(haplo.stats)
require(textshape)
require(genetics)
require(snpStats)
require(LDheatmap)
require(viridis)
require(grid)

# Load files
zkv_genotypes <- read.xlsx("Banco ZIKV - IL1B.xlsx")
zkv_metadata <- read.xlsx("metadados_zikv.xlsx")

zkv_data <- merge(zkv_metadata, zkv_genotypes, by = "Sample")


#___ Hardy-Weinberg equilibrium ___#

#' Hardy-Weinberg Equilibrium Test
#'
#' Performs an exact Hardy-Weinberg equilibrium (HWE) test for a given SNP in a specific group (e.g., cases or controls).
#'
#' @param dataset A data frame containing genotypes and phenotype.
#' @param group Group value to subset the dataset (e.g., 1 for cases, 0 for controls).
#' @param snp The name of the SNP column (e.g., "rs4848306").
#'
#' @return An object of class \code{HWExact}, containing the HWE test results.
#' @export
run_HWE <- function(dataset, group, snp) {
  genotypes <- table(dataset[dataset$Group == group, snp])
  HWEex <- HWExact(genotypes, alternative = "two.sided", pvaluetype = "selome", eps=1e-10, x.linked = FALSE, verbose = TRUE)
}

HWEex_rs306_case <- run_HWE(zkv_data, 1, "rs4848306")
HWEex_rs306_ctrl <- run_HWE(zkv_data, 0, "rs4848306")
HWEex_rs623_case <- run_HWE(zkv_data, 1, "rs1143623")
HWEex_rs623_ctrl <- run_HWE(zkv_data, 0, "rs1143623")
HWEex_rs944_case <- run_HWE(zkv_data, 1, "rs16944")
HWEex_rs944_ctrl <- run_HWE(zkv_data, 0, "rs16944")
HWEex_rs627_case <- run_HWE(zkv_data, 1, "rs1143627")
HWEex_rs627_ctrl <- run_HWE(zkv_data, 0, "rs1143627")


#___ SNP association analysis ___#
## Allele Frequency

#' Allelic Association Test
#'
#' Calculates allele frequencies in cases and controls and performs a chi-square test to assess association between allele counts and phenotype
#'
#' @param dataset A data frame containing SNP genotypes and phenotype.
#' @param case Value used to identify the case group (e.g., 1).
#' @param control Value used to identify the control group (e.g., 0).
#' @param snp The name of the SNP column (e.g., "rs4848306").
#' @param allele1 The first allele to be counted (e.g., "A").
#' @param allele2 The second allele to be counted (e.g., "G").
#'
#' @return Prints the allele count matrix and results of the chi-square test.
#' @export
run_assoc_allele <- function(dataset, case, control, snp, allele1, allele2) {
  allele1_case <- sum(str_count(dataset[dataset$Group == case, snp], allele1))
  allele1_control <- sum(str_count(dataset[dataset$Group == control, snp], allele1))
  allele2_case <- sum(str_count(dataset[dataset$Group == case, snp], allele2))
  allele2_control <- sum(str_count(dataset[dataset$Group == control, snp], allele2))
  
  alleles <- matrix(c(allele1_case, allele2_case, allele1_control, allele2_control), nrow = 2, byrow = TRUE,
                    dimnames = list(Group = c("Case", "Control"), Allele = c(allele1, allele2)))
  
  assoc_alelos_test <- chisq.test(alleles)
  print(alleles)
  print(assoc_alelos_test)
}
  
assoc_alelos_rs306 <- run_assoc_allele(zkv_data, 1, 0, "rs4848306", "A", "G")
assoc_alelos_rs623 <- run_assoc_allele(zkv_data, 1, 0, "rs1143623", "C", "G")
assoc_alelos_rs944 <- run_assoc_allele(zkv_data, 1, 0, "rs16944", "A", "G")
assoc_alelos_rs627 <- run_assoc_allele(zkv_data, 1, 0, "rs1143627", "A", "G")

## Genotype Frequency
snpdata <- setupSNP(data = zkv_data, 
                    colSNPs = 20:23,  
                    sep = "",
                    info = metadata_columns)

snpdata1 <- snpdata
snpdata1$rs16944  <- reorder(snpdata1$rs16944,  "minor")
snpdata1$rs1143627  <- reorder(snpdata1$rs1143627,  "minor")

association(Group ~ rs4848306, data = snpdata)
association(Group ~ rs1143623, data = snpdata)
association(Group ~ rs16944, data = snpdata1)
association(Group ~ rs1143627, data = snpdata1)


#___ Haplotype analysis ___#
snpsH <- c("rs4848306", "rs1143623",  "rs16944", "rs1143627")
genoH <- make.geno(snpdata, snpsH)
print(genoH)

haplo_cc_results <- haplo.cc(
  y = snpdata$Group,
  geno = genoH,
  x.adj = NULL,
  locus.label = snpsH,
  ci.prob = 0.95,
  simulate = FALSE)

print(haplo_cc_results)

## Additional test
em <- haplo.em(genoH, locus.label = snpsH)

contingency_table <- haplo.em(genoH, locus.label = snpsH) %>% 
  summary() %>% 
  as.data.frame() %>% 
  filter(posterior > 0.90) %>% 
  mutate(ID = snpdata$Sample,
         Group = snpdata$Group,
         hap4_count = (hap1 == 4) + (hap2 == 4),
         others = 2 - hap4_count) %>% 
  group_by(Group) %>% 
  summarise(Hap4 = sum(hap4_count),
            Outros = sum(others)) %>% 
  column_to_rownames("Group") 
print(contingency_table)

additional_hap_test <- chisq.test(contingency_table)
print(additional_hap_test)

#___ Compare frequencies to ABraOM Database ___# 

#' Compare Allele Frequencies with ABraOM Database
#'
#' Compares observed allele frequencies in cases and controls with expected frequencies from the ABraOM database using chi-square test for given probabilities.
#'
#' @param dataset A data frame containing SNP genotype.
#' @param snp The name of the SNP column (e.g., "rs4848306").
#' @param allele1 The first allele (e.g., "G").
#' @param allele2 The second allele (e.g., "A").
#' @param freq_db_allele1 Expected frequency of the first allele from ABraOM.
#' @param freq_db_allele2 Expected frequency of the second allele from ABraOM.
#'
#' @return Prints allele counts for cases and controls, and chi-square test results.
#' @export
run_ABraOM <- function(dataset, snp, allele1, allele2, freq_db_allele1, freq_db_allele2) {
  allele1_case <- sum(str_count(dataset[dataset$Group == 1, snp], allele1))
  allele2_case <- sum(str_count(dataset[dataset$Group == 1, snp], allele2))
  allele1_control <- sum(str_count(dataset[dataset$Group == 0, snp], allele1))
  allele2_control <- sum(str_count(dataset[dataset$Group == 0, snp], allele2))
  
  print(paste0("Case: ", allele1, "=", allele1_case, "; ", allele2, "=", allele2_case))
  case_test <- chisq.test(c(allele1_case, allele2_case), p = c(freq_db_allele1, freq_db_allele2))
  print(case_test)
  
  print(paste0("Control: ", allele1, "=", allele1_control, "; ", allele2, "=", allele2_control))
  control_test <- chisq.test(c(allele1_control, allele2_control), p = c(freq_db_allele1, freq_db_allele2))
  print(control_test)
}

rs306_ABraOM <- run_ABraOM(zkv_data, "rs4848306", "G", "A", 0.5918, 0.4082)
rs623_ABraOM <- run_ABraOM(zkv_data, "rs1143623", "C", "G", 0.7447, 0.2553)
rs944_ABraOM <- run_ABraOM(zkv_data, "rs16944", "A", "G", 0.4108, 0.5892)
rs627_ABraOM <- run_ABraOM(zkv_data, "rs1143627", "A", "G", 0.5734, 0.4266)


#___ Linkage Desequilibrium Analysis ___#
## Genetic R package (calculates p-val)
zkv_gt_snps <- zkv_genotypes %>%
  dplyr::select(-Sample)

zkv_gt_genetics <- zkv_gt_snps %>%
  mutate(rs1143627 = as.genotype(rs1143627, sep = ""),
         rs16944 = as.genotype(rs16944, sep = ""),
         rs1143623 = as.genotype(rs1143623, sep = ""),
         rs4848306 = as.genotype(rs4848306, sep = ""))

ld_genetics <- LD(zkv_gt_genetics)
print(ld_genetics)

## snpStats R package (does not calculates p-val)
zkv_gt_snstat <-  zkv_gt_snps %>%
  mutate(
    rs1143627 = case_when(
      rs1143627 == "GG" ~ 2,
      rs1143627 == "AG" ~ 1,
      rs1143627 == "AA" ~ 0),
    
    rs16944 = case_when(
      rs16944 == "AA" ~ 0,
      rs16944 == "AG" ~ 1,
      rs16944 == "GG" ~ 2),
    
    rs1143623 = case_when(
      rs1143623 == "CC" ~ 2,
      rs1143623 == "CG" ~ 1,
      rs1143623 == "GG" ~ 0),
    
    rs4848306 = case_when(
      rs4848306 == "GG" ~ 0,
      rs4848306 == "AG" ~ 1,
      rs4848306 == "AA" ~ 2))

snp_matrix <- as(as.matrix(zkv_gt_snstat), "SnpMatrix")
ld_matrix <- ld(snp_matrix, stats = c("D.prime", "R.squared"), depth = 3)
print(ld_matrix)

LDheatmap(snp_matrix, 
          genetic.distances = c(-31, -511, -1464, -3737),
          LDmeasure = "r",
          SNP.name = c("rs1143627", "rs16944", "rs1143623", "rs4848306"),
          distances = "physical",
          add.map = TRUE,
          color = magma(20),
          title = NULL,
          geneMapLabelY = 0.2,
          name="ldheatmap")

grid.edit(gPath("ldheatmap", "geneMap","SNPnames"), gp = gpar(cex=1, col="black"))

