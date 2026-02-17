library(stringr)
library(QuantQC)
library(arrow)


################
# For this analysis, it assumes for some of the file paths that you are working our 
# of the folder this R script is in for the github dir
################


## Making linker, we will need to edit this manually because we also need to add the
## plate colum (should be 1)

## We also didnt use the best string extract for the well so we missed
## ones with reinjections, we will need to fix this as well

# Fix to local path
path_raw <- '/Users/andrewleduc/Desktop/jmod_raw/'

Raw_data_ <- arrow::read_parquet(
  paste0(path_raw,'/JModPlate1Slide2.parquet')
)

# Extract well from the 
Raw_data_$well <- str_extract(Raw_data_$file_name, "(?<=_)[^_/]+(?=\\.mzML$)")

Raw_data_2 <- Raw_data_ %>% distinct(file_name,.keep_all = T)

Raw_data_2 <- Raw_data_2 %>% dplyr::select(file_name,well)
write.csv(Raw_data_2,'/Users/andrewleduc/Downloads/linker_jmod.csv')


#### Required libraries
library(QuantQC)
library(dplyr)
library(ggpubr)
library(Seurat)
library(arrow)


# Fix to local path
path_raw <- '/Users/andrewleduc/Desktop/jmod_raw/'
path_meta <- paste0(getwd(),'/../meta_data/')

#### paths for relevant raw and meta data
data_path <- paste0(path_raw,'JModPlate1Slide2.parquet')

linker_path <- paste0(path_meta,"linker_jmod.csv")


one <- paste0(path_meta,"sample_1_try2_isolated.xls")
two <- paste0(path_meta,"sample_2_isolated.xls")
all_cells <- list(mouse1 = one,
                  mouse2 = two)


#### Initialize QQC object and convert report table to peptide matrix
#### with linked meta data
r1_5day_male <- JMOD_to_QQC(data_path,linker_path, plex = 3, carrier = F)

r1_5day_male <- Miceotope_cellXpeptide_Jmod(r1_5day_male,chQVal = 1,t = 5)


r1_5day_male <- link_cellenONE_Raw(r1_5day_male,all_cells)



#### Some QC checks for LC batch effects and nPOP slide viz (optional)

#r1_5day_male <- Calculate_run_order_statistics(r1_5day_male)

# Plot MS Intensity Statistics
#PlotIntensityDrift(r1_5day_male)

#Plot Retention time Statistics
#PlotRTDrift(r1_5day_male)

PlotSlideLayout_celltype(r1_5day_male)
PlotSlideLayout_label(r1_5day_male)


#### Evaluating single cell quality compared to negative controls
#### and filtering out of failed cells

r1_5day_male <- EvaluateNegativeControls(r1_5day_male)

PlotNegCtrl(r1_5day_male)

# log10 intensity as failed cell filter
r1_5day_male <- FilterBadCells(r1_5day_male, min_intens = 9)



#### Remove extra peptides that are lowly abundant for proteins that already have
#### 5 good peptides

# for abundance
nrow(r1_5day_male@matricies@peptide)
r1_5day_male <- Trim_extra_peptides(r1_5day_male)
nrow(r1_5day_male@matricies@peptide)

# for heavy and light peptides
nrow(r1_5day_male@miceotopes@Raw_H)
r1_5day_male <- Trim_extra_peptides_miceotopes(r1_5day_male)
nrow(r1_5day_male@miceotopes@Raw_H)




#### Plot cell size vs MS intensity
PlotCellSizeVsIntensity(r1_5day_male, type = 'sample')



#### Normalize peptide data and collapse to protein level
r1_5day_male@ms_type <- 'miceotopes'
r1_5day_male <- CollapseToProtein(r1_5day_male, 1,LC_correct = T)



# Viz for LC related batch effects
PlotLC_Deviations(r1_5day_male,global_trend = F)



#write.csv(r1_5day_male@matricies@protein,paste0(getwd(),'/../data_out_temp/R_pipeline_nonorm_protein.csv'))



#### Plot number of proteins peptides and protein level data missingness
r1_5day_male@ms_type <- 'DDA'
PlotProtAndPep(r1_5day_male)
PlotDataComplete(r1_5day_male)



#### Compute and plot correlations between relative peptide levels for peptides
#### mapping to the same protein
r1_5day_male <- SharedPeptideCor(r1_5day_male)

PlotPepCor(r1_5day_male)
median(r1_5day_male@pep.cor[[1]]$Cor)

#PlotMS1vMS2(r1_5day_male)



#### KNN protein level imputation
r1_5day_male <- KNN_impute(r1_5day_male)


#### Protein level batch correction for mTRAQ label bias

# Here we employ custom normalization to also control for samples (mice)

#r1_5day_male <- BatchCorrect(r1_5day_male,run = F,labels = T)

cellenONE_meta <- r1_5day_male@meta.data
protein_mat_imputed <- r1_5day_male@matricies@protein.imputed
protein_mat <- r1_5day_male@matricies@protein

# Get meta data for batch correction
batch_label  <- cellenONE_meta %>% dplyr::filter(ID %in% colnames(protein_mat_imputed))
batch_label <- batch_label[order(match(batch_label$ID,colnames(protein_mat_imputed))),]

sc.batch_cor <- limma::removeBatchEffect(r1_5day_male@matricies@protein.imputed,batch = batch_label$label, batch2 = batch_label$sample)

sc.batch_cor_noimp <- sc.batch_cor
sc.batch_cor_noimp[is.na(protein_mat)==T] <- NA

r1_5day_male@matricies@protein.imputed <- sc.batch_cor
r1_5day_male@matricies@protein <- sc.batch_cor_noimp


## Normalize out biases
prot_mat <- r1_5day_male@matricies@protein          # rows = proteins, cols = cells/samples
hb_id    <- "P01942"
hb       <- as.numeric(prot_mat[hb_id, ])

r2 <- rep(NA_real_, nrow(prot_mat))
adj <- prot_mat

for (i in seq_len(nrow(prot_mat))) {
  # skip hemoglobin itself
  #if (rownames(prot_mat)[i] == hb_id) next

  y <- as.numeric(prot_mat[i, ])

  ok <- is.finite(y) & is.finite(hb)
  if (sum(ok) < 3) next  # not enough points to fit a line

  fit <- lm(y[ok] ~ hb[ok])

  r2_i <- summary(fit)$r.squared
  r2[i] <- r2_i

  if (is.finite(r2_i) && r2_i > 0.05) {
    # residuals on the full vector (NAs where ok=FALSE)
    res <- rep(NA_real_, length(y))
    res[ok] <- resid(fit)

    # option A: pure residuals (centered around 0)
    # adj[i, ] <- res

    # option B (usually nicer for “abundance”): add back the original mean
    adj[i, ] <- res + mean(y[ok], na.rm = TRUE)
    adj[i, ] <- adj[i, ]-mean(adj[i, ],na.rm=T)
  }
}

r1_5day_male@matricies@protein <- adj



#### Compute PCA and overlay some cell type markers and potential artifacts
r1_5day_male <- ComputePCA(r1_5day_male, impute = F)



# Hemoglobin
FeaturePCA(r1_5day_male, prot = 'P01942', imputed = F)



# Other source variance?
FeaturePCA(r1_5day_male, prot = 'P12710', imputed = F)
FeaturePCA(r1_5day_male, prot = 'P00329', imputed = F)



# Portal
FeaturePCA(r1_5day_male, prot = 'P33267', imputed = F) #Cyp2f2
FeaturePCA(r1_5day_male, prot = 'Q61176', imputed = F) # Arg1
FeaturePCA(r1_5day_male, prot = 'Q91YI0', imputed = F) # Asl



# Central

plot(r1_5day_male@matricies@protein['Q05421',],r1_5day_male@matricies@protein['P33267',])
cor(r1_5day_male@matricies@protein['Q05421',],r1_5day_male@matricies@protein['P33267',],use = 'pairwise.complete.obs')


FeaturePCA(r1_5day_male, prot = 'Q05421', imputed = F) # Cyp2e1
FeaturePCA(r1_5day_male, prot = 'P15105', imputed = F) # Glul


# plot PCA options are "Run order" "Total protein" "Condition" "Label"
PlotPCA(r1_5day_male, by = "Condition")
PlotPCA(r1_5day_male, by = "Run order")
PlotPCA(r1_5day_male, by = "Total protein")
PlotPCA(r1_5day_male, by = "Label")





#### Compute UMAP and overlay some cell type markers and potential artifacts
r1_5day_male <- ComputeUMAP(r1_5day_male)

PlotUMAP(r1_5day_male)
PlotUMAP(r1_5day_male, by = 'Total protein')
PlotUMAP(r1_5day_male, by = 'Run order')
PlotUMAP(r1_5day_male, by = 'Label')
PlotUMAP(r1_5day_male, by = 'Condition')

FeatureUMAP(r1_5day_male,prot = 'P33267')



# save whole R data file
save(r1_5day_male, file = paste0(path_meta,"/../data_out_temp/r1_5day_male.RData"))

write.csv(r1_5day_male@matricies@protein,paste0(getwd(),'/../data_out_temp/R_pipeline_protein.csv'))
write.csv(r1_5day_male@miceotopes@Alpha_pep,paste0(getwd(),'/../data_out_temp/R_pipeline_alpha_pep.csv'))





# test <- r1_5day_male@meta.data %>% filter(sample == 'mouse1')
# tes2 <- r1_5day_male@meta.data %>% filter(sample == 'mouse2')
# 
# 
# 
# hist(r1_5day_male@miceotopes@HovL_pep['QVHPDTGISSK3',intersect(test$ID,colnames(r1_5day_male@miceotopes@HovL_pep))],30)
# 
# hold1 <- intersect(test$ID,colnames(r1_5day_male@miceotopes@HovL_pep))
# hold2 <- intersect(tes2$ID,colnames(r1_5day_male@miceotopes@HovL_pep))
# 
# 
# sect <- intersect(rownames(r1_5day_male2@miceotopes@Alpha_prot),rownames(r1_5day_male2@matricies@protein_abs))
# 
# 
# cor(rowMeans(log2(r1_5day_male2@miceotopes@Alpha_prot[sect,]),na.rm=T),
#      rowMeans(log10(r1_5day_male2@matricies@protein_abs[sect,]),na.rm=T),
#     use = 'pairwise.complete.obs',method = 'spearman')
# plot(rowMeans(log2(r1_5day_male2@miceotopes@Alpha_prot[sect,hold1]*5/3),na.rm=T),
#     rowMeans(log2(r1_5day_male2@matricies@protein_abs[sect,hold2]),na.rm=T))
# 
# 
# 
# r1_5day_male <- Miceotope_protein_collapse(r1_5day_male)
# 
# r1_5day_male@miceotopes@Alpha_prot[,hold1] <- r1_5day_male@miceotopes@Alpha_prot[,hold1]*5/3
# 
# 
# plot(rowMeans(log2(r1_5day_male@miceotopes@Alpha_prot[,hold1]),na.rm=T),
#      rowMeans(log2(r1_5day_male@miceotopes@Alpha_prot[,hold2]),na.rm=T))
# abline(a=0,b=1)
# 
# sect <- intersect(rownames(r1_5day_male@miceotopes@Alpha_prot),rownames(r1_5day_male@matricies@protein))
# 
# 
# deg_norm <- QuantQC::normalize(r1_5day_male@miceotopes@Alpha_prot,log = T)
# deg_norm <- r1_5day_male@miceotopes@Alpha_prot
# 
# 
# 
# library(psych)
# cors_deg <- c()
# cors_axis <- c()
# cor_deg_axis <- c()
# numb_dp <- c()
# for(i in 1:length(sect)){
# 
#   cors_deg <- c(cors_deg,cor(deg_norm[sect[i],],r1_5day_male2@matricies@protein[sect[i],],use= 'pairwise.complete.obs'))
#   cors_axis <- c(cors_axis,cor(r1_5day_male2@matricies@protein[sect[i],],r1_5day_male2@matricies@protein['Q05421',],use= 'pairwise.complete.obs'))
# 
#   cor_deg_axis <- c(cor_deg_axis,cor(deg_norm[sect[i],],r1_5day_male2@matricies@protein['Q05421',],use= 'pairwise.complete.obs'))
# 
#   numb_dp <- c(numb_dp,pairwiseCount(deg_norm[sect[i],],r1_5day_male2@matricies@protein[sect[i],])[[1]])
# 
# }
# 
# df_deg_cor <- data.frame(prot = sect, cor_deg = cors_deg,cor_deg_axis=cor_deg_axis,cors_axis = cors_axis,numb = numb_dp)
# df_deg_cor <- df_deg_cor %>% filter(numb > 50)
# 
# 
# 
# ## Use markers to make pseudo spatial axis, look over Nature article to get more markers for this
# 
# ## Order cells by axis values
# 
# ## Plot Abundance and deg rate vs Axis
# 
# ## Try to denoise the spatial axis by binning cells
# 
# 
# 
# ## Compare to mRNA from nature paper
# 
# ## Bin these comparisons by both absolute half life of protein, and cor to relative half life
# 
# 
# 
# ## functional groupings of proteins
# 
# 
# 
# 
# QuantQC::Mice_DimPlot_turnover(r1_5day_male,reuct = 'UMAP',by = 'Total')
# 
# Mice_DimPlot_turnover()
# 
# r1_5day_male@meta.data
# 









