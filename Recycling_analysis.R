
# sdf


LF_data <- read_parquet('/Users/andrewleduc/Desktop/Github/Kevin_PTI/report.parquet')
LF_data$Run[LF_data$Run == '2026-01-14_2ngL2_20Th-10ms_28min_nce30'] <- 'L2'
LF_data$Run[LF_data$Run == '2026-01-14_2ngL1_20Th-10ms_28min_nce30'] <- 'L1'
LF_data$Run


LF_data <- LF_data %>% filter(Ms1.Area != 0)

sum(LF_data$Run ==unique(LF_data$Run)[1])
sum(LF_data$Run ==unique(LF_data$Run)[2])

LF_data$seqcharge <- paste0(LF_data$Stripped.Sequence,LF_data$Precursor.Charge)

# Intensity_check <- LF_data %>% group_by(seqcharge,Run) %>% summarise(Ms1.Area = median(Ms1.Area,na.rm=T))
# Intensity_check <- reshape2::dcast(Intensity_check,seqcharge ~ Run, value.var = 'Ms1.Area')
# Intensity_check$seqcharge <- NULL
# Intensity_check <- as.matrix(Intensity_check)
# for(i in 1:nrow(Intensity_check)){
#   Intensity_check[i,] <- Intensity_check[i,]/mean(Intensity_check[i,])
# }
#
# Intensity_check <- Intensity_check[complete.cases(Intensity_check),]
# Intensity_check <- melt(Intensity_check)
# ggplot(Intensity_check,aes(x = Var2, y = log2(value))) + geom_boxplot() + ylab('log2 intensity intersected')
# median((Intensity_check$value)[Intensity_check$Var2== unique(Intensity_check$Var2)[1]])/
#   median((Intensity_check$value)[Intensity_check$Var2== unique(Intensity_check$Var2)[2]])

LF_data$kcount <- str_count(LF_data$Stripped.Sequence, "K")
LF_data <- LF_data %>% filter(kcount == 2)


LF_data$hcount <- str_count(LF_data$Precursor.Id,"K_6C13-6")





HL_dat <- LF_data %>% filter(hcount == 1)
HH_dat <- LF_data %>% filter(hcount == 2)

sect <- intersect(HL_dat$seqcharge,HH_dat$seqcharge)
HL_dat <- HL_dat %>% filter(seqcharge %in% sect)
HH_dat <- HH_dat %>% filter(seqcharge %in% sect)

HL_dat <- HL_dat %>% group_by(seqcharge,Run) %>% summarise(Precursor.Quantity = median(Precursor.Quantity,na.rm=T))

HL_dat <- reshape2::dcast(HL_dat,seqcharge ~ Run, value.var = 'Precursor.Quantity')
HH_dat <- reshape2::dcast(HH_dat,seqcharge ~ Run, value.var = 'Precursor.Quantity')
sect2 <- intersect(colnames(HL_dat),colnames(HH_dat))

HL_dat <- HL_dat[,sect2]
HH_dat <- HH_dat[,sect2]
rownames(HL_dat) <- HL_dat$seqcharge
rownames(HH_dat) <- HH_dat$seqcharge
HL_dat$seqcharge <- NULL
HH_dat$seqcharge <- NULL
HL_dat <- as.matrix(HL_dat)
HH_dat <- as.matrix(HH_dat)

pct_light <- (HL_dat/HH_dat)/(HL_dat/HH_dat+2) # Marko
pct_Heavy <- 2/(HL_dat/HH_dat+2)

tot_I <- HL_dat+HH_dat

plot(pct_Heavy[,2],log(tot_I[,2]))
plot(pct_Heavy[,1],log(tot_I[,1]))

colnames(pct_Heavy)


pct_Heavy2 <- pct_Heavy[complete.cases(pct_Heavy),]
tot_I2 <- tot_I[complete.cases(tot_I),]
plot(pct_Heavy2[,2],log(tot_I2[,2]))
plot(pct_Heavy2[,1],log(tot_I2[,1]))

pct_Heavy_plot <- melt(pct_Heavy2)

ggplot(pct_Heavy_plot,aes(x = Var2,y = value)) + geom_boxplot()+
  xlab('') + ylab('% Heavy') + theme_bw(base_size = 18)


median(pct_Heavy2[,1]) #t=3
median(pct_Heavy2[,2]) # t=5


# IC for guess
a_init <- log(2)/0.7
b_init <- log(2)/20

time_obs <- c(0, 3, 5)
recycle_value_obs <- c(0,.51, .54)


# Use nls() to fit the model and tune a and b
fit <- nls(recycle_value_obs ~ 1 - .5 * exp(-a * time_obs) - .5 * exp(-b * time_obs),
           start = list(a = a_init, b = b_init),
           control = nls.control(maxiter = 1000),
           algorithm = "port")



# Extract the fitted parameters
params <- coef(fit)
a_fitted <- as.numeric(params["a"])
b_fitted <- as.numeric(params["b"])


t <- seq(0,20,by = .1)

L_aa_fitted <-  .5 * exp(-a_fitted * t) + .5 * exp(-b_fitted * t)
plot(t,1-L_aa_fitted)


#### Required libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(QuantQC)
library(stringr)
library(foreach)
library(doParallel)
library(reshape2)

library(lme4)
library(splines)
library(gridExtra)
library(mgcv)




##################
# raw data path
##################

path_raw <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/raw/02_Recycle_searches/'

path_dat <- '/Users/andrewleduc/Desktop/Github/Miceotoptes_single_cell/dat/'


##################
# Functions for recycle correction analysis
##################
TLS <- function(vect1,vect2){

  vect1[vect1 == -Inf] <- NA
  vect1[vect1 == Inf] <- NA



  vect2[vect2 == -Inf] <- NA
  vect2[vect2 == Inf] <- NA

  int_x <- mean(vect1,na.rm=T)
  int_y <- mean(vect2,na.rm=T)

  vect1 <- vect1-int_x
  vect2 <- vect2-int_y

  mat <- cbind(vect1,vect2)

  mat <- mat[complete.cases(mat),]

  TLS_mat <- svd(mat)$v

  slope <- TLS_mat[1,1]/TLS_mat[2,1]

  int <- c(int_x,int_y)

  return(list(slope,int))

}

NA_false_ratios = function(ratio_mat, t, anno, recycle_params,prep){

  recycle_params <- recycle_params %>% filter(Time == t)

  colnames(ratio_mat) <- paste0(colnames(ratio_mat),prep)

  ratio_mat <- ratio_mat[,intersect(colnames(ratio_mat),meta_data$ID)]

  for(i in 1:ncol(ratio_mat)){

    meta_temp <- meta_data %>% filter(ID == colnames(ratio_mat)[i])
    Final_values_temp <- recycle_params %>% filter(age == meta_temp$age)
    if(meta_temp$Cell_Type %in% c('Basal','Secratory','Cilliated','Immune')){
      Final_values_temp <- Final_values_temp %>% filter(tissue_save =='Epithelial')
    }else{
      Final_values_temp <- Final_values_temp %>% filter(tissue_save =='NotEpithelial')
    }

    for (j in 1:nrow(ratio_mat)) {

      # Get the current L and L0 values
      L_measured <- (ratio_mat[j, i])

      if(is.na(L_measured) == F){
        if(L_measured < Final_values_temp$avg){
          ratio_mat[j, i] <- NA
        }


      }
    }

  }
  return(ratio_mat)
}

# Function that solves for Kr (alpha is variable used) math is in methods or can read from function
L_function <- function(alpha, L0, a, b, t) {
  c = .5
  beta = L0*alpha
  #Lt <- (exp(-(a*t) - b*t)*(-(a*beta*c*exp(a*t)) - b*beta*c*exp(b*t) +
  #                            a*beta*c*exp(a*t + b*t) + b*beta*c*exp(a*t + b*t) + a*b*exp(a*t + b*t)*L0))/(a*b)


  Lt <- (a*beta*c - 2*alpha*beta*c + b*beta*c + alpha*beta*c*exp((-a + alpha)*t) - b*beta*c*exp((-a + alpha)*t) -
           a*beta*c*exp((alpha - b)*t) + alpha*beta*c*exp((alpha - b)*t) - a*alpha*L0 + alpha^2*L0 + a*b*L0 - alpha*b*L0)/
    ((-a + alpha)*(alpha - b)*exp(alpha*t))

  return(Lt)
}

# RSS objective function for optim to solve for Kr
objective_function <- function(alpha, L_measured, L0, a, b, t) {

  L_function <- function(alpha, L0, a, b, t) {
    c = .5
    beta = L0*alpha
    #Lt <- (exp(-(a*t) - b*t)*(-(a*beta*c*exp(a*t)) - b*beta*c*exp(b*t) +
    #                            a*beta*c*exp(a*t + b*t) + b*beta*c*exp(a*t + b*t) + a*b*exp(a*t + b*t)*L0))/(a*b)


    Lt <- (a*beta*c - 2*alpha*beta*c + b*beta*c + alpha*beta*c*exp((-a + alpha)*t) - b*beta*c*exp((-a + alpha)*t) -
             a*beta*c*exp((alpha - b)*t) + alpha*beta*c*exp((alpha - b)*t) - a*alpha*L0 + alpha^2*L0 + a*b*L0 - alpha*b*L0)/
      ((-a + alpha)*(alpha - b)*exp(alpha*t))

    return(Lt)
  }

  L_predicted <- L_function(alpha, L0, a, b, t)
  residuals <- L_predicted - (L_measured)
  return(sum(residuals^2))
}

NA_false_ratios = function(ratio_mat, t, anno, recycle_params,prep){

  recycle_params <- recycle_params %>% filter(Time == t)

  colnames(ratio_mat) <- paste0(colnames(ratio_mat),prep)

  ratio_mat <- ratio_mat[,intersect(colnames(ratio_mat),meta_data$ID)]

  for(i in 1:ncol(ratio_mat)){

    meta_temp <- anno %>% filter(ID == colnames(ratio_mat)[i])
    Final_values_temp <- recycle_params %>% filter(age == meta_temp$age)

    if(meta_temp$Cell_Type %in% c('Basal','Secratory','Cilliated','Immune')){
      Final_values_temp <- Final_values_temp %>% filter(Recycle_type =='Epithelial')
    }else{
      Final_values_temp <- Final_values_temp %>% filter(Recycle_type =='NotEpithelial')
    }

    for (j in 1:nrow(ratio_mat)) {

      # Get the current L and L0 values
      L_measured <- (ratio_mat[j, i])

      if(is.na(L_measured) == F){

        if((1-L_measured) > Final_values_temp$avg){
          ratio_mat[j, i] <- NA
        }


      }
    }

  }
  return(ratio_mat)
}

Recycle_adjust_par_bar <- function(L_mat, L0_mat, t, recycle_params, objective_function, ncores = 2) {
  # Load required packages
  library(dplyr)
  library(foreach)
  library(doSNOW)


  # Subset to the common columns in meta_data


  # Prepare an output matrix with appropriate dimnames
  out_mat <- matrix(NA, nrow = nrow(L_mat), ncol = ncol(L_mat),
                    dimnames = list(rownames(L_mat), colnames(L_mat)))

  # Set up a SNOW cluster and register it
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)

  # Set up a progress bar that tracks progress over columns
  pb <- txtProgressBar(max = ncol(L_mat), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  # Parallelize the outer loop over columns using foreach:
  # Each iteration returns a vector of optimized alpha values for one column.
  result_list <- foreach(i = seq_len(ncol(L_mat)),
                         .packages = "dplyr",
                         .options.snow = opts) %dopar% {
                           col_res <- rep(NA, nrow(L_mat))  # initialize result vector for column i
                           current_ID <- colnames(L_mat)[i]



                           # Process each row (peptide) in the current column
                           for (j in seq_len(nrow(L_mat))) {
                             L_measured <- L_mat[j, i]
                             L0_value   <- L0_mat[j, i]
                             if (!is.na(L_measured) && !is.na(L0_value)) {
                               # Run the 1D optimization (using L-BFGS-B here)
                               result <- optim(
                                 par = 0.3,
                                 fn = objective_function,
                                 L_measured = L_measured,
                                 L0 = L0_value,
                                 a = recycle_params[1],
                                 b = recycle_params[2],
                                 t = t,
                                 method = "L-BFGS-B",
                                 lower = 0,
                                 upper = 10
                               )
                               col_res[j] <- result$par
                             }
                           }
                           # Return the vector for this column
                           col_res
                         }

  # Close the progress bar and stop the cluster
  close(pb)
  stopCluster(cl)

  # Combine the results from each column into a matrix.
  # result_list is a list of column vectors; cbind them together.
  out_mat <- do.call(cbind, result_list)
  colnames(out_mat) <- colnames(L_mat)
  rownames(out_mat) <- rownames(L_mat)

  return(out_mat)
}




##################
# Apply recycling corrections to single cells
##################


meta = r1_5day_male@meta.data
#mouse 3 day

days_3 <- meta %>% filter(sample == 'mouse1')
days_3 <- days_3 %>% filter(ID %in% colnames(r1_5day_male@miceotopes@Raw_L))

#mouse 5 day

days_5 <- meta %>% filter(sample == 'mouse2')
days_5 <- days_5 %>% filter(ID %in% colnames(r1_5day_male@miceotopes@Raw_L))

r1 <- r1_5day_male@miceotopes@Raw_L[,days_3$ID]/(r1_5day_male@miceotopes@Raw_L[,days_3$ID]+r1_5day_male@miceotopes@Raw_H[,days_3$ID])
r2 <- r1_5day_male@miceotopes@Raw_L[,days_5$ID]/(r1_5day_male@miceotopes@Raw_L[,days_5$ID]+r1_5day_male@miceotopes@Raw_H[,days_5$ID])


#na_r1 <- NA_false_ratios(r1,3,meta_data,Final_values,'_prep1')
#na_r2 <- NA_false_ratios(r2,5,meta_data,Final_values,'_prep1')




adj_r1 <- Recycle_adjust_par_bar(r1_5day_male@miceotopes@Raw_L[,days_3$ID],(r1_5day_male@miceotopes@Raw_L[,days_3$ID]+r1_5day_male@miceotopes@Raw_H[,days_3$ID])
                                 ,3,meta_data,c(a_fitted,b_fitted),objective_function)

adj_r1[adj_r1==10] <- NA

adj_r2 <- Recycle_adjust_par_bar(r1_5day_male@miceotopes@Raw_L[,days_5$ID],(r1_5day_male@miceotopes@Raw_L[,days_5$ID]+r1_5day_male@miceotopes@Raw_H[,days_5$ID])
                                 ,5,meta_data,c(a_fitted,b_fitted),objective_function)
adj_r2[adj_r2==10] <- NA




#### Looking at not corrected values for comparison
r1 <- -log(r1[,colnames(na_r1)])/5
r4 <- -log(r4[,colnames(na_r4)])/10

median(log(2)/r1,na.rm=T)
median(log(2)/r2,na.rm=T)

median((log(2)/adj_r1),na.rm=T)
median((log(2)/adj_r2),na.rm=T)


adj_r1_ <- adj_r1[rowSums(is.na(adj_r1)==F) > 20,]
adj_r4_ <- adj_r4[rowSums(is.na(adj_r4)==F) > 20,]

ss <- intersect(rownames(adj_r1_),rownames(adj_r2_))


adj_r1[adj_r1==0] <- NA
adj_r2[adj_r2==0] <- NA


#### Compare rates for 5 and 10 day mice with corrected data
df_corrected_compare <- data.frame(
  logx = rowMeans(log2(cbind(adj_r2[ss,young_p2$ID],adj_r1[ss,young_p1$ID])), na.rm = TRUE),
  logy = rowMeans(log2(cbind(adj_r3[ss,young_p3$ID],adj_r4[ss,young_p4$ID])), na.rm = TRUE)
)

cor(df_corrected_compare$logx,df_corrected_compare$logy,use = 'pairwise.complete.obs')

df_corrected_compare$X <- 2 ^ df_corrected_compare$logx
df_corrected_compare$Y <- 2 ^ df_corrected_compare$logy

ggplot(df_corrected_compare, aes(X, Y)) +
  geom_point(alpha = .5, size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  geom_abline(intercept = -.07, slope = 1/1.1,color = 'red')+
  scale_x_log10()+scale_y_log10()+
  coord_cartesian(xlim = c(.05,1),ylim = c(.05,1))+
  labs(x = "5 day degradation rate, Days^-1",
       y = "10 day degradation rate, Days^-1") +
  theme_classic(base_size = 15) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.line     = element_blank()
  )

1/TLS(df$logx,df$logy)[[1]]


#### Compare rates for 5 and 10 day mice without correction applied for reference
df_uncor <- data.frame(
  logx = rowMeans(log2(cbind(r1[ss,young_p1$ID],r2[ss,young_p2$ID])), na.rm = TRUE),
  logy = rowMeans(log2(cbind(r3[ss,young_p3$ID],r4[ss,young_p4$ID])), na.rm = TRUE)
)

cor(df_uncor$logx,df_uncor$logy,use = 'pairwise.complete.obs')

df_uncor$X <- 2 ^ df_uncor$logx
df_uncor$Y <- 2 ^ df_uncor$logy

ggplot(df_uncor, aes(X, Y)) +
  geom_point(alpha = .5, size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  geom_abline(intercept = -.38, slope = 1/1.6, linetype = "solid",color = 'red')+
  scale_x_log10()+scale_y_log10()+
  coord_cartesian(xlim = c(.03,.2),ylim = c(.03,.2))+
  labs(x = "5 day degradation rate, Days^-1",
       y = "10 day degradation rate, Days^-1") +
  theme_classic(base_size = 15) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # add frame
    axis.line     = element_blank()                                       # hide dup. axes
  )

1/TLS(df$logx,df$logy)[[1]]


#### Compare magnitudes of corrected and uncorrected half life

df_compare <- data.frame(
  uncorrected = rowMeans(r1[ss,], na.rm = TRUE),
  corrected = rowMeans(adj_r1[ss,], na.rm = TRUE)
)
df_compare <- melt(df_compare)

ggplot(df_compare, aes(y = log(2)/value, fill = variable)) +
  geom_histogram(alpha = 0.5, bins = 30, colour = "black",
                 orientation = "y") +           # <- horizontal bars
  facet_wrap(~ variable, ncol = 2) +
  labs(
    x = "Count in each bin",
    y = "Value",
    title = ""
  ) +
  dot_plot +
  theme(legend.position = "none")+
  scale_y_log10() + ylab('Clearance half life (days)')+ xlab('# of proteins')+
  scale_fill_manual(values = c('grey80','black'))


#### Compare old and young clearance rates

df_old_vs_young <- data.frame(
  logx = rowMeans((cbind(adj_r1[ss,young_p1$ID])), na.rm = TRUE),
  logy = rowMeans((cbind(adj_r1[ss,old_p1$ID])), na.rm = TRUE)
)

cor(df_old_vs_young$logx,df_old_vs_young$logy,use = 'pairwise.complete.obs')

ggplot(df_old_vs_young, aes(logx, logy)) +
  geom_point(alpha = .5, size = 1.2) +
  geom_abline(intercept = 0,  slope = 1,  linetype = "dashed") +
  geom_abline(intercept = -.1, slope = 1,  colour = "red") +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = c(.05, 1), ylim = c(.05, 1)) +
  labs(x = "2 months", y = "24 months") +
  theme_classic(base_size = 15) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # add frame
    axis.line     = element_blank()                                       # hide dup. axes
  )



#### Save data in QQC objects


r1_5day_male@miceotopes@Alpha_pep <- adj_r1
r2_5day_female@miceotopes@Alpha_pep <- adj_r2


r1_5day_male <- Miceotope_protein_collapse(r1_5day_male)
r2_5day_female <- Miceotope_protein_collapse(r2_5day_female)


##################
# Save updated QQC objects
##################

save(r1_5day_male, file = paste0(path_dat,'03_QuantQC_objects/r1_5day_male.RData'))
save(r2_5day_female, file = paste0(path_dat,'03_QuantQC_objects/r2_5day_female.RData'))

##################
# Save Single cell matricies for deg
##################

p1_alpha <-  QuantQC::Normalize_reference_vector(r1_5day_male@miceotopes@Alpha_prot, log = T)
p2_alpha <-  QuantQC::Normalize_reference_vector(r2_5day_female@miceotopes@Alpha_prot , log = T)

p1_alpha_abs <-  r1_5day_male@miceotopes@Alpha_prot
p2_alpha_abs <-  r2_5day_female@miceotopes@Alpha_prot


### Save deg rates for all IDed proteins
colnames(p1_alpha) <- paste0(colnames(p1_alpha),'_prep1')
colnames(p2_alpha) <- paste0(colnames(p2_alpha),'_prep2')
colnames(p3_alpha) <- paste0(colnames(p3_alpha),'_prep3')
colnames(p4_alpha) <- paste0(colnames(p4_alpha),'_prep4')
p1_alpha <- as.data.frame(p1_alpha)
p1_alpha$prot <- rownames(p1_alpha)
p2_alpha <- as.data.frame(p2_alpha)
p2_alpha$prot <- rownames(p2_alpha)
p3_alpha <- as.data.frame(p3_alpha)
p3_alpha$prot <- rownames(p3_alpha)
p4_alpha <- as.data.frame(p4_alpha)
p4_alpha$prot <- rownames(p4_alpha)
p_1and2 <- p1_alpha %>% merge(p2_alpha, by = 'prot',all = TRUE)
p_123 <- p_1and2 %>% merge(p3_alpha, by = 'prot',all = TRUE)
p_all_alpha <- p_123 %>% merge(p4_alpha, by = 'prot',all = TRUE)
rownames(p_all_alpha) <- p_all_alpha$prot
p_all_alpha$prot <- NULL
p_all_alpha <- as.matrix(p_all_alpha)

write.csv(p_all_alpha,paste0(path_dat,'04_Gene_X_SingleCell_and_annotations/clearance_relative.csv'))


colnames(p1_alpha_abs) <- paste0(colnames(p1_alpha_abs),'_prep1')
colnames(p2_alpha_abs) <- paste0(colnames(p2_alpha_abs),'_prep2')
colnames(p3_alpha_abs) <- paste0(colnames(p3_alpha_abs),'_prep3')
colnames(p4_alpha_abs) <- paste0(colnames(p4_alpha_abs),'_prep4')
p1_alpha_abs <- as.data.frame(p1_alpha_abs)
p1_alpha_abs$prot <- rownames(p1_alpha_abs)
p2_alpha_abs <- as.data.frame(p2_alpha_abs)
p2_alpha_abs$prot <- rownames(p2_alpha_abs)
p3_alpha_abs <- as.data.frame(p3_alpha_abs)
p3_alpha_abs$prot <- rownames(p3_alpha_abs)
p4_alpha_abs <- as.data.frame(p4_alpha_abs)
p4_alpha_abs$prot <- rownames(p4_alpha_abs)
p_1and2 <- p1_alpha_abs %>% merge(p2_alpha_abs, by = 'prot',all = TRUE)
p_123 <- p_1and2 %>% merge(p3_alpha_abs, by = 'prot',all = TRUE)
p_all_alpha_abs <- p_123 %>% merge(p4_alpha_abs, by = 'prot',all = TRUE)
rownames(p_all_alpha_abs) <- p_all_alpha_abs$prot
p_all_alpha_abs$prot <- NULL
p_all_alpha_abs <- as.matrix(p_all_alpha_abs)

write.csv(p_all_alpha_abs,paste0(path_dat,'04_Gene_X_SingleCell_and_annotations/clearance_absolute.csv'))
























