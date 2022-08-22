#################
# Load packages #
#################

options(warn=-1)
library(caret)
library(cowplot)
library(patchwork)
library(iml)
library(tidyverse)
library(ggembl)
library(vegan)
#library(ggpubr)
#library(ggsignif)
library(mlr3)
library(mlr3learners)
library("mlr3tuning")
library(ranger)
library(mlr3viz)
library(pROC)
library("paradox")
library('missForest')
library(scales)

#############
# Load data #
#############

mpa2Profiles <- read_tsv("data/metaphlan_profiles.tsv", col_names = T)
engraftmentData <- read_tsv("data/engraftment.tsv", col_names = T)
metaFmts <- read_tsv("data/metadata_fmts_clean_stage_II.tsv", col_names = T)
metaFmts <- metaFmts %>% select(fmt_id, single_pre_sample, single_donor_sample, display_dataset)
metaSamples <- read_tsv("data/metadata_samples_clean_stage_II.tsv", col_names = T)
metaSamples <- metaSamples %>% select(dataset, sample_name, profile_id)

clinicalMeta <- read_csv("data/clinical_metadata_for_resub.csv")
# not used currently
#mdPhenotypes <- read_tsv("metadata/MD.phenotypes.tsv")
traitar <- read_tsv("data/traitar_output_646SGBs_tab.txt")
prevalences <- read_csv("data/environments-aitor.csv")
alpha_divs <- read_tsv("data/alpha_diversities.tsv")


engraftmentData <- left_join(engraftmentData, metaFmts, by = 'fmt_id')
engraftmentData <- left_join(engraftmentData, metaSamples, by = c('post_sample' = 'sample_name', 'dataset')) %>% rename('post_profile_id' = 'profile_id')
engraftmentData <- left_join(engraftmentData, metaSamples, by = c('single_pre_sample' = 'sample_name', 'dataset')) %>% rename('pre_profile_id' = 'profile_id')
engraftmentData <- left_join(engraftmentData, metaSamples, by = c('single_donor_sample' = 'sample_name', 'dataset')) %>% rename('donor_profile_id' = 'profile_id')

# Keep only principle time point
engraftmentData <- engraftmentData %>% filter(is_principal_timepoint)

# Clean up columns
engraftmentData <- engraftmentData %>% select(fmt_id,
                                              clade,
                                              dataset, 
                                              time_point,
                                              detected_mp_pre,
                                              detected_mp_post,
                                              detected_mp_donor,
                                              display_dataset, 
                                              pre_profile_id, 
                                              post_profile_id, 
                                              donor_profile_id,
                                              single_pre_sample,
                                              single_donor_sample) 


engraftmentData$dummyData <- F
dummyDataDonor <- list()
dummyDataRecipient <- list()
set.seed(1)

#######################################################################
# Simulate FMT triads with fake donors/recipients (For Figure 4D/E/F) #
#######################################################################

mmm <- read_tsv("data/metadata_fmts.tsv")
for(ds in unique(engraftmentData$dataset)){
  # find fmt-IDs with only unique donor individuals
  ii <-  engraftmentData %>% 
    filter(dataset == ds) %>% 
    filter(dummyData == FALSE) %>% 
    left_join(mmm, by = 'fmt_id') %>% 
    select(fmt_id, donor_name) %>% 
    distinct(donor_name, .keep_all = T) %>% 
    pull(fmt_id)
  for (fmt_id in ii){
    for (other_fmt_id in ii[ii != fmt_id]){
      fmt_ids <- c(fmt_id, other_fmt_id)
      stopifnot(length(unique(fmt_ids)) == 2)
      tmp <- engraftmentData %>% 
        filter(dummyData == FALSE) %>% 
        filter(fmt_id == fmt_ids[1]) %>% 
        select(-detected_mp_donor, -donor_profile_id)
      tmp2 <- engraftmentData %>% 
        filter(dummyData == FALSE) %>% 
        filter(fmt_id == fmt_ids[2]) %>% 
        select(clade, detected_mp_donor, donor_profile_id)
      test <- full_join(tmp, tmp2, by = 'clade')
      test$dummyData <- str_c("dummyDonor___", fmt_ids[1], "___", fmt_ids[2])
      test <- test %>% relocate(colnames(engraftmentData))
      # ... MARK THIS AS DUMMY DATA but otherwise make sure to normalize everything as usual.
      dummyDataDonor[[length(dummyDataDonor) + 1]] <- test
      
      ###
      # ... Now do the same as above but for recipients :)
      tmp <- engraftmentData %>% filter(fmt_id == fmt_ids[1]) %>% select(-detected_mp_pre, -pre_profile_id)
      tmp2 <- engraftmentData %>% filter(fmt_id == fmt_ids[2]) %>% select(clade, detected_mp_pre, pre_profile_id)
      test <- full_join(tmp, tmp2, by = 'clade')
      # ... MARK THIS AS DUMMY DATA but otherwise make sure to normalize everything as usual.
      test$dummyData <- str_c("dummyRecipient___",fmt_ids[1], "___", fmt_ids[2])
      # Bring columns into same order
      test <- test %>% relocate(colnames(engraftmentData))
      dummyDataRecipient[[length(dummyDataRecipient) + 1]] <- test 
    }
  }
}

dummyDataDonor <- do.call('rbind', dummyDataDonor)
dummyDataRecipient <- do.call('rbind', dummyDataRecipient)

engraftmentData <- rbind(engraftmentData, dummyDataDonor, dummyDataRecipient)

dummyDataBoth <- list()
set.seed(1)
for(ds in unique(engraftmentData$dataset)){
  ii <-  engraftmentData %>% 
    filter(dataset == ds) %>% 
    filter(dummyData == FALSE) %>% 
    left_join(mmm, by = 'fmt_id') %>% 
    select(fmt_id, donor_name) %>% 
    distinct(donor_name, .keep_all = T) %>% 
    pull(fmt_id)
  for (fmt_id in ii){
    for (other_fmt_id in ii[ii != fmt_id]){
      fmt_ids <- c(fmt_id, other_fmt_id)
      stopifnot(length(unique(fmt_ids)) == 2)
      tmp <- engraftmentData %>% filter(dummyData == FALSE) %>% filter(fmt_id == fmt_ids[1]) %>% select(-detected_mp_donor, -donor_profile_id, -detected_mp_pre, -pre_profile_id)
      tmp2 <- engraftmentData %>% filter(dummyData == FALSE) %>% filter(fmt_id == fmt_ids[2]) %>% select(clade, detected_mp_donor, donor_profile_id, detected_mp_pre, pre_profile_id)
      test <- full_join(tmp, tmp2, by = 'clade')
      test$dummyData <- str_c("dummyBoth___", fmt_ids[1], "___", fmt_ids[2])
      test <- test %>% relocate(colnames(engraftmentData))
      dummyDataBoth[[length(dummyDataBoth) + 1]] <- test
    }
  }
}


dummyDataBoth <- do.call('rbind', dummyDataBoth)
engraftmentData <- rbind(engraftmentData, dummyDataBoth)


mpa2ProfilesLong <- pivot_longer(mpa2Profiles, cols = c(-dataset, -profile_id), names_to = "clade")
mpa2ProfilesLong <- mpa2ProfilesLong %>% filter(str_detect(string = clade, pattern = "t__"))
ttt <- str_split(mpa2ProfilesLong$clade, "[|]", simplify = T)
ttt <- ttt %>% as.data.frame() %>% rename(phylum = V2,
                                          class = V3,
                                          order = V4,
                                          family = V5,
                                          genus = V6,
                                          species = V8) %>% select(-V7)
mpa2ProfilesLong <- cbind(mpa2ProfilesLong %>% select(-clade), ttt) %>% as_tibble()

# Get BC and Jaccard Distances between donors and recipients.
mpa2ProfilesForDist <- mpa2Profiles %>% select(-dataset) %>% as.data.frame() %>% distinct(profile_id, .keep_all = T)
rownames(mpa2ProfilesForDist) <- mpa2ProfilesForDist$profile_id
mpa2ProfilesForDist$profile_id <- NULL

mpa2ProfilesForDist <- mpa2ProfilesForDist[, apply(mpa2ProfilesForDist, 2, mean) != 0]
mpa2ProfilesForDist <- mpa2ProfilesForDist[apply(mpa2ProfilesForDist, 1, mean) != 0, ]

bray <- mpa2ProfilesForDist %>% vegdist(., method ='bray')
jaccard <- mpa2ProfilesForDist %>% vegdist(., method ='jaccard')

bray <- as.data.frame(as.matrix(bray))
bray$sampleID1 <- rownames(bray)
brayLong <- pivot_longer(bray, cols = -sampleID1)
colnames(brayLong) <- c("sampleID1", 'sampleID2', 'brayDistance')

jaccard <- as.data.frame(as.matrix(jaccard))
jaccard$sampleID1 <- rownames(jaccard)
jaccardLong <- pivot_longer(jaccard, cols = -sampleID1)
colnames(jaccardLong) <- c("sampleID1", 'sampleID2', 'jaccardDistance')

# Merge strain with metaphlan data
engraftmentData <- left_join(engraftmentData, mpa2ProfilesLong %>% select(profile_id, species, dataset, value), by = c("pre_profile_id" = "profile_id", c('clade' = 'species'), 'dataset')) %>% rename(preFMTAbundance = value)
engraftmentData <- left_join(engraftmentData, mpa2ProfilesLong %>% select(profile_id, species, dataset, value), by = c("post_profile_id" = "profile_id", c('clade' = 'species'), 'dataset')) %>% rename(postFMTAbundance = value)
engraftmentData <- left_join(engraftmentData, mpa2ProfilesLong , by = c("donor_profile_id" = "profile_id", c('clade' = 'species'), 'dataset')) %>% rename(donorAbundance = value)
 
engraftmentData <- left_join(engraftmentData, 
                             brayLong, 
                             by = c("pre_profile_id" = 'sampleID1', 
                                    'donor_profile_id' = 'sampleID2'))
engraftmentData <- left_join(engraftmentData, 
                             jaccardLong, 
                             by = c("pre_profile_id" = 'sampleID1', 
                                    'donor_profile_id' = 'sampleID2'))

# Restrict data to those data points that are seen in either recipient pre-FMT or donor
# Remaining data points ('double zeroes') are too easy to predict (are practically always 0)
engraftmentData <- engraftmentData %>% 
    filter(!is.na(detected_mp_pre)) %>% 
    filter(!is.na(detected_mp_post)) %>% 
    filter(!is.na(detected_mp_donor))
engraftmentData <- engraftmentData %>% 
    mutate(doubleZeroes = !(preFMTAbundance > 0 | donorAbundance > 0))

# Keep only columns of interest.
engraftmentData <- engraftmentData %>% select(dummyData, 
                                              dataset, 
                                              display_dataset, 
                                              single_pre_sample, 
                                              single_donor_sample, 
                                              pre_profile_id, 
                                              donor_profile_id, 
                                              clade,  
                                              time_point, 
                                              phylum, 
                                              class, 
                                              order, 
                                              family, 
                                              genus, 
                                              fmt_id, 
                                              preFMTAbundance, 
                                              postFMTAbundance, 
                                              donorAbundance, 
                                              brayDistance, 
                                              jaccardDistance, 
                                              doubleZeroes) %>% rename(species = clade)
# Make sure columns have the appropriate data type
engraftmentData <- engraftmentData %>% mutate(species = as.factor(species), 
                                              dataset = as.factor(dataset), 
                                              order = as.factor(order), 
                                              family = as.factor(family), 
                                              genus = as.factor(genus))
# get donorTopreFMTAbundanceRatio
engraftmentData <- engraftmentData %>% mutate(donorTopreFMTAbundanceRatio = donorAbundance/preFMTAbundance)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
colnames(clinicalMeta) <- map_chr(colnames(clinicalMeta), function(x) gsub(pattern = "[-%° \\(\\)\\.]", replacement = "_", x = x))
colnames(traitar) <- map_chr(colnames(traitar), function(x) gsub(pattern = "[-%° \\(\\)\\.]", replacement = "_", x = x))
traitar <- traitar %>% rename(species = SGB)
traitar <- left_join(data.frame(species = unique(engraftmentData$species)),
                     traitar)

fa <- c("disease_parsed", 
        "infectious_disease",
        "Infused_material_parsed")
clinicalMeta <- clinicalMeta %>% mutate(across(all_of(fa), as.factor))
traitar <- traitar %>% mutate(across(-species, as.factor))
# Not all SGBs have traitar traits. Impute them.
traitar[, -which(colnames(traitar) == "species")] <- missForest(xmis = traitar[, -which(colnames(traitar) == "species")],  verbose = T)$ximp



engraftmentData <- left_join(engraftmentData, clinicalMeta, by = c('dataset' = 'study'))
# mdPhenotypes is mostly na, so don't use it for now
#engraftmentData <- left_join(engraftmentData, mdPhenotypes, by = 'species')
engraftmentData <- left_join(engraftmentData, traitar, by = 'species')

prevalences$`SGB ID` <- map_chr(prevalences$`SGB ID`, function(x) str_c("t", str_split(string = x, pattern = "__")[[1]][1], sep = "__"))
prevalences$Prevalence <- map_dbl(prevalences$Prevalence, function(x) {
  tmp <-  str_replace(x, "[%]", "")
  tmp <-  str_replace(tmp, "[,]", ".")
  return(as.numeric(tmp))
})

prevalences <- pivot_wider(prevalences, id_cols = `SGB ID`, names_from = Environment, values_from = Prevalence) %>% select(`SGB ID`, Stool_W)

# Use metaFmts and metaSamples to link profiles to sample_names. 
# Then Generate alpha-div as features and add alpha-div of preFMT and of Donor to feature set. 
alpha_divs <- left_join(alpha_divs, metaSamples, 
                        by = c("dataset", "profile_id" = "profile_id")) %>% 
                                select(dataset, 
                                       sample_name, 
                                       alpha_diversity)
engraftmentData <- left_join(engraftmentData, alpha_divs, by = c("single_pre_sample" = "sample_name", "dataset")) %>% 
  rename(alpha_div_pre = alpha_diversity) %>% 
  left_join(alpha_divs, by = c("single_donor_sample" = "sample_name", "dataset")) %>% 
  rename(alpha_div_donor = alpha_diversity) 
engraftmentData <- engraftmentData %>% select(-single_pre_sample, -single_donor_sample)

prevalences$Stool_W[is.na(prevalences$Stool_W)] <- 0
engraftmentData <- left_join(engraftmentData, prevalences, by = c('species' = "SGB ID"))

# Very few are missing from prevalences. From what I can see all extremely rare species.
engraftmentData$Stool_W[is.na(engraftmentData$Stool_W)] <- 0

# Mutate all column types to factors (if not numerical, basically)
engraftmentData$dataset <- as.factor(engraftmentData$dataset)
colClasses <- lapply(engraftmentData, class)
names(colClasses) <- NULL
colClasses <- unlist(colClasses)
columnNames <- colnames(engraftmentData)[colClasses %in% c("character", "logical")]
engraftmentData <- engraftmentData %>% mutate(across(all_of(columnNames), as.factor))
engraftmentData$dataset <- as.factor(engraftmentData$dataset)

colClasses <- lapply(engraftmentData, class)
names(colClasses) <- NULL
colClasses <- unlist(colClasses)

features <- unique(c("phylum",
                     "class",
                     "order", 
                     "family", 
                     "genus", 
                     'species', 
                     'time_point', 
                     "donorAbundance", 
                     'preFMTAbundance',
                     'postFMTAbundance',
                     'donorTopreFMTAbundanceRatio',
                     "alpha_div_pre",
                     "alpha_div_donor",
                     "brayDistance",
                     'jaccardDistance',
                     "Stool_W", # Prevalence in westernized gut
                     #clinicalMeta %>% select(-study) %>% colnames(),
                     "infectious_disease"
                      #mdPhenotypes %>% select(-species) %>% colnames()
              ))
####################################################
# End of data reading, merging, preprocessing etc. #
####################################################

engraftmentData <- engraftmentData %>%
    mutate(doubleZeroes = as.logical(as.character(doubleZeroes)))
# Crucial
save <- engraftmentData
         
# For debugging                       
save.image("tmp/R_ML_image_3.rsave")

engr <- read_tsv("data/engraftment_metrics_species.tsv")
dat <- engr %>% 
select(clade, received_total_strains_support, received_total_strains_rate) %>%
filter(received_total_strains_support > 15) %>%
arrange(desc(received_total_strains_rate)) %>%
left_join(read_tsv("data/SGB_taxonomys.tab", col_names = F) %>%
         rename(clade = X1,
               species = X2) %>%
         mutate(species = map_chr(species, function(x) str_split(x, "[|]")[[1]][7])), by = 'clade') %>%
left_join(engraftmentData %>% 
         select(species, donorAbundance) %>%
         group_by(species) %>%
         summarize(meanDonorAbundance = mean(donorAbundance),
                  donorPrevalence = mean(donorAbundance > 0)) %>%
         rename(clade = species), by = 'clade')
c1 <- round(cor(dat$received_total_strains_rate,
         dat$meanDonorAbundance, method = 'spearman'), 3)
c2 <- round(cor(dat$received_total_strains_rate,
         dat$donorPrevalence, method = 'spearman'), 3)
c1.t <- round(cor.test(dat$received_total_strains_rate,
         dat$meanDonorAbundance, method = 'spearman')$p.value, 3)
c2.t <- round(cor.test(dat$received_total_strains_rate,
         dat$donorPrevalence, method = 'spearman')$p.value, 3)                                  
p1 <- ggplot(dat, aes(x = received_total_strains_rate, y = meanDonorAbundance)) +
                                  geom_point() +
                                  theme_embl() +
                                  xlab("strain engraftment rate") +
                                  ylab("mean abundance in donors") +
                                  annotate(geom = 'text',
                                           x = 0.2, y = 4,
                                            label = str_c("Spearmon's rho: ", c1, "\np-value: ", "<0.001"))
p2 <- ggplot(dat, aes(x = received_total_strains_rate, y = donorPrevalence)) +
                                  geom_point() +
                                  theme_embl() +
                                  xlab("strain engraftment rate") +
                                  ylab("prevalence in donors") +
                                  annotate(geom = 'text',
                                           x = 0.2, y = 1.2,
                                            label = str_c("Spearmon's rho: ", c2, "\np-value: ", c2.t))
                                  
options(repr.plot.width=11, repr.plot.height=5)
library(patchwork)
# This Figure isn't in the final manuscript
ggsave(plot = p1 + p2, filename = "figures/strain_engraftment_rate_vs_donor_abundance_and_prevalence.pdf", width = 8, height = 3.5)


#############
# CV AUROCS #
#############

results <- list()
set.seed(1337)
numIterations <- 5
folds <- 5

for (ds in save %>% pull(display_dataset) %>% unique()){
  engraftmentData <- save %>%
    filter(dummyData == "FALSE") %>%  
    filter(!doubleZeroes) %>% 
    filter(display_dataset == ds)
  
  rfModel <- lrn("classif.ranger", importance = "permutation")
  rfModel$predict_type <- 'prob'
  engraftmentClassifier <- TaskClassif$new(id = "engraftmentClassifier", backend = engraftmentData %>% 
                                             select(all_of(c(features, 
                                                             "postFMTAbundance", 
                                                             'display_dataset'))) %>% 
                                             mutate(postFMTAbundance = ifelse(postFMTAbundance > 0, 
                                                                              "1", 
                                                                              "0")) %>% 
                                             mutate(postFMTAbundance = as.factor(postFMTAbundance)), 
                                           target = "postFMTAbundance", positive = "1")
  cv <- rsmp('repeated_cv', folds = folds, repeats = numIterations)
  cv$instantiate(engraftmentClassifier)
  
  for (i in 1:(folds * numIterations)){
    train_set <- cv$train_set(i)
    test_set <- cv$test_set(i)
    rfModel$train(task = engraftmentClassifier, row_ids = train_set)
    pred <- rfModel$predict(engraftmentClassifier, row_ids = test_set)
    pred <- cbind(pred$row_ids, pred$truth, pred$prob[, 1], pred$response)
    colnames(pred) <- c('row_id', 'truth', 'prob', 'response')
    pred <- as.data.frame(pred)
    pred$display_dataset <- ds
    pred$repeats <- cv$repeats(i)
    r <- cbind(pred, engraftmentData[test_set, ] %>% select(fmt_id)) %>%
      group_by(fmt_id) %>% 
      nest()
    results[[length(results) + 1]] <- r
  }
}

# For debugging
save.image("tmp/R_ML_image_4.rsave")

# Get confusion matrices/accuracy for classifier
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
CMs <- map(do.call('rbind', results) %>% 
            unnest() %>% 
            group_by(fmt_id, row_id, display_dataset) %>%
            summarize(truth = mean(truth),
                      prob = mean(prob),
                      response = Mode(response)) %>%
           group_by(display_dataset) %>%
           nest() %>%
           pull(data), function(x){
    # At some point in the scripts above a factor with levels c(0, 1)
    # Gets cast to numeric which turns it into c(1,2). So turn back for clarity.
    #x$truth <- ifelse(x$truth == 2, 1, 0)
    #x$response <- ifelse(x$response == 2, 1, 0)
    #filter(preFMTAbundance != 0 | donorAbundance != 0)
    tt <- confusionMatrix(factor(x %>% pull(truth), levels = c(1,2)),
               factor(x %>% pull(response), levels = c(1,2)))
    return(list(tt$table, tt$overall[['Accuracy']]))
    
})


CMs <- tibble(confTab = map(CMs, function(x) x[[1]] %>% as.data.frame()),
              accuracy = map(CMs, function(x) x[[2]]))
CMs$ds <- (save %>% pull(display_dataset) %>% unique())                        
CMs <- CMs %>%
mutate(plot = pmap(list(confTab, accuracy, ds), function(a, b, d) {
  ggplot(data =  a %>%
         rename(Truth = Reference), mapping = aes(x = Truth, y = Prediction)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = round(Freq/sum(a$Freq), 3)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  theme_embl() + 
    theme(legend.position = "none") +
    ggtitle(str_c("Dataset: ", d, "\nAccuracy: ", round(b, 3))) %>%
    return()
}))
options(repr.plot.width=14, repr.plot.height=17)
ggsave(plot = wrap_plots(CMs$plot, ncol = 5), filename = "figures/ED_Figure_6.pdf", width = 14, height = 17)
wrap_plots(CMs$plot, ncol = 5)


#b <- do.call('rbind', resultsForCVManyFeatures) %>%
#unnest() %>% 
#group_by(fmt_id, 
#         row_id, 
#         display_dataset) %>% 
#summarize(prob = mean(prob), 
#          truth = unique(truth)) %>%
#group_by(fmt_id, display_dataset) %>%
#nest() %>% 
#filter(map_lgl(data, function(x) return(length(unique(x$truth)) != 1))) %>%
#      mutate(r = map(data, function(x) return(suppressMessages(roc(x$truth, x$prob))))) %>%
#      mutate(r = map(r, function(x) {
#        tmp <- (data.frame(x$sensitivities, 
#                           x$specificities, 
#                           ds, 
#                           length(test_set), 
#                           round(x$auc, 3), 
#                           iteration = i, 
#                           pon = 1:length(x$sensitivities)))
#        colnames(tmp) <- c("Sensitivity", "Specificity", "testDataset", "testInstances", 'auc', "iteration", "pos")
#        return(tmp)
#      })) %>%
#select(fmt_id, display_dataset, r) %>% unnest() %>%
#                     group_by(fmt_id) %>%
#                     summarize(auc = mean(auc))


#a <- do.call('rbind', resultsForCVLittleFeatures) %>%
#unnest() %>% 
#group_by(fmt_id, 
#         row_id, 
#         display_dataset) %>% 
#summarize(prob = mean(prob), 
#          truth = unique(truth)) %>%
#group_by(fmt_id, display_dataset) %>%
#nest() %>% 
#filter(map_lgl(data, function(x) return(length(unique(x$truth)) != 1))) %>%
#      mutate(r = map(data, function(x) return(suppressMessages(roc(x$truth, x$prob))))) %>%
#      mutate(r = map(r, function(x) {
#        tmp <- (data.frame(x$sensitivities, 
#                           x$specificities, 
#                           ds, 
#                           length(test_set), 
#                           round(x$auc, 3), 
#                           iteration = i, 
#                           pon = 1:length(x$sensitivities)))
#        colnames(tmp) <- c("Sensitivity", "Specificity", "testDataset", "testInstances", 'auc', "iteration", "pos")
#        return(tmp)
#      })) %>%
#select(fmt_id, display_dataset, r) %>% unnest() %>%
#                     group_by(fmt_id) %>%
#                     summarize(auc = mean(auc))

#full_join(a,  b, by = 'fmt_id', suffix = c(".few", ".many")) %>%
#ggplot(aes(x = auc.few, y = auc.many)) +
#geom_point() +
#geom_abline(slope = 1, intercept = 0) +
#theme_embl() +
#xlab("AUC (big model)") +
#ylab("AUC (small model)")
#ggsave(filename = "figures/CV_AUCs_small_big_model.pdf", width = 4, height = 4)


resultsForCV <- results
results <- do.call('rbind', results)
results <- results %>% 
unnest() %>% 
group_by(fmt_id, 
         row_id, 
         display_dataset) %>% 
summarize(prob = mean(prob), 
          truth = unique(truth)) %>%
group_by(fmt_id, display_dataset) %>%
nest() %>% 
filter(map_lgl(data, function(x) return(length(unique(x$truth)) != 1))) %>%
      mutate(r = map(data, function(x) return(suppressMessages(roc(x$truth, x$prob))))) %>%
      mutate(r = map(r, function(x) {
        tmp <- (data.frame(x$sensitivities, 
                           x$specificities, 
                           ds, 
                           length(test_set), 
                           round(x$auc, 3), 
                           iteration = i, 
                           pon = 1:length(x$sensitivities)))
        colnames(tmp) <- c("Sensitivity", "Specificity", "testDataset", "testInstances", 'auc', "iteration", "pos")
        return(tmp)
      })) %>%
select(fmt_id, display_dataset, r) %>% unnest()
cvaurocs_curve_per_FMT_ID <- results

cvaurocs <- results %>% mutate(display_dataset = testDataset) %>%
  ggplot(aes(x = 1-Sensitivity, Specificity, group = fmt_id)) + 
  geom_line(alpha = 0.1) + 
  theme_embl() +
  xlab("FPR") + 
  ylab("TPR")

ggsave(plot = cvaurocs, filename = "figures/Figure_4A_Panel_right.pdf", width = 2.5, height = 2.5)

results <- list()
resultsPred <- list()
abundancePredictions <- list()
models <- list()
set.seed(1337)
numIterations <- 1
aucsDummyDonor <- list()
aucsDummyRecipient <- list()
aucsDummyBoth <-
richnessDummyDonor <- list()
richnessPredictedUnalteredInstance <- list()
abundancePredictedUnalteredInstance <- list()
CumprevProteoDummyDonor <- list()
cumAbundanceProteoDummyDonor <- list()
CumprevOralDummyDonor <- list()
cumAbundanceOralDummyDonor <- list()
CumprevTopPredict1DummyDonor <- list()
cumAbundanceTopPredict1DummyDonor <- list()

cumAbundanceFirmicutesDummyDonor <- list()
cumAbundanceBacteroidetesDummyDonor <- list()
cumPrevFirmicutesDummyDonor <- list()
cumPrevBacteroidetesDummyDonor <- list() 

feature_importances <- list()

oralAbPrev <- read_tsv("data/oral_abundances_and_prevalences.tsv")
oralAbPrev <- read_tsv("data/SGB_taxonomys.tab", col_names = F) %>%
mutate(speciesJoin = map_chr(X2, function(x) str_split(x, "[|]")[[1]][7]),
   species = map_chr(X2, function(x) str_split(x, "[|]")[[1]][8])) %>% 
                      select(species, speciesJoin) %>%
                     right_join(oralAbPrev, by = c('speciesJoin' = 'species')) %>%
                     # 3 out of >2k species are lost...
                     filter(!is.na(species)) %>%
                     filter(prevalence > 0.2) %>%
                     pull(species)
                     
topPredict1RankSpecies <- read_tsv("data/Ranking_PREDICT1_mpa4beta.tsv") %>% 
                  select(`...1`, `AVG rank`) %>% 
                  arrange(`AVG rank`) %>% 
                  head(80) %>% 
                  rename(species = `...1`) %>%
                  mutate(species = map_chr(species, function(x) str_split(x, "[|]")[[1]][2])) %>%
                  pull(species)

numIterations <- 1 # Should always be 1...
for (iteration in 1:numIterations){
  engraftmentData <- save %>%
    filter(dummyData == "FALSE") %>%
    mutate(donorTopreFMTAbundanceRatio = ifelse(is.na(donorTopreFMTAbundanceRatio), 
                                                0, 
                                                donorTopreFMTAbundanceRatio))
    #filter(!doubleZeroes)
   for (ds in (engraftmentData %>% pull(display_dataset) %>% unique())){
    rfModel_reg <- lrn("regr.ranger", importance = 'permutation')
    rfModel_class <- lrn("classif.ranger", importance = 'permutation')
    #rfModel_class <- lrn("classif.ranger")
    rfModel_class$predict_type <- 'prob'
    engraftmentR <- TaskRegr$new(id = "engraftmentR", backend = engraftmentData %>% 
                                               select(all_of(c(features, 
                                                               "postFMTAbundance", 
                                                               'display_dataset'))),
                                             target = "postFMTAbundance")
    engraftmentC <- TaskClassif$new(id = "engraftmentC", backend = engraftmentData %>% 
                                           select(all_of(c(features, 
                                                           "postFMTAbundance", 
                                                           'display_dataset'))) %>%
                                           mutate(postFMTAbundance = ifelse(postFMTAbundance > 0, 
                                                                            "1", 
                                                                            "0")) %>% 
                                           mutate(postFMTAbundance = as.factor(postFMTAbundance)),
                                         target = "postFMTAbundance")
      
        # Set training and test set
    train_set <- engraftmentR$data() %>% mutate(rowId = 1:dim(engraftmentR$data())[1]) %>% 
      filter(display_dataset != ds) %>% 
      pull(rowId) 
    test_set <- engraftmentR$data() %>% mutate(rowId = 1:dim(engraftmentR$data())[1]) %>% 
      filter(display_dataset == ds) %>% 
      pull(rowId)   
       
      a <- save %>%
      filter(str_detect(dummyData, "dummyDonor")) %>%
      mutate(donorTopreFMTAbundanceRatio = ifelse(is.na(donorTopreFMTAbundanceRatio), 
                                            0, 
                                            donorTopreFMTAbundanceRatio)) %>%
      # filter(!doubleZeroes) %>%
      filter(display_dataset == ds) %>%
      select(all_of(c(features, "dummyData", "postFMTAbundance", 'display_dataset'))) %>%
      mutate(postFMTAbundanceN = postFMTAbundance) %>%
      mutate(postFMTAbundance = ifelse(postFMTAbundance > 0, "1", "0")) %>%
      mutate(postFMTAbundance = as.factor(postFMTAbundance)) %>% 
      group_by(dummyData) %>%
      nest()
       
    engraftmentClassifierDummyDonor <- purrr::map(a$data, function(x) {
      TaskClassif$new(id = "engraftmentClassifierMockDonor", 
                      backend = x, 
                      target = "postFMTAbundance", 
                      positive = "1")
    })
    names(engraftmentClassifierDummyDonor) <- a$dummyData
    engraftmentRegressorDummyDonor <- purrr::map(a$data, function(x) {
      TaskRegr$new(id = "engraftmentRegressorMockDonor", 
                      backend = x, 
                      target = "postFMTAbundanceN")
    })
    names(engraftmentRegressorDummyDonor) <- a$dummyData
       
       
      a <- save %>%
      filter(str_detect(dummyData, "dummyRecipient")) %>%
      mutate(donorTopreFMTAbundanceRatio = ifelse(is.na(donorTopreFMTAbundanceRatio), 
                                            0, 
                                            donorTopreFMTAbundanceRatio)) %>%
      # filter(!doubleZeroes) %>%
      filter(display_dataset == ds) %>%
      select(all_of(c(features, "dummyData", "postFMTAbundance", 'display_dataset'))) %>%
      mutate(postFMTAbundanceN = postFMTAbundance) %>%
      mutate(postFMTAbundance = ifelse(postFMTAbundance > 0, "1", "0")) %>%
      mutate(postFMTAbundance = as.factor(postFMTAbundance)) %>% 
      group_by(dummyData) %>%
      nest()
    engraftmentClassifierDummyRecipient <- purrr::map(a$data, function(x) {
      TaskClassif$new(id = "engraftmentClassifierMockRecipient", 
                      backend = x, 
                      target = "postFMTAbundance", 
                      positive = "1")
    })
    names(engraftmentClassifierDummyRecipient) <- a$dummyData
    engraftmentRegressorDummyRecipient <- purrr::map(a$data, function(x) {
      TaskRegr$new(id = "engraftmentClassifierMockRecipient", 
                   backend = x, 
                   target = "postFMTAbundanceN")
    })
    names(engraftmentRegressorDummyRecipient) <- a$dummyData
       
       
      a <- save %>%
      filter(str_detect(dummyData, "dummyBoth")) %>%
      mutate(donorTopreFMTAbundanceRatio = ifelse(is.na(donorTopreFMTAbundanceRatio), 
                                            0, 
                                            donorTopreFMTAbundanceRatio)) %>%
      # filter(!doubleZeroes) %>%
      filter(display_dataset == ds) %>%
      select(all_of(c(features, "dummyData", "postFMTAbundance", 'display_dataset'))) %>%
      mutate(postFMTAbundanceN = postFMTAbundance) %>%
      mutate(postFMTAbundance = ifelse(postFMTAbundance > 0, "1", "0")) %>%
      mutate(postFMTAbundance = as.factor(postFMTAbundance)) %>% 
      group_by(dummyData) %>%
      nest()
    engraftmentClassifierDummyBoth <- purrr::map(a$data, function(x) {
      TaskClassif$new(id = "engraftmentClassifierMockBoth", 
                      backend = x, 
                      target = "postFMTAbundance", 
                      positive = "1")
    })
    names(engraftmentClassifierDummyBoth) <- a$dummyData
    engraftmentRegressorDummyBoth <- purrr::map(a$data, function(x) {
      TaskRegr$new(id = "engraftmentClassifierMockBoth", 
                   backend = x, 
                   target = "postFMTAbundanceN")
    })
    names(engraftmentRegressorDummyBoth) <- a$dummyData       
       
       
    if (length(test_set) < 2){
      next
    }
       
rfModel_reg$train(task = engraftmentR, row_ids = train_set)
rfModel_class$train(task = engraftmentC, row_ids = train_set)
predR <- rfModel_reg$predict(engraftmentR, row_ids = test_set)
pred <- rfModel_class$predict(engraftmentC, row_ids = test_set)
       
predBefore <- pred
       
models[[length(models) + 1 ]] <- rfModel_class
       
resultsPred[[length(resultsPred) + 1]] <- cbind(pred %>% 
as.data.table() %>%
as.data.frame(), engraftmentData[test_set, ] %>% select(preFMTAbundance, donorAbundance, fmt_id))

predR <- predR %>% 
   as.data.table()  %>%
   as.data.frame() %>% 
   left_join(pred %>%
        as.data.table() %>%
        as.data.frame() %>%
        select(row_ids, response) %>%
         rename(response_C = response), by = 'row_ids') %>%
   mutate(response = ifelse(response_C == 0, 0, response)) %>%
 select(-response_C)
predR$ds <- ds
abundancePredictions[[length(abundancePredictions) + 1]] <- cbind(predR, 
                                                                  engraftmentData[test_set, ] %>% 
                                                                  select(preFMTAbundance, donorAbundance, fmt_id))
    
    stopifnot(all(levels(pred$truth) %in% c(1,0)))
    richnessPredictedUnalteredInstance[[length(richnessPredictedUnalteredInstance) + 1]] <- engraftmentData %>%
      mutate(row_ids = 1:dim(.)[1]) %>% 
      filter(row_ids %in% pred$row_ids) %>%
      select(row_ids, fmt_id, phylum, species) %>%
      #mutate(row_ids = as.character(row_ids))
      left_join(cbind(pred$row_ids, pred %>% as.data.table() %>% as.data.frame() %>% pull(prob.1), as.vector(pred$truth)) %>% 
                  as.data.frame() %>% 
                  as_tibble() %>% 
                  rename(row_ids = V1, propability = V2, truth = V3) %>%
                  mutate(row_ids = as.numeric(row_ids)) %>%
                  mutate(propability = as.numeric(propability)), by = 'row_ids') %>%
      group_by(fmt_id) %>% 
      nest() %>% 
      mutate(richness = map_dbl(data, function(x){
        #return(sum(x$propability > optimalThresh)) # This is the real richness
        return(sum(x$propability)) # This is the sum of the prediction probabilities
      })) %>%
      mutate(cumPrevalenceProteobacteria = map_dbl(data, function(x) {
        return(x %>% filter(phylum == "p__Proteobacteria") %>% pull(propability) %>% sum())
      })) %>%
      mutate(cumPrevalencePredict1 = map_dbl(data, function(x) {
        return(x %>% filter(species %in% topPredict1RankSpecies) %>% pull(propability) %>% sum())
      })) %>%
      mutate(cumPrevalenceOral = map_dbl(data, function(x) {
        return(x %>% filter(species %in% oralAbPrev) %>% pull(propability) %>% sum())
      })) %>%
      mutate(cumPrevalenceFirmicutes = map_dbl(data, function(x) {
        return(x %>% filter(phylum == "p__Firmicutes") %>% pull(propability) %>% sum())
      })) %>%
      mutate(cumPrevalenceBacteroidetes = map_dbl(data, function(x) {
        return(x %>% filter(phylum == "p__Bacteroidetes") %>% pull(propability) %>% sum())
      })) %>%       
      select(fmt_id, 
             richness, 
             cumPrevalenceProteobacteria,
             #cumAbundanceProteobacteria,
             cumPrevalencePredict1,
             #cumAbundancePredict1,
             cumPrevalenceOral,
             cumPrevalenceFirmicutes,
             cumPrevalenceBacteroidetes)
    
    abundancePredictedUnalteredInstance[[length(abundancePredictedUnalteredInstance) + 1]] <- engraftmentData %>%
      mutate(row_ids = 1:dim(.)[1]) %>% 
      filter(row_ids %in% predR$row_ids) %>%
      select(row_ids, fmt_id, phylum, species) %>%
      #mutate(row_ids = as.character(row_ids))
      left_join(cbind(predR$row_ids, predR$response, as.vector(predR$truth)) %>% 
                  as.data.frame() %>% 
                  as_tibble() %>% 
                  rename(row_ids = V1, response = V2, truth = V3) %>%
                  mutate(row_ids = as.numeric(row_ids)) %>%
                  mutate(response = as.numeric(response)), by = 'row_ids') %>%
      group_by(fmt_id) %>% 
      nest() %>% 
      mutate(cumAbundanceProteobacteria = map_dbl(data, function(x) {
        return(x %>% filter(phylum == "p__Proteobacteria") %>% pull(response) %>% sum())
      })) %>%
      mutate(cumAbundanceFirmicutes = map_dbl(data, function(x) {
        return(x %>% filter(phylum == "p__Firmicutes") %>% pull(response) %>% sum())
      })) %>%       
      mutate(cumAbundanceBacteroidetes = map_dbl(data, function(x) {
        return(x %>% filter(phylum == "p__Bacteroidetes") %>% pull(response) %>% sum())
      })) %>%            
      mutate(cumAbundancePredict1 = map_dbl(data, function(x) {
        return(x %>% filter(species %in% topPredict1RankSpecies) %>% pull(response) %>% sum())
      })) %>%
      mutate(cumAbundanceOral = map_dbl(data, function(x) {
        return(x %>% filter(species %in% oralAbPrev) %>% pull(response) %>% sum())
      })) %>%
      select(fmt_id, 
             #richness, 
             #cumPrevalenceProteobacteria,
             cumAbundanceProteobacteria,
             cumAbundanceFirmicutes,
             cumAbundanceBacteroidetes,
             #cumPrevalencePredict1,
             cumAbundancePredict1,
             #cumPrevalenceOral)
             cumAbundanceOral)
    pred <- cbind(pred$row_ids, pred$truth, pred %>% as.data.table() %>% as.data.frame() %>% pull(prob.1))
    colnames(pred) <- c('row_id', 'truth', 'prob')
    pred <- as.data.frame(pred)
    predDummyDonor <- map(engraftmentClassifierDummyDonor, function(x){
      tmp <- rfModel_class$predict(x, row_ids = 1:x$backend$nrow)
      #tmp$set_threshold(optimalThresh)
      tmp <- tmp %>%
        as.data.table() %>%
        as.data.frame()
      tmp <- cbind(tmp, x$data() %>% select(phylum, class, order, family, genus, species))
      return(tmp)
    })
    predDummyRecipient <- map(engraftmentClassifierDummyRecipient, function(x){
      tmp <- rfModel_class$predict(x, row_ids = 1:x$backend$nrow)
      #tmp$set_threshold(optimalThresh)
      tmp <- tmp %>%
        as.data.table() %>%
        as.data.frame()
      tmp <- cbind(tmp, x$data() %>% select(phylum, class, order, family, genus, species))
      return(tmp)
    })
    predDummyBoth <- map(engraftmentClassifierDummyBoth, function(x){
      tmp <- rfModel_class$predict(x, row_ids = 1:x$backend$nrow)
      #tmp$set_threshold(optimalThresh)
      tmp <- tmp %>%
        as.data.table() %>%
        as.data.frame()
      tmp <- cbind(tmp, x$data() %>% select(phylum, class, order, family, genus, species))
      return(tmp)
    })
    r <- suppressMessages(roc(pred$truth, pred$prob))
       
    predAbundanceDummyDonor <- map2(engraftmentClassifierDummyDonor, 
                                    engraftmentRegressorDummyDonor, function(x, y){

      pred <- rfModel_class$predict(x, row_ids = 1:x$backend$nrow)
      #pred$set_threshold(optimalThresh)                                        
      predR <- rfModel_reg$predict(y, row_ids = 1:x$backend$nrow)
      predR <- predR %>%
                                        as.data.table() %>%
                                        as.data.frame() %>%
                                        left_join(pred %>%
         as.data.table() %>%
         as.data.frame() %>%
         select(row_ids, response) %>%
          rename(response_C = response), by = 'row_ids') %>%
    mutate(response = ifelse(response_C == 0, 0, response)) %>%
    select(-response_C)                        
      predR <- cbind(predR, x$data() %>% select(phylum, class, order, family, genus, species))
      return(predR)
    })
    predAbundanceDummyRecipient <- map2(engraftmentClassifierDummyRecipient, 
                                    engraftmentRegressorDummyRecipient, function(x, y){

      pred <- rfModel_class$predict(x, row_ids = 1:x$backend$nrow)
      #pred$set_threshold(optimalThresh)                                        
      predR <- rfModel_reg$predict(y, row_ids = 1:x$backend$nrow)
      predR <- predR %>%
                                        as.data.table() %>%
                                        as.data.frame() %>%
                                        left_join(pred %>%
         as.data.table() %>%
         as.data.frame() %>%
         select(row_ids, response) %>%
          rename(response_C = response), by = 'row_ids') %>%
    mutate(response = ifelse(response_C == 0, 0, response)) %>%
    select(-response_C)                        
      predR <- cbind(predR, x$data() %>% select(phylum, class, order, family, genus, species))
      return(predR)
    })
    
    #rDummyDonor <- suppressMessages(roc(predDummyDonor$truth, predDummyDonor$prob[, 1]))
    rDummyDonor <- map(predDummyDonor, function(x) suppressMessages(roc(x$truth, x$prob.1)))
    rDummyRecipient <- map(predDummyRecipient, function(x) suppressMessages(roc(x$truth, x$prob.1)))
    rDummyBoth <- map(predDummyBoth, function(x) suppressMessages(roc(x$truth, x$prob.1)))
    
    aucsDummyDonor[[length(aucsDummyDonor) + 1]] <- rDummyDonor
    names(aucsDummyDonor)[length(aucsDummyDonor)] <- str_c(ds, "__", iteration)
    aucsDummyRecipient[[length(aucsDummyRecipient) + 1]] <- rDummyRecipient
    names(aucsDummyRecipient)[length(aucsDummyRecipient)] <- str_c(ds, "__", iteration)
    aucsDummyBoth[[length(aucsDummyBoth) + 1]] <- rDummyBoth
    names(aucsDummyBoth)[length(aucsDummyBoth)] <- str_c(ds, "__", iteration)
    # CumprevOralDummyDonor <- list()
    # CumprevTopPredict1DummyDonor <- list()
    # /shares/CIBIO-Storage/CM/scratch/projects/nkarcher_FMT_meta/data/oral_abundances_and_prevalences.tsv                  
    # It could be some species in the file below are not in our files 
    # read_tsv("/shares/CIBIO-Storage/CM/scratch/projects/nkarcher_FMT_meta/data/Ranking_PREDICT1_mpa4beta.tsv") %>% 
                      # select(`...1`, `AVG rank`) %>% 
                      # arrange(`AVG rank`) %>% 
                      # head(80) %>% 
                      # rename(species = `...1`) %>%
                      # mutate(species = map_chr(species, function(x) str_split(x, "[|]")[[1]][2])) %>%
                      # pull(species)
                      
    #asdsaddsa
    #richnessDummyDonor[[length(richnessDummyDonor) + 1]] <- map(predDummyDonor, function(x) sum(x$confusion[1, ])) # This is the real richness 
    richnessDummyDonor[[length(richnessDummyDonor) + 1]] <- map(predDummyDonor, function(x) {
        sum(x$prob.1)
    }) # This is the sum of the prediction probabilities
    names(richnessDummyDonor[[length(richnessDummyDonor)]]) <- names(engraftmentClassifierDummyDonor)
    names(richnessDummyDonor)[length(richnessDummyDonor)] <- str_c(ds, "__", iteration)
    
    CumprevProteoDummyDonor[[length(CumprevProteoDummyDonor) + 1]] <- map(predDummyDonor, function(x) {
        sum(x %>% filter(phylum == "p__Proteobacteria") %>% pull(prob.1))
    })
    names(CumprevProteoDummyDonor[[length(CumprevProteoDummyDonor)]]) <- names(engraftmentClassifierDummyDonor)
    names(CumprevProteoDummyDonor)[length(CumprevProteoDummyDonor)] <- str_c(ds, "__", iteration)

                      
    cumPrevFirmicutesDummyDonor[[length(cumPrevFirmicutesDummyDonor) + 1]] <- map(predDummyDonor, function(x) {
        sum(x %>% filter(phylum == "p__Firmicutes") %>% pull(prob.1))
    })
    names(cumPrevFirmicutesDummyDonor[[length(cumPrevFirmicutesDummyDonor)]]) <- names(engraftmentClassifierDummyDonor)
    names(cumPrevFirmicutesDummyDonor)[length(cumPrevFirmicutesDummyDonor)] <- str_c(ds, "__", iteration)
                      
    cumPrevBacteroidetesDummyDonor[[length(cumPrevBacteroidetesDummyDonor) + 1]] <- map(predDummyDonor, function(x) {
        sum(x %>% filter(phylum == "p__Bacteroidetes") %>% pull(prob.1))
    })
    names(cumPrevBacteroidetesDummyDonor[[length(cumPrevBacteroidetesDummyDonor)]]) <- names(engraftmentClassifierDummyDonor)
    names(cumPrevBacteroidetesDummyDonor)[length(cumPrevBacteroidetesDummyDonor)] <- str_c(ds, "__", iteration)       
                      

    
    CumprevOralDummyDonor[[length(CumprevOralDummyDonor) + 1]] <- map(predDummyDonor, function(x) {
        x %>% filter(species %in% oralAbPrev) %>% pull(prob.1) %>% sum()
    })
    names(CumprevOralDummyDonor[[length(CumprevOralDummyDonor)]]) <- names(engraftmentClassifierDummyDonor)
    names(CumprevOralDummyDonor)[length(CumprevOralDummyDonor)] <- str_c(ds, "__", iteration)
    
    CumprevTopPredict1DummyDonor[[length(CumprevTopPredict1DummyDonor) + 1]] <- map(predDummyDonor, function(x) {
        x %>% filter(species %in% topPredict1RankSpecies) %>% pull(prob.1) %>% sum()
    })
    names(CumprevTopPredict1DummyDonor[[length(CumprevTopPredict1DummyDonor)]]) <- names(engraftmentClassifierDummyDonor)
    names(CumprevTopPredict1DummyDonor)[length(CumprevTopPredict1DummyDonor)] <- str_c(ds, "__", iteration)
    
    cumAbundanceProteoDummyDonor[[length(cumAbundanceProteoDummyDonor) + 1]] <- map(predAbundanceDummyDonor, function(x){
        x %>% filter(phylum == "p__Proteobacteria") %>% pull(response) %>% sum() %>% return()
    })
    names(cumAbundanceProteoDummyDonor[[length(cumAbundanceProteoDummyDonor)]]) <- names(engraftmentRegressorDummyDonor)
    names(cumAbundanceProteoDummyDonor)[length(cumAbundanceProteoDummyDonor)] <- str_c(ds, "__", iteration)
                                                
    cumAbundanceOralDummyDonor[[length(cumAbundanceOralDummyDonor) + 1]] <- map(predAbundanceDummyDonor, function(x){
        x %>% filter(species %in% oralAbPrev) %>% pull(response) %>% sum() %>% return()
    })
    names(cumAbundanceOralDummyDonor[[length(cumAbundanceOralDummyDonor)]]) <- names(engraftmentRegressorDummyDonor)
    names(cumAbundanceOralDummyDonor)[length(cumAbundanceOralDummyDonor)] <- str_c(ds, "__", iteration)
                      
    cumAbundanceTopPredict1DummyDonor[[length(cumAbundanceTopPredict1DummyDonor) + 1]] <- map(predAbundanceDummyDonor, function(x){
        x %>% filter(species %in% topPredict1RankSpecies) %>% pull(response) %>% sum() %>% return()
    })
    names(cumAbundanceTopPredict1DummyDonor[[length(cumAbundanceTopPredict1DummyDonor)]]) <- names(engraftmentRegressorDummyDonor)
    names(cumAbundanceTopPredict1DummyDonor)[length(cumAbundanceTopPredict1DummyDonor)] <- str_c(ds, "__", iteration)                  
                      
    cumAbundanceFirmicutesDummyDonor[[length(cumAbundanceFirmicutesDummyDonor) + 1]] <- map(predAbundanceDummyDonor, function(x){
        x %>% filter(phylum == "p__Firmicutes") %>% pull(response) %>% sum() %>% return()
    })
    names(cumAbundanceFirmicutesDummyDonor[[length(cumAbundanceFirmicutesDummyDonor)]]) <- names(engraftmentRegressorDummyDonor)
    names(cumAbundanceFirmicutesDummyDonor)[length(cumAbundanceFirmicutesDummyDonor)] <- str_c(ds, "__", iteration)
                      
    cumAbundanceBacteroidetesDummyDonor[[length(cumAbundanceBacteroidetesDummyDonor) + 1]] <- map(predAbundanceDummyDonor, function(x){
        x %>% filter(phylum == "p__Bacteroidetes") %>% pull(response) %>% sum() %>% return()
    })
    names(cumAbundanceBacteroidetesDummyDonor[[length(cumAbundanceBacteroidetesDummyDonor)]]) <- names(engraftmentRegressorDummyDonor)
    names(cumAbundanceBacteroidetesDummyDonor)[length(cumAbundanceBacteroidetesDummyDonor)] <- str_c(ds, "__", iteration)
         
                                                                                   
    r <- data.frame(r$sensitivities, r$specificities, ds, length(test_set), round(r$auc, 3))
    r$iteration <- iteration
    r$pon <- 1:dim(r)[1]
    # r$infectious <- bool
    colnames(r) <- c("Sensitivity", "Specificity", "testDataset", "testInstances", 'auc', "iteration", "pos")
    results[[length(results) + 1]] <- r
    feature_importances[[length(feature_importances) + 1]] <- list(rfModel_class$importance(), rfModel_reg$importance())
    names(feature_importances)[length(feature_importances)] <- str_c(ds,  "__", iteration)
  }
}

# For debugging
save.image("tmp/R_ML_image_6_new.rsave")

lodoroc_no_subsampling <- do.call('rbind', resultsPred) %>% 
    left_join(save %>% select(display_dataset, fmt_id) %>% distinct()) %>%
    filter(preFMTAbundance != 0 | donorAbundance != 0) %>%
    group_by(display_dataset, fmt_id) %>%
    nest() %>%
    mutate(auc = map_dbl(data, function(x) suppressMessages(roc(x$truth, x$prob.1)$auc))) %>%
    select(display_dataset, fmt_id, auc)
                         
tmp <- rbind(lodoroc_no_subsampling %>%
                            arrange(display_dataset, fmt_id, auc) %>%
                            mutate(type = "LODO"),
                            cvaurocs_curve_per_FMT_ID %>% 
                            #rename(display_dataset = testDataset) %>%
                            arrange(display_dataset, fmt_id, auc) %>%
                            mutate(type = "CV") %>%
                            select(display_dataset, fmt_id, auc, type) %>%
                            distinct(.keep_all = TRUE)) %>%
                            group_by(display_dataset) %>%
                            nest() %>%
                            mutate(N = map_dbl(data, function(x) return(length(unique(x$fmt_id))))) %>%
                            mutate(display_dataset = map2_chr(display_dataset, N, 
                                                              function(a, n) {
                                                                  return(str_c(a, " (n=", n, ")", sep = "", collapse = ""))
                                                              })) %>%
                            unnest()     
tmp <- tmp %>%                         
                            mutate(display_dataset = factor(display_dataset, 
                                levels = tmp %>%
                                filter(type == "LODO") %>%
                                group_by(display_dataset) %>% 
                                summarize(m = median(auc)) %>% 
                                arrange(desc(m)) %>% 
                                pull(display_dataset)))                                               
pvals <- tmp %>%
ungroup() %>%
select(display_dataset, fmt_id, auc, type)  %>%
pivot_wider(id_cols = c(display_dataset, fmt_id), names_from = type, values_from = auc) %>%
group_by(display_dataset) %>%
nest() %>%
mutate(wilcox_pval = map_dbl(data, function(x) return(wilcox.test(x = x$LODO, y = x$CV)$p.value))) %>%
mutate(wilcox_pval = ifelse(wilcox_pval < 0.05, "*", ""))                                               


aurocboxplot <- ggplot(tmp) +
  geom_boxplot(aes(x = display_dataset, y = auc, fill = type), outlier.color =  NA) +
  #geom_jitter(aes(x = display_dataset, y = auc, group = type), width = 0.1, height = 0, alpha = 0.5) +
  geom_point(aes(x = display_dataset, y = auc, fill = type), position=position_jitterdodge(dodge.width=1, jitter.width = 0.2), alpha = 0.5, size = 0.5) +
  #geom_point(data = cvroc %>% group_by(display_dataset) %>% mutate(auroc = median(auroc)), aes(x = display_dataset, y = auroc), color = '#43a2ca') +
  theme_embl() + 
  ylim(c(0.4, 1)) +
  xlab("Dataset") +
  ylab("AUROC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("#a6cee3", "#b2df8a")) + 
  scale_color_manual(values = c("#a6cee3", "#b2df8a")) +
  geom_text(data = pvals,
           aes(x = display_dataset, y = 1, label = wilcox_pval)) +
  ylim(c(0.5, 1))

#ggsave(plot = aurocboxplot, filename = "figures/aurocboxplot.pdf", width = 6, height = 3)
ggsave(plot = aurocboxplot, filename = "figures/Figure_4A_boxplots.pdf", width = 6.5, height = 3)
options(repr.plot.width=8, repr.plot.height=5)


# Get AUCs without removing double 0s (for reviewer 1).
lodorocs <- do.call('rbind', resultsPred) %>% 
    left_join(save %>% select(display_dataset, fmt_id) %>% distinct()) %>%
    #filter(preFMTAbundance != 0 | donorAbundance != 0) %>%
    group_by(display_dataset, fmt_id) %>%
    nest() %>%
    mutate(auc = map_dbl(data, function(x) suppressMessages(roc(x$truth, x$prob.1)$auc))) %>%
                     select(-data)
                         #head()
lodorocs$display_dataset <- factor(as.vector(lodorocs$display_dataset),
                                  levels = lodorocs %>%
                                  group_by(display_dataset) %>%
                                  summarize(m = mean(auc)) %>%
                                  arrange(desc(m)) %>%
                                  pull(display_dataset))  
options(repr.plot.width=7.5, repr.plot.height=3.25)  
lodorocs %>% ggplot(aes(x = display_dataset, y = auc)) + geom_boxplot() +
theme_embl() +
theme(axis.text.x = element_text(angle = 45, hjust = 1))                          

lodorocs <- do.call('rbind', resultsPred) %>% 
    left_join(save %>% select(display_dataset, fmt_id) %>% distinct()) %>%
    filter(preFMTAbundance != 0 | donorAbundance != 0) %>%
    group_by(display_dataset, fmt_id) %>%
    nest() %>%
    mutate(roc = map(data, function(x) suppressMessages(roc(x$truth, x$prob.1)))) %>%
    select(display_dataset, fmt_id, roc)  %>% 
    mutate(auc = map2(roc, display_dataset, function(x, ds) {
        tmp <- (data.frame(x$sensitivities, 
                           x$specificities, 
                           ds, 
                           #length(test_set), 
                           round(x$auc, 3))) 
                           #iteration = i, 
                           #pon = 1:length(x$sensitivities)))
        colnames(tmp) <- c("Sensitivity", "Specificity", "display_dataset", 'auc')
        return(tmp)
      })) %>%
    select(fmt_id, auc) %>%
    unnest() %>%
    select(Sensitivity, Specificity, display_dataset, fmt_id) %>%
      ggplot(aes(x = 1-Sensitivity, Specificity, group = fmt_id)) + 
      geom_line(alpha = 0.1) + 
      theme_embl() +
      xlab("FPR") + 
      ylab("TPR")
options(repr.plot.width=9.5, repr.plot.height=5)                     
# This figure isn't part of the final manuscript.
ggsave(plot = lodorocs + cvaurocs, 
       filename = "figures/ROC_curve_per_FMT_ID_new_LODO_plus_CV.pdf", width = 5, height = 2.5)                     





featImpMat <- tibble(data = feature_importances,
                    label = names(feature_importances))
#ss <- featImpMat
featImpMat$data <- map(featImpMat$data, function(x) return(x[[1]]))
featImpMat <- featImpMat %>% mutate(data = map(data, function(x) {
  return(data.frame(featureName = names(x),
                    feature = x))
}))
featImpMat <- featImpMat %>% tidyr::unnest()
featImpMat$dataset <- map_chr(featImpMat$label, function(x) str_split(x, "__")[[1]][1])
featImpMat$iteration <- map_chr(featImpMat$label, function(x) str_split(x, "__")[[1]][2])

featImpMat$featureName[featImpMat$featureName == "Antibiotics_parsed__minimum_days_since_last_dose_"] <- "antibiotics"
featImpMat$featureName[featImpMat$featureName == "donorTopreFMTAbundanceRatio"] <- "donorPreFMTRatio"
featImpMat$featureName[featImpMat$featureName == "Coccus___pairs_or_chains_predominate"] <- "coccus__pairs"
featImpMat$featureName[featImpMat$featureName == "Stool_W"] <- "prevalence"
featImpMat$featureName[featImpMat$featureName == "alpha_div_pre"] <- "alphaDiversityPreFMT"
featImpMat$featureName[featImpMat$featureName == "alpha_div_donor"] <- "alphaDiversityDonor"                                
                                
featImpMatR <- featImpMat %>%
            rename(m = feature)
featImpMat <- featImpMat %>% 
  group_by(featureName) %>% 
  summarize(sdd = sd(feature), m = mean(feature))
# Fix some feature names for easier display


#featImpMat$featureName[featImpMat$featureName == "DonorPreRatio"] <- "donor_pre_fmt_ratio"
featImpMat$featureName <- factor(featImpMat$featureName, levels = featImpMat %>% arrange(desc(m)) %>% pull(featureName))
featImpMat$featureName <- factor(as.vector(featImpMat$featureName), levels = rev(levels(featImpMat$featureName)))

# #####################################################
# # Important! I'm removing lowly associated features #
# #####################################################
# featImpMat <- featImpMat %>% filter(m >= 0.003)

groups <- list(c('prevalence',
                 "phylum", 
                 'class',
                 'order',
                 'family',
                 'genus',
                 'species'), # species information (TODO: Add rest)
               c("preFMTAbundance", 'alphaDiversityPreFMT'), # preFMT
               c("donorAbundance", "alphaDiversityDonor"), # donor 
               c("donorPreFMTRatio", "jaccardDistance", "brayDistance", "time_point"),
               c("display_dataset",
                 'infectious_disease')) # instance-specific  
names(groups) <- c("Species", "Pre-FMT", "Donor", "Instance", 'Cohort')
# add rest
#groups[[length(groups) + 1]] <- featImpMat$featureName[! featImpMat$featureName %in% unlist(groups)]
#names(groups)[length(groups)] <- 'rest'
n <- names(groups)
groups <- do.call('rbind', map2(groups, names(groups), function(a, b){
  return(data.frame(group = b,
                    featureName = a))
}))
featImpMat <- left_join(featImpMat, groups, by = 'featureName')
featImpMat$group <- factor(as.vector(featImpMat$group), levels = n)
                                
featImpMat <- featImpMat %>% 
  group_by(group) %>% 
  nest() %>% 
  mutate(data = map(data, function(x) {
    return(x %>% mutate(featureName = factor(featureName, levels = x %>% arrange(desc(m)) %>% pull(featureName))))
    })) %>% unnest()

featImpMat <- featImpMat %>% arrange(m)
featImpMat$featureName <- factor(featImpMat$featureName, levels = featImpMat$featureName)

featImpMatR %>% group_by(featureName) %>% tally()

featImpMatR <- left_join(featImpMatR, groups, by = 'featureName')
featImpMatR$group <- factor(featImpMatR$group, levels = levels(featImpMat$group))
featImpMatR$featureName <- factor(featImpMatR$featureName, levels = levels(featImpMat$featureName))                                    

p3 <- ggplot() + 
  geom_bar(data = featImpMat, aes(x = featureName, y = m), stat = 'identity', alpha = 0.5)  + 
  geom_jitter(data = featImpMatR, aes(x = featureName, y = m), height = 0,  alpha = 0.4, size = 0.1) + 
  geom_errorbar(data = featImpMat, aes(ymin=m-sdd, ymax=m+sdd, x = featureName, y = m), width=.2, alpha = 0.5)
p3 <- p3 + theme_embl()
p3 <- p3 + xlab("Feature") + ylab("Feature Importance")
p3 <- p3 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p3 <- p3 + facet_grid(.~group, scales = 'free_x', space = 'free')
#p3 <- p3 + coord_flip()
#p3 <- p3 + theme(legend.position = 'bottom',
                 #axis.text.x = element_text(angle = 45, hjust = 1))
p3 <- p3 + guides(fill = guide_legend(nrow = 3))
#ggsave(plot = p3, filename = "figures/feature_imporance_all.pdf", width = 5.75, height = 2.6)
ggsave(plot = p3, filename = "figures/Figures_4B.pdf", width = 5.75, height = 2.6)

# save.image(str_c("2__", version[1], "__", version[2], ".rdata"))

###
###  Evaluation of simulations
###
# Keep only prevalent species
#prof <- read_csv("/home/nicolai/mp_CM/CIBIO-Storage/CM/scratch/projects/nkarcher_FMT_meresults/general_population_profiles.csv")
prof <- read_csv("data/prevalences_fmt_byenv.tsv")

prof <- as.data.frame(prof)
rownames(prof) <- prof$clade
prof$clade <- NULL
prevs <- apply(prof, 1, function(x) mean(x > 0))
prevalentSpecies <- names(prevs)[prevs > 0.1]





# Get confusion matrices/accuracy for classifier
CMs <- map(resultsPred, function(x){
    # At some point in the scripts above a factor with levels c(0, 1)
    # Gets cast to numeric which turns it into c(1,2). So turn back for clarity.
    #x$truth <- ifelse(x$truth == 2, 1, 0)
    #x$response <- ifelse(x$response == 2, 1, 0)
    x <- x %>%
    filter(preFMTAbundance != 0 | donorAbundance != 0)
    tt <- confusionMatrix(factor(x %>% pull(truth), levels = c(0,1)),
               factor(x %>% pull(response), levels = c(0,1)))
    return(list(tt$table, tt$overall[['Accuracy']]))
    
})





CMs <- tibble(confTab = map(CMs, function(x) x[[1]] %>% as.data.frame()),
              accuracy = map(CMs, function(x) x[[2]]))
CMs$ds <- (save %>% pull(display_dataset) %>% unique())                        
CMs <- CMs %>%
mutate(plot = pmap(list(confTab, accuracy, ds), function(a, b, d) {
  ggplot(data =  a %>%
         rename(Truth = Reference), mapping = aes(x = Truth, y = Prediction)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = round(Freq/sum(a$Freq), 3)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  theme_embl() + 
    theme(legend.position = "none") +
    ggtitle(str_c("Dataset: ", d, "\nAccuracy: ", round(b, 3))) %>%
    return()
}))
options(repr.plot.width=14, repr.plot.height=17)
ggsave(plot = wrap_plots(CMs$plot, ncol = 5), filename = "figures/ED_Figure_7.pdf", width = 14, height = 17)
wrap_plots(CMs$plot, ncol = 5)

# Visualize regression model performance


b <- do.call('rbind', abundancePredictions) %>%
left_join(engraftmentData %>% 
         select(display_dataset, fmt_id) %>%
         distinct(.keep_all = TRUE), by = 'fmt_id')
#filter(display_dataset == "GollR_2020")
b <- b %>% group_by(ds) %>% nest() %>% mutate(plots = map2(data, ds, function(a, dss) {
    a <- a %>%
    filter(preFMTAbundance > 0 | donorAbundance > 0)
    tmp <- a
    tmp$truth <- log10(tmp$truth + 1E-5)
    tmp$response <- log10(tmp$response + 1E-5)
    linFit <- lm(response ~ truth, data = tmp)
    interc <- linFit$coefficients[1]
    slope <- linFit$coefficients[2]
    p <- ggplot(a %>%
                mutate(prediction = response), aes(x = log10(truth+1E-5), y = log10(prediction +1E-5))) + 
    geom_point(alpha = 0.1) +
    geom_abline(slope = 1, intercept = 0, linetype = 'solid') +
    #geom_abline(slope = slope, intercept = interc, linetype = 'dashed') +
    #geom_smooth(method = 'lm') + 
    theme_embl() +
    ggtitle(dss) +
    annotate(geom = "text", x = -1.1, y = -3, label = str_c("Spearman cor.: ", round(cor(log10(a$truth+1E-5), 
                                                                     log10(a$response+1E-5), 
                                                                     method = "spearman"), 3))) +
    annotate(geom = "text", x = -0.925, y = -3.5, label = str_c("Spearman cor. *: ", round(cor(log10(a %>% 
                                                                                                filter(truth != 0 &
                                                                                                      response != 0) %>%
                                                                                               mutate(truth = truth + 1E-5,
                                                                                                     response = response + 1E-5) %>%
                                                                                               pull(truth)), 
                                                                                          log10(a %>% 
                                                                                                filter(truth != 0 &
                                                                                                      response != 0) %>%
                                                                                               mutate(truth = truth + 1E-5,
                                                                                                     response = response + 1E-5) %>%
                                                                                               pull(response)),
                                                                     method = "spearman"), 3))) +    
    xlim(c(-5, 1)) + 
    ylim(c(-5, 1)) +
    #annotate(geom = "text", x = -3.5, y = 1.1, label = str_c("R2: ", round(summary(linFit)$r.squared, 3)))    
    NULL
    xdens <- axis_canvas(p, axis = "x") +
    geom_histogram(data = a, 
                aes(x = log10(truth+1E-5)),
                alpha = 0.7, size = 0.15, bins = 30) 
    ydens <- axis_canvas(p, axis = "y", coord_flip = TRUE) +
    geom_histogram(data = a, 
                aes(x = log10(response+1E-5)),
                alpha = 0.7, size = 0.15, bins = 30)  +
    coord_flip()
    #facet_wrap(.~ds, nrow = 5) +
    #geom_text(data = a %>% filter(ds == "VermaS_2021")
    #          group_by(ds) %>%
    #          summarize(c = round(cor(log10(truth+1E-6), log10(response + 1E-6)), 3)),
    #         aes(x = -3, y = 1, label = c), inherit.aes = F)
    options(repr.plot.width=16, repr.plot.height=16)
    p1 <- insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
    p1 <- insert_yaxis_grob(p1, ydens, grid::unit(0.2, "null"), position = 'right')
    return(ggdraw(p1))
}))
ggsave(plot = wrap_plots(b$plots, ncol = 5), filename = "figures/Figure_4G.pdf", width = 16, height = 16)
options(repr.plot.width=16, repr.plot.height=16)
wrap_plots(b$plots, ncol = 5)





# Visualize regression model perfor## Here I'm only plotting one big plot and I downsample points to equal depth for fair evaluation.
options(scipen=999)
set.seed(1)
b <- do.call('rbind', abundancePredictions) %>%
    left_join(engraftmentData %>% 
             select(display_dataset, fmt_id) %>%
             distinct(.keep_all = TRUE), by = 'fmt_id')
datasetCounts <- b %>%
     filter(preFMTAbundance > 0 | donorAbundance > 0) %>%
pull(display_dataset) %>% 
table() %>%
as.data.frame() %>%
arrange(Freq)
p2 <- b %>% nest() %>% mutate(plots = map2(data, ds, function(a, dss) {
    a <- a %>%
    filter(preFMTAbundance > 0 | donorAbundance > 0) %>%
    group_by(display_dataset) %>%
    nest() %>%
    mutate(data = map(data, function(x) x %>%
                      sample_n(datasetCounts$Freq[1]) %>%
                     return())) %>%
    unnest()
    tmp <- a
    tmp$truth <- log10(tmp$truth + 1E-5)
    tmp$response <- log10(tmp$response + 1E-5)
    linFit <- lm(response ~ truth, data = tmp)
    interc <- linFit$coefficients[1]
    slope <- linFit$coefficients[2]
    p <- ggplot() + 
    geom_point(data = a %>%
               mutate(prediction = response) %>%
               mutate(`Post-FMT species\nrelative abundance (%)` = truth+1E-5,
                     `predicted post-FMT species\nrelative abundance (%)` = prediction+1E-5), 
               aes(x = `Post-FMT species\nrelative abundance (%)`, y = `predicted post-FMT species\nrelative abundance (%)`), 
               alpha = 0.05) +
    geom_abline(slope = 1, intercept = 0, linetype = 'solid') +
    #geom_abline(slope = slope, intercept = interc, linetype = 'dashed') +
    #geom_smooth(method = 'lm') + 
    theme_embl() +
    #ggtitle(dss) +
    annotate(geom = "text", x = 0.1, y = 0.001, label = str_c("Spearman cor.: ", round(cor(log10(a$truth+1E-5), 
                                                                     log10(a$response+1E-5), 
                                                                     method = "spearman"), 3))) +
    annotate(geom = "text", x = 0.1, y = 0.0001, label = str_c("Spearman cor. *: ", round(cor(log10(a %>% 
                                                                                                filter(truth != 0 &
                                                                                                      response != 0) %>%
                                                                                               mutate(truth = truth + 1E-5,
                                                                                                     response = response + 1E-5) %>%
                                                                                               pull(truth)), 
                                                                                          log10(a %>% 
                                                                                                filter(truth != 0 &
                                                                                                      response != 0) %>%
                                                                                               mutate(truth = truth + 1E-5,
                                                                                                     response = response + 1E-5) %>%
                                                                                               pull(response)),
                                                                     method = "spearman"), 3))) +    

    scale_x_log10(breaks = c(1E-4, 1E-2, 1E0,1E2), 
                  #labels = comma,
                  limits = c(1E-5, 100)) + 
    scale_y_log10(breaks = c(1E-4, 1E-2, 1E0,1E2), 
                  #labels = comma,
                  limits = c(1E-5, 100)) + 
    #annotate(geom = "text", x = -3.5, y = 1.1, label = str_c("R2: ", round(summary(linFit)$r.squared, 3)))    
    NULL
    xdens <- axis_canvas(p, axis = "x") +
    geom_histogram(data = a %>% 
                   mutate(prediction = response) %>%
                   mutate(`Post-FMT species abundance (%)` = truth+1E-5,
                   `predicted post-FMT species abundance (%)` = prediction+1E-5), 
                aes(x = `Post-FMT species abundance (%)`),
                alpha = 0.7, size = 0.15, bins = 30) +
    scale_x_log10()
    #scale_y_log10()
    ydens <- axis_canvas(p, axis = "y", coord_flip = TRUE) +
    geom_histogram(data = a %>% 
                   mutate(prediction = response) %>%
                   mutate(`Post-FMT species abundance (%)` = truth+1E-5,
                   `predicted post-FMT species abundance` = prediction+1E-5), 
                aes(x = `predicted post-FMT species abundance`), 
                alpha = 0.7, size = 0.15, bins = 30)  +
    coord_flip() +
    scale_x_log10() 
    #scale_y_log10()
    #facet_wrap(.~ds, nrow = 5) +
    #geom_text(data = a %>% filter(ds == "VermaS_2021")
    #          group_by(ds) %>%
    #          summarize(c = round(cor(log10(truth+1E-6), log10(response + 1E-6)), 3)),
    #         aes(x = -3, y = 1, label = c), inherit.aes = F)
    #options(repr.plot.width=16, repr.plot.height=16)
    p1 <- insert_xaxis_grob(p, xdens, grid::unit(.2, "null"), position = "top")
    p1 <- insert_yaxis_grob(p1, ydens, grid::unit(0.2, "null"), position = 'right')
    return(ggdraw(p1))
}))
options(repr.plot.width=4, repr.plot.height=4)
p2$plots[[1]]
# This figure isn'r part of the final manuscript
ggsave(plot = p2$plots[[1]], filename = "figures/abundancePredictionScatterNonZeroPreFMTAbundanceORNonZeroDonorAbundanceOnePlot.pdf", width = 3.6, height = 3.25)
#ggsave(plot = wrap_plots(b$plots, ncol = 5), filename = "figures/abundancePredictionScatterNonZeroPreFMTAbundanceORNonZeroDonorAbundance.pdf", width = 16, height = 16)

tmp <- richnessPredictedUnalteredInstance
richnessPredictedUnalteredInstance <- do.call('rbind', richnessPredictedUnalteredInstance) 
abundancePredictedUnalteredInstance <- do.call("rbind", abundancePredictedUnalteredInstance)
results <- do.call('rbind', results)

realValues <- save %>% 
  filter(dummyData == "FALSE") %>% 
  group_by(fmt_id, display_dataset) %>%
  nest() %>%
  summarize(postFMTRichness = map_dbl(data, function(x) sum(x$postFMTAbundance > 0) %>% return()), 
            PreFMTRichness = map_dbl(data, function(x) sum(x$preFMTAbundance > 0) %>% return()),
            DonorRichness = map_dbl(data, function(x) sum(x$donorAbundance > 0) %>% return()),
            proteoPrevalence = map_dbl(data, function(x) x %>% filter(phylum == "p__Proteobacteria") %>%
                                       mutate(postFMTAbundance = postFMTAbundance > 0) %>%
                                      pull(postFMTAbundance) %>%
                                      sum() %>%
                                      return()), 
            donorProteoPrevalence = map_dbl(data, function(x) x %>% filter(phylum == "p__Proteobacteria") %>%
                                       mutate(donorAbundance = donorAbundance > 0) %>%
                                      pull(donorAbundance) %>%
                                      sum() %>%
                                      return()), 
            proteoAbundance = map_dbl(data, function(x) x %>% filter(phylum == "p__Proteobacteria") %>%
                                       #mutate(postFMTAbundance = postFMTAbundance > 0) %>%
                                      pull(postFMTAbundance) %>%
                                      sum() %>%
                                      return()),
            donorProteoAbundance = map_dbl(data, function(x) x %>% filter(phylum == "p__Proteobacteria") %>%
                                      #mutate(donorAbundance = donorAbundance > 0) %>%
                                      pull(donorAbundance) %>%
                                      sum() %>%
                                      return()), 
            FirmicuteAbundance = map_dbl(data, function(x) x %>% filter(phylum == "p__Firmicutes") %>%
                                       #mutate(postFMTAbundance = postFMTAbundance > 0) %>%
                                      pull(postFMTAbundance) %>%
                                      sum() %>%
                                      return()),
            donorFirmicuteAbundance = map_dbl(data, function(x) x %>% filter(phylum == "p__Firmicutes") %>%
                                      #mutate(donorAbundance = donorAbundance > 0) %>%
                                      pull(donorAbundance) %>%
                                      sum() %>%
                                      return()),      
            FirmicutePrevalence = map_dbl(data, function(x) x %>% filter(phylum == "p__Firmicutes") %>%
                                      mutate(postFMTAbundance = postFMTAbundance > 0) %>%
                                      pull(postFMTAbundance) %>%
                                      sum() %>%
                                      return()),
            donorFirmicutePrevalence = map_dbl(data, function(x) x %>% filter(phylum == "p__Firmicutes") %>%
                                      mutate(donorAbundance = donorAbundance > 0) %>%
                                      pull(donorAbundance) %>%
                                      sum() %>%
                                      return()),          
            BacteroidetesAbundance = map_dbl(data, function(x) x %>% filter(phylum == "p__Bacteroidetes") %>%
                                       #mutate(postFMTAbundance = postFMTAbundance > 0) %>%
                                      pull(postFMTAbundance) %>%
                                      sum() %>%
                                      return()),
            donorBacteroidetesAbundance = map_dbl(data, function(x) x %>% filter(phylum == "p__Bacteroidetes") %>%
                                      #mutate(donorAbundance = donorAbundance > 0) %>%
                                      pull(donorAbundance) %>%
                                      sum() %>%
                                      return()),      
            BacteroidetesPrevalence = map_dbl(data, function(x) x %>% filter(phylum == "p__Bacteroidetes") %>%
                                      mutate(postFMTAbundance = postFMTAbundance > 0) %>%
                                      pull(postFMTAbundance) %>%
                                      sum() %>%
                                      return()),
            donorBacteroidetesPrevalence = map_dbl(data, function(x) x %>% filter(phylum == "p__Bacteroidetes") %>%
                                      mutate(donorAbundance = donorAbundance > 0) %>%
                                      pull(donorAbundance) %>%
                                      sum() %>%
                                      return()),                                         
            predict1Prevalence = map_dbl(data, function(x) x %>%
                                         filter(species %in% topPredict1RankSpecies) %>%
                                         mutate(postFMTAbundance = postFMTAbundance > 0) %>%
                                         pull(postFMTAbundance) %>%
                                         sum() %>%
                                         return()),  
            donorPredict1Prevalence = map_dbl(data, function(x) x %>%
                                         filter(species %in% topPredict1RankSpecies) %>%
                                         mutate(donorAbundance = donorAbundance > 0) %>%
                                         pull(donorAbundance) %>%
                                         sum() %>%
                                         return()),                                    
            predict1Abundance = map_dbl(data, function(x) x %>%
                                         filter(species %in% topPredict1RankSpecies) %>%
                                         #mutate(postFMTAbundance = postFMTAbundance > 0) %>%
                                         pull(postFMTAbundance) %>%
                                         sum() %>%
                                         return()), 
            donorPredict1Abundance = map_dbl(data, function(x) x %>%
                                         filter(species %in% topPredict1RankSpecies) %>%
                                         #mutate(donorAbundance = donorAbundance > 0) %>%
                                         pull(donorAbundance) %>%
                                         sum() %>%
                                         return()),                                      
            oralPrevalence =map_dbl(data, function(x) x %>%
                                         filter(species %in% oralAbPrev) %>%
                                         mutate(postFMTAbundance = postFMTAbundance > 0) %>%
                                         pull(postFMTAbundance) %>%
                                         sum() %>%
                                         return()), 
            donorOralPrevalence =map_dbl(data, function(x) x %>%
                                         filter(species %in% oralAbPrev) %>%
                                         mutate(donorAbundance = donorAbundance > 0) %>%
                                         pull(donorAbundance) %>%
                                         sum() %>%
                                         return()),                                    
            oralAbundance = map_dbl(data, function(x) x %>%
                                         filter(species %in% oralAbPrev) %>%
                                         #mutate(postFMTAbundance = postFMTAbundance > 0) %>%
                                         pull(postFMTAbundance) %>%
                                         sum() %>%
                                         return()),
            donorOralAbundance =map_dbl(data, function(x) x %>%
                                         filter(species %in% oralAbPrev) %>%
                                         #mutate(donorAbundance = donorAbundance > 0) %>%
                                         pull(donorAbundance) %>%
                                         sum() %>%
                                         return()))                                   

# Some datasets only two or less donors. Don't show these.
aucsDummyDonor[which(map_lgl(aucsDummyDonor, function(x) length(x) <= 2))] <- NULL
aucsDummyRecipient[which(map_lgl(aucsDummyRecipient, function(x) length(x) <= 2))] <- NULL
aucsDummyBoth[which(map_lgl(aucsDummyBoth, function(x) length(x) <= 2))] <- NULL
 
aucsDummyDonor <- do.call('rbind', map2(names(aucsDummyDonor), aucsDummyDonor, function(a,b) {
  return(data.frame(r = a, auc = map_dbl(b, function(b) b$auc)))
})) %>%
  mutate(type = "dummyDonor") %>% as.data.frame()
aucsDummyRecipient <- do.call('rbind', map2(names(aucsDummyRecipient), aucsDummyRecipient, function(a,b) {
  return(data.frame(r = a, auc = map_dbl(b, function(b) b$auc)))
})) %>%
  mutate(type = "dummyRecipient") %>% as.data.frame()
aucsDummyBoth <- do.call('rbind', map2(names(aucsDummyBoth), aucsDummyBoth, function(a,b) {
  return(data.frame(r = a, auc = map_dbl(b, function(b) b$auc)))
})) %>%
  mutate(type = "dummyBoth") %>% as.data.frame()
dummy <- rbind(aucsDummyRecipient, aucsDummyDonor, aucsDummyBoth) %>%
  mutate(r = map_chr(r, function(x) str_split(x, "__")[[1]][1]))

all <- rbind(dummy, results %>%
               select(testDataset, auc) %>%
               rename(r = testDataset) %>%
               distinct() %>%
               mutate(type = 'real')) %>%
  mutate(auc = as.numeric(auc))
all$type[all$type=="dummyDonor"] <- "Donor exchanged"
all$type[all$type=="dummyRecipient"] <- "Recipient exchanged"
all$type[all$type=="real"] <- "Neither exchanged"
all$type[all$type=="dummyBoth"] <- "Both exchanged"
all <- all %>% rename(Group = type)
# For some reason some aucs are << 0.5. Probably an issue with how the aurocs are computed. Fix like this.
all$auc[all$auc < 0.5] <- 1 - all$auc[all$auc < 0.5]
# ggplot(all, aes(x = r, y = auc, fill = type)) +
#   geom_boxplot(position = 'dodge') +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1 ))
tmp <- all %>% filter(Group == 'Neither exchanged') %>% select(-Group)
all <- all %>% left_join(tmp, by = 'r') %>% mutate(auc_difference = auc.x-auc.y)
ttt <- all %>% group_by(r) %>% do(., {
  tmp <- .
  a <- tmp %>% filter(Group == "Donor exchanged") %>% pull(auc_difference) %>% median()
  b <- tmp %>% filter(Group == "Recipient exchanged") %>% pull(auc_difference) %>% median()
  data.frame(r = tmp$r[1],
             m = a-b)
})
tt <- ttt %>%
  arrange(m) %>%
  pull(r)
all <- all %>% mutate(r = factor(r, levels = tt))
tt <- left_join(all %>% 
                  select(r) %>% 
                  rename(study = r) %>% 
                  filter(study %in%  map_chr(unique(aucsDummyDonor$r), function(x) str_split(x, "__")[[1]][1])), clinicalMeta %>% select(study, infectious_disease)) %>% distinct()
                                                                                
tt$infectious_disease[tt$study == "IaniroG_2020"] <- FALSE
tt$infectious_disease[tt$study == "This_study_Cdiff"] <- TRUE
tt$infectious_disease[tt$study == "This_study_MDR"] <- TRUE
tt$infectious_disease[tt$study == "This_study_IBD"] <- FALSE
tt$infectious_disease[tt$study == "DavarD_2021"] <- FALSE

tt$study <- factor(tt$study, levels = levels(all$r))
# For later...                                             
p2  <- ggplot(tt %>%
                mutate(group = ifelse(as.vector(infectious_disease),
                                      "Infectious disease", 
                                      "Non-infectious disease")), aes(x = study, y = 1, fill = group)) + 
  geom_bar(stat = 'identity') + 
  theme_embl() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = c('#bfd3e6', '#8c6bb1'))

ab_data <- data.frame(r = c("KumarR_2017", #1
                            "ZhaoH_2020", #2
                            "PodlesnyD_2020", #3
                            "MossE_2017",
                            "SmillieC_2018", #5
                            "LeoS_2020",
                            "BaruchE_2020", #7
                            "HouriganS_2019",
                            "This_study_Cdiff", #9
                            "This_study_MDR",
                            "BarYosephH_2020", # 11
                            "GollR_2020",
                            "LiS_2016", #13
                            "This_study_IBD",
                            "SuskindD_2015", #15
                            "IaniroG_2020", 
                            "DammanC_2015",#17
                            "VaughnB_2016", 
                            "KongL_2020",#19
                            "DavarD_2021",
                           "AggarwalaV_2021",
                           "KoopenAM_2021",
                           "VermaS_2021",
                           "WatsonAR_2021"),
                      antibiotics = c("Treatment close to FMT", #1
                                      "No antibiotics", #2 
                                      "Treatment close to FMT", #3
                                      "Treatment close to FMT",
                                      "Treatment close to FMT", #5
                                      "FMT pre-treatment",
                                      "FMT pre-treatment", #7
                                      "Treatment close to FMT",
                                      "FMT pre-treatment", #9
                                      "Treatment close to FMT",
                                      "Treatment close to FMT", #11
                                      "No antibiotics",
                                      "No antibiotics", #13
                                      "No antibiotics",
                                      "FMT pre-treatment",#15
                                      "No antibiotics",
                                      "No antibiotics",#17
                                      "No antibiotics",
                                      "No antibiotics",#19
                                      "No antibiotics",
                                     "FMT pre-treatment",
                                     "No antibiotics",
                                     "Treatment close to FMT",
                                     "FMT pre-treatment"))
ab_data$antibiotics <- factor(ab_data$antibiotics, levels = c('No antibiotics', 
                                                              'Treatment close to FMT',
                                                              'FMT pre-treatment'))
ab_data$r <- factor(ab_data$r, levels = levels(all$r))
ab_data <- ab_data %>% filter(r %in% unique(map_chr(aucsDummyDonor$r, function(x) str_split(x, "__")[[1]][1])))

pa <- all %>%
               filter(Group != "Neither exchanged") %>%
               filter(Group != "Both exchanged") %>%
               filter(auc_difference > -0.15) %>%
left_join(ab_data, by = "r") %>%
filter(Group != "Recipient exchanged") %>%
mutate(antibiotics = ifelse(antibiotics == "No antibiotics", "no antibiotics", "antibiotics")) %>%
mutate(antibiotics = factor(antibiotics, levels = c("no antibiotics", 'antibiotics')))
pa <- left_join(pa, pa %>% group_by(antibiotics) %>% tally()) %>%
mutate(antibiotics = map2_chr(antibiotics, n, function(a, b) str_c(a, " N(", b, ")"))) %>%
ggplot(aes(x = antibiotics, y = auc_difference)) + 
geom_boxplot(outlier.color = NA) +
theme_embl() + 
ylab("Difference in AUC\nupon donor exchange") +
xlab("") +
geom_jitter(alpha = 0.1) +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
annotate('text', x = 1.5, y = 0.05, label = "p<2e-16")

pVal1 <- all %>%
               filter(Group != "Neither exchanged") %>%
               filter(Group != "Both exchanged") %>%
               filter(auc_difference > -0.15) %>%
left_join(ab_data, by = "r") %>%
filter(Group != "Recipient exchanged") %>%
mutate(antibiotics = ifelse(antibiotics == "No antibiotics", "no antibiotics", "antibiotics")) %>%
mutate(antibiotics = factor(antibiotics, levels = c("no antibiotics", 'antibiotics')))
pVal1 <- wilcox.test(x = pVal1 %>% 
                    filter(antibiotics == "antibiotics") %>%
                    pull(auc_difference),
                    y = pVal1 %>%
                    filter(antibiotics == "no antibiotics") %>%
                    pull(auc_difference))

pb <- all %>%
               filter(Group != "Neither exchanged") %>%
               filter(Group != "Both exchanged") %>%
               filter(auc_difference > -0.15) %>%
left_join(tt %>%
                mutate(group = ifelse(as.vector(infectious_disease),
                                      "Infectious disease", 
                                      "Non-infectious disease")) %>%
mutate(r=study), by = "r") %>%
filter(Group != "Recipient exchanged") %>%
#mutate(group = ifelse(antibiotics == "No antibiotics", "no antibiotics", "antibiotics")) %>%
mutate(group = factor(group, levels = c("Non-infectious disease", 'Infectious disease')))
pb <- left_join(pb, pb %>% group_by(group) %>% tally()) %>%
mutate(group = map2_chr(group, n, function(a, b) str_c(a, " N(", b, ")"))) %>%
ggplot(aes(x = group, y = auc_difference)) + 
geom_boxplot(outlier.color = NA) +
theme_embl() + 
ylab("Difference in AUC\nupon donor exchange") +
xlab("") +
geom_jitter(alpha = 0.1) +
theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
annotate('text', x = 1.5, y = 0.05, label = "p<2e-16")


pVal2 <- all %>%
               filter(Group != "Neither exchanged") %>%
               filter(Group != "Both exchanged") %>%
               filter(auc_difference > -0.15) %>%
left_join(tt %>%
                mutate(group = ifelse(as.vector(infectious_disease),
                                      "Infectious disease", 
                                      "Non-infectious disease")) %>%
mutate(r=study), by = "r") %>%
filter(Group != "Recipient exchanged") %>%
#mutate(group = ifelse(antibiotics == "No antibiotics", "no antibiotics", "antibiotics")) %>%
mutate(group = factor(group, levels = c("Non-infectious disease", 'Infectious disease')))
pVal2 <- wilcox.test(x = pVal2 %>% 
                    filter(group == "Infectious disease") %>%
                    pull(auc_difference),
                    y = pVal2 %>%
                    filter(group == "Non-infectious disease") %>%
                    pull(auc_difference))

ggsave(plot = pa + pb, filename = "figures/Figure_4C_raw.pdf",
       width = 7, height = 3.5)
p2.1  <- ggplot(ab_data, aes(x = r, y = 1, fill = antibiotics)) + 
  geom_bar(stat = 'identity') + 
  theme_embl() + 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  scale_fill_manual(values = c('#1f77b4', '#ff7f0e', "#2ca02c"))
all <- all %>% 
                left_join(all %>%
                group_by(r) %>%
                tally(), by = 'r') %>%
                # Since r is a factor this reorders the rows with respect to the levels of r.
             arrange(r) %>%
             mutate(r = map2_chr(r, n, function(a,b) {
                 return(str_c(a, "(n=", b, ")", sep = "", collapse = ""))
             })) %>%
             mutate(r = factor(r, levels = unique(r))) %>%
             select(-n)
p1 <- ggplot( all %>%
               filter(Group != "Neither exchanged") %>%
               filter(Group != "Both exchanged") %>%
               filter(auc_difference > -0.15)
             , aes(x = r, y = auc_difference, fill = Group)) +
  geom_abline(slope = 0, intercept = 0, alpha = 0.2) + 
  geom_boxplot(position = 'dodge', outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width = 0.25), alpha = 0.3, size = 0.005) + 
  theme_embl() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1 )) +
  xlab("") +
  ylab("Difference in AUROC") + 
  theme(legend.position = "top",
        legend.box = 'horizontal') +
  ylim(c(-0.15, 0.06)) +
  geom_text(data = all %>% 
               filter(Group != "Neither exchanged") %>%
               filter(Group != "Both exchanged") %>%
                pivot_wider(id_cols = r, names_from = Group, values_from = auc_difference) %>%
                mutate(wilcox_pval = map2_dbl(`Recipient exchanged`, `Donor exchanged`, function(a,b){
                    return(wilcox.test(x = a, y =b )$p.value)
                })) %>%
                select(r, wilcox_pval) %>%
                mutate(wilcox_pval = ifelse(wilcox_pval < 0.05, "*", "")),
                            aes(x = r, y = 0.05, label = wilcox_pval), inherit.aes = F) +      
  NULL   
options(repr.plot.width=12, repr.plot.height=12)
# This figure isn't part of the final manuscript.
ggsave(filename = "figures/test_probal_sum_new_without_both_exchanged.pdf", plot = p2 / p2.1/ p1 + plot_layout(heights = c(1, 1, 6.5)), width = 7*1.125, height = 3.5*1.125)

############################
### richness dummy donor ###
############################
richnessDummyDonor[which(map_lgl(richnessDummyDonor, function(x) length(x) == 0))] <- NULL
richnessDummyDonor <- do.call('rbind', map2(richnessDummyDonor, names(richnessDummyDonor), function(x, nb) {
  return(data.frame(acceptor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][2]),
                    donor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][3]),
                    postFMTRichness = unlist(x),
                    dataset = nb))
}))
richnessDummyDonor <- richnessDummyDonor %>% relocate(dataset, .after = donor_fmt_id)
rownames(richnessDummyDonor) <- NULL  
richnessDummyDonor <- as_tibble(richnessDummyDonor)
richnessDummyDonor <- inner_join(richnessDummyDonor %>% 
                                  mutate(predictedPostFMTRichness = postFMTRichness) %>%
                                  mutate(acceptor_fmt_id = as.double(as.vector(acceptor_fmt_id))) %>%
                                  select(acceptor_fmt_id,
                                         donor_fmt_id, 
                                         dataset, 
                                         predictedPostFMTRichness), 
                                richnessPredictedUnalteredInstance %>% 
                                  mutate(predictedPostFMTRichnessUnalteredInstance = richness) %>%
                                  select(fmt_id, predictedPostFMTRichnessUnalteredInstance),
                                by = c('acceptor_fmt_id' = 'fmt_id'))
richnessDummyDonor$acceptor_fmt_id <- factor(richnessDummyDonor$acceptor_fmt_id, 
                                             levels = richnessDummyDonor %>% 
                                               group_by(acceptor_fmt_id) %>% 
                                               summarize(m  = mean(predictedPostFMTRichness)) %>% 
                                               arrange(desc(m)) %>%
                                               pull(acceptor_fmt_id))    
    
    
    
richnessDummyDonor <- richnessDummyDonor %>% 
  mutate(differenceRichness = predictedPostFMTRichness - predictedPostFMTRichnessUnalteredInstance)
richnessDummyDonor$donor_fmt_id <- factor(as.vector(richnessDummyDonor$donor_fmt_id), levels = richnessDummyDonor %>% 
  group_by(donor_fmt_id) %>% 
  summarize(differenceRichness = mean(differenceRichness)) %>% 
  arrange(desc(differenceRichness)) %>%
  pull(donor_fmt_id))
    





#require(ggsigni)
options(repr.plot.width=3,repr.plot.height=3)
rbind(richnessDummyDonor %>% 
group_by(dataset) %>% 
filter(predictedPostFMTRichness == max(predictedPostFMTRichness)) %>%
mutate(kind = 'top donors'), 
richnessDummyDonor %>% 
group_by(dataset) %>% 
filter(predictedPostFMTRichness == min(predictedPostFMTRichness)) %>%
mutate(kind = 'bottom donors')) %>%
select(donor_fmt_id, kind) %>%
left_join(realValues %>% 
          select(fmt_id, DonorRichness) %>%
         rename(donor_fmt_id = fmt_id) %>%
         mutate(donor_fmt_id = as.factor(donor_fmt_id))) %>%
mutate(kind = factor(kind, levels = c("top donors",
                                                     "bottom donors"))) %>%
ggplot(aes(x = kind, y = DonorRichness)) +
geom_boxplot() +
theme_embl() +
geom_jitter(height = 0, width = 0.3, alpha = 0.3) +
# geom_signif(
#     comparisons = list(c("top donors", "bottom donors")) #map_signif_level = TRUE, textsize = 6
# ) +
ylim(c(50, 430))
ggsave(filename = "figures/supplementary_figure_14.pdf", width = 3, height = 3)

donorOrder <- richnessDummyDonor %>% 
group_by(donor_fmt_id) %>%
summarize(med = median(differenceRichness)) %>%
arrange(desc(med)) %>%
pull(donor_fmt_id)
pp <- ggplot(richnessDummyDonor %>%
             mutate(donor_fmt_id = factor(as.vector(donor_fmt_id),
                                         levels = donorOrder)), aes(x = donor_fmt_id, y = differenceRichness)) + 
  geom_boxplot(outlier.colour = NA) +
  #geom_point(position = position_jitter(), alpha = 0.5, size = 0.5) +
  geom_jitter(alpha = 0.5, size = 0.2, width = 0.5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ylab("Difference in\npredicted post-FMT\nrichness upon donor exchange") +
  xlab("Donor individual")
pp3 <- ggplot(realValues %>% 
              filter(fmt_id %in% richnessDummyDonor$donor_fmt_id) %>%
              mutate(fmt_id = factor(fmt_id, levels = levels(richnessDummyDonor$donor_fmt_id))), 
              aes(x = fmt_id, 
                  y = DonorRichness)) + 
  geom_point(size = 1) + 
  geom_smooth(aes(group = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("richness of\nimplanted donor") +
  ylim(c(0, max(realValues$DonorRichness) * 1.05))
ggsave(plot = pp3/ pp + plot_layout(heights = c(1.25,3.75)), filename = "figures/Figure4D.pdf", width = 7.25, height = 4.1)
options(repr.plot.width=14, repr.plot.height=8)    

########################################################
### cumulative prevalence proteobacteria dummy donor ###
########################################################
CumprevProteoDummyDonor[which(map_lgl(CumprevProteoDummyDonor, function(x) length(x) == 0))] <- NULL
CumprevProteoDummyDonor <- do.call('rbind', map2(CumprevProteoDummyDonor, names(CumprevProteoDummyDonor), function(x, nb) {
  return(data.frame(acceptor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][2]),
                    donor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][3]),
                    postFMTProteoCumPrev = unlist(x),
                    dataset = nb))
}))
CumprevProteoDummyDonor[which(map_lgl(CumprevProteoDummyDonor, function(x) length(x) == 0))] <- NULL
CumprevProteoDummyDonor <- CumprevProteoDummyDonor %>% relocate(dataset, .after = donor_fmt_id)   
rownames(CumprevProteoDummyDonor) <- NULL                                      
CumprevProteoDummyDonor <- as_tibble(CumprevProteoDummyDonor)
CumprevProteoDummyDonor <- inner_join(CumprevProteoDummyDonor %>% 
                                  mutate(predictedPostFMTProteoPrev = postFMTProteoCumPrev) %>%
                                  mutate(acceptor_fmt_id = as.double(as.vector(acceptor_fmt_id))) %>%
                                  select(acceptor_fmt_id,
                                         donor_fmt_id, 
                                         dataset, 
                                         predictedPostFMTProteoPrev), 
                                richnessPredictedUnalteredInstance %>% 
                                  mutate(predictedCumPrevalenceProteobacteriaUnalteredInstance = cumPrevalenceProteobacteria) %>%
                                  select(fmt_id, predictedCumPrevalenceProteobacteriaUnalteredInstance),
                                by = c('acceptor_fmt_id' = 'fmt_id'))
CumprevProteoDummyDonor <- CumprevProteoDummyDonor %>% 
  mutate(differenceProteoPrevalence = predictedPostFMTProteoPrev - predictedCumPrevalenceProteobacteriaUnalteredInstance)
CumprevProteoDummyDonor$donor_fmt_id <- factor(as.vector(CumprevProteoDummyDonor$donor_fmt_id), levels = CumprevProteoDummyDonor %>% 
  group_by(donor_fmt_id) %>% 
  summarize(differenceProteoPrevalence = mean(differenceProteoPrevalence)) %>% 
  arrange(desc(differenceProteoPrevalence)) %>%
  pull(donor_fmt_id)) 





donorOrder <- CumprevProteoDummyDonor %>% 
group_by(donor_fmt_id) %>%
summarize(med = median(differenceProteoPrevalence)) %>%
arrange(desc(med)) %>%
pull(donor_fmt_id)
pp <- ggplot(CumprevProteoDummyDonor %>%
             mutate(donor_fmt_id = factor(as.vector(donor_fmt_id),
                                         levels = donorOrder)), aes(x = donor_fmt_id, y = differenceProteoPrevalence)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ylab("Difference in\npredicted post-FMT\nProteobacterial prevalence upon donor exchange") +
  xlab("Donor individual")
pp3 <- ggplot(realValues %>% 
              filter(fmt_id %in% CumprevProteoDummyDonor$donor_fmt_id) %>%
              mutate(fmt_id = factor(fmt_id, levels = levels(CumprevProteoDummyDonor$donor_fmt_id))), 
              aes(x = fmt_id, 
                  y = donorProteoPrevalence)) + 
  geom_point(size = 1) + 
  geom_smooth(aes(group = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("Proteobacterial prevalence of\nimplanted donor")+
  ylim(c(0, max(realValues$donorProteoPrevalence) * 1.05))
options(repr.plot.width=14, repr.plot.height=8)    
# This figure isn't part of the final manuscript
ggsave(plot = pp, filename = "figures/BestDonorsPlot1_ProteoBacteriaPrev.pdf", width = 4.25, height = 2.1)

###########################################
### cumulative prevalence Predict1 bact ###
###########################################    
CumprevTopPredict1DummyDonor[which(map_lgl(CumprevTopPredict1DummyDonor, function(x) length(x) == 0))] <- NULL
CumprevTopPredict1DummyDonor <- do.call('rbind', map2(CumprevTopPredict1DummyDonor, names(CumprevTopPredict1DummyDonor), function(x, nb) {
  return(data.frame(acceptor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][2]),
                    donor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][3]),
                    postFMTPredict1CumPrev = unlist(x),
                    dataset = nb))
}))
CumprevTopPredict1DummyDonor[which(map_lgl(CumprevTopPredict1DummyDonor, function(x) length(x) == 0))] <- NULL
CumprevTopPredict1DummyDonor <- CumprevTopPredict1DummyDonor %>% relocate(dataset, .after = donor_fmt_id)   
rownames(CumprevTopPredict1DummyDonor) <- NULL                                      
CumprevTopPredict1DummyDonor <- as_tibble(CumprevTopPredict1DummyDonor)
CumprevTopPredict1DummyDonor <- inner_join(CumprevTopPredict1DummyDonor %>% 
                                  mutate(predictedPostFMTPredict1Prev = postFMTPredict1CumPrev) %>%
                                  mutate(acceptor_fmt_id = as.double(as.vector(acceptor_fmt_id))) %>%
                                  select(acceptor_fmt_id,
                                         donor_fmt_id, 
                                         dataset, 
                                         predictedPostFMTPredict1Prev), 
                                richnessPredictedUnalteredInstance %>% 
                                  mutate(predictedCumPrevalencePredict1UnalteredInstance = cumPrevalencePredict1) %>%
                                  select(fmt_id, predictedCumPrevalencePredict1UnalteredInstance),
                                by = c('acceptor_fmt_id' = 'fmt_id'))
CumprevTopPredict1DummyDonor <- CumprevTopPredict1DummyDonor %>% 
  mutate(differencePredict1Prevalence = predictedPostFMTPredict1Prev - predictedCumPrevalencePredict1UnalteredInstance)
CumprevTopPredict1DummyDonor$donor_fmt_id <- factor(as.vector(CumprevTopPredict1DummyDonor$donor_fmt_id), 
                                                    levels = CumprevTopPredict1DummyDonor %>% 
  group_by(donor_fmt_id) %>% 
  summarize(differencePredict1Prevalence = mean(differencePredict1Prevalence)) %>% 
  arrange(desc(differencePredict1Prevalence)) %>%
  pull(donor_fmt_id))   
                                           





donorOrder <- CumprevTopPredict1DummyDonor %>% 
group_by(donor_fmt_id) %>%
summarize(med = median(differencePredict1Prevalence)) %>%
arrange(desc(med)) %>%
pull(donor_fmt_id)
pp <- ggplot(CumprevTopPredict1DummyDonor %>%
             mutate(donor_fmt_id = factor(as.vector(donor_fmt_id),
                                         levels = donorOrder)), aes(x = donor_fmt_id, y = differencePredict1Prevalence)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ylab("Difference in\npredicted post-FMT\nPredict1 bacterial prevalence upon donor exchange") +
  xlab("Donor individual")
pp3 <- ggplot(realValues %>% 
              filter(fmt_id %in% CumprevTopPredict1DummyDonor$donor_fmt_id) %>%
              mutate(fmt_id = factor(fmt_id, levels = levels(CumprevTopPredict1DummyDonor$donor_fmt_id))), 
              aes(x = fmt_id, 
                  y = donorPredict1Prevalence)) + 
  geom_point(size = 1) + 
  geom_smooth(aes(group = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("Predict1 prevalence of\nimplanted donor")+
  ylim(c(0, max(realValues$donorPredict1Prevalence) * 1.05))
options(repr.plot.width=14, repr.plot.height=8)    
# This figure isn't part of the final manuscript
ggsave(plot = pp, filename = "figures/BestDonorsPlot1_Pedict1Prev.pdf", width = 4.25, height = 2.1)

###########################################
### cumulative prevalence Oral bact ###
###########################################    
CumprevOralDummyDonor[which(map_lgl(CumprevOralDummyDonor, function(x) length(x) == 0))] <- NULL
CumprevOralDummyDonor <- do.call('rbind', map2(CumprevOralDummyDonor, names(CumprevOralDummyDonor), function(x, nb) {
  return(data.frame(acceptor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][2]),
                    donor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][3]),
                    postFMTOralCumPrev = unlist(x),
                    dataset = nb))
}))
CumprevOralDummyDonor[which(map_lgl(CumprevOralDummyDonor, function(x) length(x) == 0))] <- NULL
CumprevOralDummyDonor <- CumprevOralDummyDonor %>% relocate(dataset, .after = donor_fmt_id)   
rownames(CumprevOralDummyDonor) <- NULL                                      
CumprevOralDummyDonor <- as_tibble(CumprevOralDummyDonor)
CumprevOralDummyDonor <- inner_join(CumprevOralDummyDonor %>% 
                                  mutate(predictedPostFMTOralPrev = postFMTOralCumPrev) %>%
                                  mutate(acceptor_fmt_id = as.double(as.vector(acceptor_fmt_id))) %>%
                                  select(acceptor_fmt_id,
                                         donor_fmt_id, 
                                         dataset, 
                                         predictedPostFMTOralPrev), 
                                richnessPredictedUnalteredInstance %>% 
                                  mutate(predictedCumPrevalenceOralUnalteredInstance = cumPrevalenceOral) %>%
                                  select(fmt_id, predictedCumPrevalenceOralUnalteredInstance),
                                by = c('acceptor_fmt_id' = 'fmt_id'))
CumprevOralDummyDonor <- CumprevOralDummyDonor %>% 
  mutate(differenceOralPrevalence = predictedPostFMTOralPrev - predictedCumPrevalenceOralUnalteredInstance)
CumprevOralDummyDonor$donor_fmt_id <- factor(as.vector(CumprevOralDummyDonor$donor_fmt_id), levels = CumprevOralDummyDonor %>% 
  group_by(donor_fmt_id) %>% 
  summarize(differenceOralPrevalence = mean(differenceOralPrevalence)) %>% 
  arrange(desc(differenceOralPrevalence)) %>%
  pull(donor_fmt_id))               





donorOrder <- CumprevOralDummyDonor %>% 
group_by(donor_fmt_id) %>%
summarize(med = median(differenceOralPrevalence)) %>%
arrange(desc(med)) %>%
pull(donor_fmt_id)
pp <- ggplot(CumprevOralDummyDonor %>%
             mutate(donor_fmt_id = factor(as.vector(donor_fmt_id),
                                         levels = donorOrder)), aes(x = donor_fmt_id, y = differenceOralPrevalence)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ylab("Difference in\npredicted post-FMT\nOral bacterial prevalence upon donor exchange") +
  xlab("Donor individual")
pp3 <- ggplot(realValues %>% 
              filter(fmt_id %in% CumprevOralDummyDonor$donor_fmt_id) %>%
              mutate(fmt_id = factor(fmt_id, levels = levels(CumprevOralDummyDonor$donor_fmt_id))), 
              aes(x = fmt_id, 
                  y = donorOralPrevalence)) + 
  geom_point(size = 1) + 
  geom_smooth(aes(group = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("prevalence of\nOral bacteria in implanted donor")+
  ylim(c(0, max(realValues$donorOralPrevalence) * 1.05))
options(repr.plot.width=14, repr.plot.height=8)    
# This figure isn't part of the final manuscript
ggsave(plot = pp, filename = "figures/BestDonorsPlot1_OralPrev.pdf", width = 4.25, height = 2.1)

########################################
### cumulative prevalence Firmicutes ###
########################################
cumPrevFirmicutesDummyDonor[which(map_lgl(cumPrevFirmicutesDummyDonor, function(x) length(x) == 0))] <- NULL
cumPrevFirmicutesDummyDonor <- do.call('rbind', map2(cumPrevFirmicutesDummyDonor, names(cumPrevFirmicutesDummyDonor), function(x, nb) {
  return(data.frame(acceptor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][2]),
                    donor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][3]),
                    postFMTFirmicutesCumPrev = unlist(x),
                    dataset = nb))
}))
#cumPrevFirmicutesDummyDonor[which(map_lgl(cumPrevFirmicutesDummyDonor, function(x) length(x) == 0))] <- NULL
cumPrevFirmicutesDummyDonor <- cumPrevFirmicutesDummyDonor %>% relocate(dataset, .after = donor_fmt_id)   
rownames(cumPrevFirmicutesDummyDonor) <- NULL                                      
cumPrevFirmicutesDummyDonor <- as_tibble(cumPrevFirmicutesDummyDonor)
cumPrevFirmicutesDummyDonor <- inner_join(cumPrevFirmicutesDummyDonor %>% 
                                  mutate(predictedPostFMTFirmicutesCumPrev = postFMTFirmicutesCumPrev) %>%
                                  mutate(acceptor_fmt_id = as.double(as.vector(acceptor_fmt_id))) %>%
                                  select(acceptor_fmt_id,
                                         donor_fmt_id, 
                                         dataset, 
                                         predictedPostFMTFirmicutesCumPrev), 
                                richnessPredictedUnalteredInstance %>% 
                                  mutate(predictedCumPrevalenceFirmicutesUnalteredInstance = cumPrevalenceFirmicutes) %>%
                                  select(fmt_id, predictedCumPrevalenceFirmicutesUnalteredInstance),
                                by = c('acceptor_fmt_id' = 'fmt_id'))
cumPrevFirmicutesDummyDonor <- cumPrevFirmicutesDummyDonor %>% 
  mutate(differenceFirmicutesPrevalence = predictedPostFMTFirmicutesCumPrev - predictedCumPrevalenceFirmicutesUnalteredInstance)
cumPrevFirmicutesDummyDonor$donor_fmt_id <- factor(as.vector(cumPrevFirmicutesDummyDonor$donor_fmt_id), levels = cumPrevFirmicutesDummyDonor %>% 
  group_by(donor_fmt_id) %>% 
  summarize(differenceFirmicutesPrevalence = mean(differenceFirmicutesPrevalence)) %>% 
  arrange(desc(differenceFirmicutesPrevalence)) %>%
  pull(donor_fmt_id))      





donorOrder <- cumPrevFirmicutesDummyDonor %>% 
group_by(donor_fmt_id) %>%
summarize(med = median(differenceFirmicutesPrevalence)) %>%
arrange(desc(med)) %>%
pull(donor_fmt_id)
pp <- ggplot(cumPrevFirmicutesDummyDonor %>%
             mutate(donor_fmt_id = factor(as.vector(donor_fmt_id),
                                         levels = donorOrder)), aes(x = donor_fmt_id, y = differenceFirmicutesPrevalence)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ylab("Difference in\npredicted post-FMT\nFirmicute species prevalence upon donor exchange") +
  xlab("Donor individual")
pp3 <- ggplot(realValues %>% 
              filter(fmt_id %in% cumPrevFirmicutesDummyDonor$donor_fmt_id) %>%
              mutate(fmt_id = factor(fmt_id, levels = levels(cumPrevFirmicutesDummyDonor$donor_fmt_id))), 
              aes(x = fmt_id, 
                  y = donorFirmicutePrevalence)) + 
  geom_point(size = 1) + 
  geom_smooth(aes(group = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("prevalence of\nFirmicute species bacteria in implanted donor")+
  ylim(c(0, max(realValues$donorFirmicutePrevalence) * 1.05))
options(repr.plot.width=14, repr.plot.height=8)    
# This figure isn't part of the final manuscript
ggsave(plot = pp, filename = "figures/BestDonorsPlot1_FirmicutesPrev.pdf", width = 4.25, height = 2.1)





########################################
### cumulative prevalence Bacteroidetes ###
########################################
cumPrevBacteroidetesDummyDonor[which(map_lgl(cumPrevBacteroidetesDummyDonor, function(x) length(x) == 0))] <- NULL
cumPrevBacteroidetesDummyDonor <- do.call('rbind', map2(cumPrevBacteroidetesDummyDonor, names(cumPrevBacteroidetesDummyDonor), function(x, nb) {
  return(data.frame(acceptor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][2]),
                    donor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][3]),
                    postFMTBacteroidetesCumPrev = unlist(x),
                    dataset = nb))
}))
#cumPrevBacteroidetesDummyDonor[which(map_lgl(cumPrevBacteroidetesDummyDonor, function(x) length(x) == 0))] <- NULL
cumPrevBacteroidetesDummyDonor <- cumPrevBacteroidetesDummyDonor %>% relocate(dataset, .after = donor_fmt_id)   
rownames(cumPrevBacteroidetesDummyDonor) <- NULL                                      
cumPrevBacteroidetesDummyDonor <- as_tibble(cumPrevBacteroidetesDummyDonor)
cumPrevBacteroidetesDummyDonor <- inner_join(cumPrevBacteroidetesDummyDonor %>% 
                                  mutate(predictedPostFMTBacteroidetesCumPrev = postFMTBacteroidetesCumPrev) %>%
                                  mutate(acceptor_fmt_id = as.double(as.vector(acceptor_fmt_id))) %>%
                                  select(acceptor_fmt_id,
                                         donor_fmt_id, 
                                         dataset, 
                                         predictedPostFMTBacteroidetesCumPrev), 
                                richnessPredictedUnalteredInstance %>% 
                                  mutate(predictedCumPrevalenceBacteroidetesUnalteredInstance = cumPrevalenceBacteroidetes) %>%
                                  select(fmt_id, predictedCumPrevalenceBacteroidetesUnalteredInstance),
                                by = c('acceptor_fmt_id' = 'fmt_id'))
cumPrevBacteroidetesDummyDonor <- cumPrevBacteroidetesDummyDonor %>% 
  mutate(differenceBacteroidetesPrevalence = predictedPostFMTBacteroidetesCumPrev - predictedCumPrevalenceBacteroidetesUnalteredInstance)
cumPrevBacteroidetesDummyDonor$donor_fmt_id <- factor(as.vector(cumPrevBacteroidetesDummyDonor$donor_fmt_id), levels = cumPrevBacteroidetesDummyDonor %>% 
  group_by(donor_fmt_id) %>% 
  summarize(differenceBacteroidetesPrevalence = mean(differenceBacteroidetesPrevalence)) %>% 
  arrange(desc(differenceBacteroidetesPrevalence)) %>%
  pull(donor_fmt_id))      





donorOrder <- cumPrevBacteroidetesDummyDonor %>% 
group_by(donor_fmt_id) %>%
summarize(med = median(differenceBacteroidetesPrevalence)) %>%
arrange(desc(med)) %>%
pull(donor_fmt_id)
pp <- ggplot(cumPrevBacteroidetesDummyDonor %>%
             mutate(donor_fmt_id = factor(as.vector(donor_fmt_id),
                                         levels = donorOrder)), aes(x = donor_fmt_id, y = differenceBacteroidetesPrevalence)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ylab("Difference in\npredicted post-FMT\nBacteroidete species prevalence upon donor exchange") +
  xlab("Donor individual")
pp3 <- ggplot(realValues %>% 
              filter(fmt_id %in% cumPrevBacteroidetesDummyDonor$donor_fmt_id) %>%
              mutate(fmt_id = factor(fmt_id, levels = levels(cumPrevBacteroidetesDummyDonor$donor_fmt_id))), 
              aes(x = fmt_id, 
                  y = donorBacteroidetesPrevalence)) + 
  geom_point(size = 1) + 
  geom_smooth(aes(group = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("prevalence of\nBacteroidete species bacteria in implanted donor")+
  ylim(c(0, max(realValues$donorBacteroidetesAbundance) * 1.05))
options(repr.plot.width=14, repr.plot.height=8) 
# This figure isn't part of the final manuscript   
ggsave(plot = pp, filename = "figures/BestDonorsPlot1_BacteroidetesPrev.pdf", width = 4.25, height = 2.1)

############################################
### cumulative abundance Proteo bacteria ###
############################################
cumAbundanceProteoDummyDonor[which(map_lgl(cumAbundanceProteoDummyDonor, function(x) length(x) == 0))] <- NULL
cumAbundanceProteoDummyDonor <- do.call('rbind', map2(cumAbundanceProteoDummyDonor, names(cumAbundanceProteoDummyDonor), function(x, nb) {
  return(data.frame(acceptor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][2]),
                    donor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][3]),
                    postFMTProteoCumAb = unlist(x),
                    dataset = nb))
}))
cumAbundanceProteoDummyDonor[which(map_lgl(cumAbundanceProteoDummyDonor, function(x) length(x) == 0))] <- NULL
cumAbundanceProteoDummyDonor <- cumAbundanceProteoDummyDonor %>% relocate(postFMTProteoCumAb, .after = donor_fmt_id)   
rownames(cumAbundanceProteoDummyDonor) <- NULL                                      
cumAbundanceProteoDummyDonor <- as_tibble(cumAbundanceProteoDummyDonor)
cumAbundanceProteoDummyDonor <- inner_join(cumAbundanceProteoDummyDonor %>% 
                                  mutate(predictedPostFMTProteoCumAb = postFMTProteoCumAb) %>%
                                  mutate(acceptor_fmt_id = as.double(as.vector(acceptor_fmt_id))) %>%
                                  select(acceptor_fmt_id,
                                         donor_fmt_id, 
                                         dataset, 
                                         predictedPostFMTProteoCumAb), 
                                abundancePredictedUnalteredInstance %>% 
                                  mutate(predictedCumAbundanceProteobacteria = cumAbundanceProteobacteria) %>%
                                  select(fmt_id, predictedCumAbundanceProteobacteria),
                                by = c('acceptor_fmt_id' = 'fmt_id'))
cumAbundanceProteoDummyDonor <- cumAbundanceProteoDummyDonor %>% 
  mutate(differenceProteoAbundance = predictedPostFMTProteoCumAb - predictedCumAbundanceProteobacteria)
cumAbundanceProteoDummyDonor$donor_fmt_id <- factor(as.vector(cumAbundanceProteoDummyDonor$donor_fmt_id), levels = cumAbundanceProteoDummyDonor %>% 
  group_by(donor_fmt_id) %>% 
  summarize(differenceProteoAbundance = mean(differenceProteoAbundance)) %>% 
  arrange(desc(differenceProteoAbundance)) %>%
  pull(donor_fmt_id)) 
                                           





donorOrder <- cumAbundanceProteoDummyDonor %>% 
group_by(donor_fmt_id) %>%
summarize(med = median(differenceProteoAbundance)) %>%
arrange(desc(med)) %>%
pull(donor_fmt_id)
pp <- ggplot(cumAbundanceProteoDummyDonor %>%
             mutate(donor_fmt_id = factor(as.vector(donor_fmt_id),
                                         levels = donorOrder)), aes(x = donor_fmt_id, y = differenceProteoAbundance)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ylab("Difference in\npredicted post-FMT\ proteobacterial abundance upon donor exchange") +
  xlab("Donor individual")
pp3 <- ggplot(realValues %>% 
              filter(fmt_id %in% cumAbundanceProteoDummyDonor$donor_fmt_id) %>%
              mutate(fmt_id = factor(fmt_id, levels = levels(cumAbundanceProteoDummyDonor$donor_fmt_id))), 
              aes(x = fmt_id, 
                  y = donorProteoAbundance)) + 
  geom_point(size = 1) + 
  geom_smooth(aes(group = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("abundance of\nProteobacteria in implanted donor")+
  ylim(c(0, max(realValues$donorProteoAbundance) * 1.05))
options(repr.plot.width=14, repr.plot.height=8)    
# This figure isn't part of the final manuscript
ggsave(plot = pp, filename = "figures/BestDonorsPlot1_ProteoAbundance.pdf", width = 4.25, height = 2.1)

#######################################
### cumulative abundance Firmicutes ###
#######################################
cumAbundanceFirmicutesDummyDonor[which(map_lgl(cumAbundanceFirmicutesDummyDonor, function(x) length(x) == 0))] <- NULL
cumAbundanceFirmicutesDummyDonor <- do.call('rbind', map2(cumAbundanceFirmicutesDummyDonor, names(cumAbundanceFirmicutesDummyDonor), function(x, nb) {
  return(data.frame(acceptor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][2]),
                    donor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][3]),
                    postFMTFirmicutesCumAb = unlist(x),
                    dataset = nb))
}))
cumAbundanceFirmicutesDummyDonor[which(map_lgl(cumAbundanceFirmicutesDummyDonor, function(x) length(x) == 0))] <- NULL
cumAbundanceFirmicutesDummyDonor <- cumAbundanceFirmicutesDummyDonor %>% relocate(postFMTFirmicutesCumAb, .after = donor_fmt_id)   
rownames(cumAbundanceFirmicutesDummyDonor) <- NULL                                      
cumAbundanceFirmicutesDummyDonor <- as_tibble(cumAbundanceFirmicutesDummyDonor)
cumAbundanceFirmicutesDummyDonor <- inner_join(cumAbundanceFirmicutesDummyDonor %>% 
                                  mutate(predictedPostFMTFirmicutesCumAb = postFMTFirmicutesCumAb) %>%
                                  mutate(acceptor_fmt_id = as.double(as.vector(acceptor_fmt_id))) %>%
                                  select(acceptor_fmt_id,
                                         donor_fmt_id, 
                                         dataset, 
                                         predictedPostFMTFirmicutesCumAb), 
                                abundancePredictedUnalteredInstance %>% 
                                  mutate(predictedCumAbundanceFirmicutes = cumAbundanceFirmicutes) %>%
                                  select(fmt_id, predictedCumAbundanceFirmicutes),
                                by = c('acceptor_fmt_id' = 'fmt_id'))
cumAbundanceFirmicutesDummyDonor <- cumAbundanceFirmicutesDummyDonor %>% 
  mutate(differenceFirmicutesAbundance = predictedPostFMTFirmicutesCumAb - predictedCumAbundanceFirmicutes)
cumAbundanceFirmicutesDummyDonor$donor_fmt_id <- factor(as.vector(cumAbundanceFirmicutesDummyDonor$donor_fmt_id), levels = cumAbundanceFirmicutesDummyDonor %>% 
  group_by(donor_fmt_id) %>% 
  summarize(differenceFirmicutesAbundance = mean(differenceFirmicutesAbundance)) %>% 
  arrange(desc(differenceFirmicutesAbundance)) %>%
  pull(donor_fmt_id))               
                                               





donorOrder <- cumAbundanceFirmicutesDummyDonor %>% 
group_by(donor_fmt_id) %>%
summarize(med = median(differenceFirmicutesAbundance)) %>%
arrange(desc(med)) %>%
pull(donor_fmt_id)
pp <- ggplot(cumAbundanceFirmicutesDummyDonor %>%
             mutate(donor_fmt_id = factor(as.vector(donor_fmt_id),
                                         levels = donorOrder)), aes(x = donor_fmt_id, y = differenceFirmicutesAbundance)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ylab("Difference in\npredicted post-FMT\ Firmicutes abundance upon donor exchange") +
  xlab("Donor individual")
pp3 <- ggplot(realValues %>% 
              filter(fmt_id %in% cumAbundanceFirmicutesDummyDonor$donor_fmt_id) %>%
              mutate(fmt_id = factor(fmt_id, levels = levels(cumAbundanceFirmicutesDummyDonor$donor_fmt_id))), 
              aes(x = fmt_id, 
                  y = donorFirmicuteAbundance)) + 
  geom_point(size = 1) + 
  geom_smooth(aes(group = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("abundance of\nFirmicutes in implanted donor")+
  ylim(c(0, max(realValues$donorFirmicuteAbundance) * 1.05))
options(repr.plot.width=14, repr.plot.height=8)    
# This figure isn't part of the final manuscript
ggsave(plot = pp, filename = "figures/BestDonorsPlot1_FirmicuteAbundance.pdf", width = 4.25, height = 2.1)





#######################################
### cumulative abundance Bacteroidetes ###
#######################################
cumAbundanceBacteroidetesDummyDonor[which(map_lgl(cumAbundanceBacteroidetesDummyDonor, function(x) length(x) == 0))] <- NULL
cumAbundanceBacteroidetesDummyDonor <- do.call('rbind', map2(cumAbundanceBacteroidetesDummyDonor, names(cumAbundanceBacteroidetesDummyDonor), function(x, nb) {
  return(data.frame(acceptor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][2]),
                    donor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][3]),
                    postFMTBacteroidetesCumAb = unlist(x),
                    dataset = nb))
}))
cumAbundanceBacteroidetesDummyDonor[which(map_lgl(cumAbundanceBacteroidetesDummyDonor, function(x) length(x) == 0))] <- NULL
cumAbundanceBacteroidetesDummyDonor <- cumAbundanceBacteroidetesDummyDonor %>% relocate(postFMTBacteroidetesCumAb, .after = donor_fmt_id)   
rownames(cumAbundanceBacteroidetesDummyDonor) <- NULL                                      
cumAbundanceBacteroidetesDummyDonor <- as_tibble(cumAbundanceBacteroidetesDummyDonor)
cumAbundanceBacteroidetesDummyDonor <- inner_join(cumAbundanceBacteroidetesDummyDonor %>% 
                                                 mutate(predictedPostFMTBacteroidetesCumAb = postFMTBacteroidetesCumAb) %>%
                                                 mutate(acceptor_fmt_id = as.double(as.vector(acceptor_fmt_id))) %>%
                                                 select(acceptor_fmt_id,
                                                        donor_fmt_id, 
                                                        dataset, 
                                                        predictedPostFMTBacteroidetesCumAb), 
                                               abundancePredictedUnalteredInstance %>% 
                                                 mutate(predictedCumAbundanceBacteroidetes = cumAbundanceBacteroidetes) %>%
                                                 select(fmt_id, predictedCumAbundanceBacteroidetes),
                                               by = c('acceptor_fmt_id' = 'fmt_id'))
cumAbundanceBacteroidetesDummyDonor <- cumAbundanceBacteroidetesDummyDonor %>% 
  mutate(differenceBacteroidetesAbundance = predictedPostFMTBacteroidetesCumAb - predictedCumAbundanceBacteroidetes)
cumAbundanceBacteroidetesDummyDonor$donor_fmt_id <- factor(as.vector(cumAbundanceBacteroidetesDummyDonor$donor_fmt_id), levels = cumAbundanceBacteroidetesDummyDonor %>% 
                                                          group_by(donor_fmt_id) %>% 
                                                          summarize(differenceBacteroidetesAbundance = mean(differenceBacteroidetesAbundance)) %>% 
                                                          arrange(desc(differenceBacteroidetesAbundance)) %>%
                                                          pull(donor_fmt_id))  
                                                  





donorOrder <- cumAbundanceBacteroidetesDummyDonor %>% 
group_by(donor_fmt_id) %>%
summarize(med = median(differenceBacteroidetesAbundance)) %>%
arrange(desc(med)) %>%
pull(donor_fmt_id)
pp <- ggplot(cumAbundanceBacteroidetesDummyDonor  %>%
             mutate(donor_fmt_id = factor(as.vector(donor_fmt_id),
                                         levels = donorOrder)), aes(x = donor_fmt_id, y = differenceBacteroidetesAbundance)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ylab("Difference in\npredicted post-FMT\ Bacteroidetes abundance upon donor exchange") +
  xlab("Donor individual")
pp3 <- ggplot(realValues %>% 
                filter(fmt_id %in% cumAbundanceBacteroidetesDummyDonor$donor_fmt_id) %>%
                mutate(fmt_id = factor(fmt_id, levels = levels(cumAbundanceBacteroidetesDummyDonor$donor_fmt_id))), 
              aes(x = fmt_id, 
                  y = donorBacteroidetesAbundance)) + 
  geom_point(size = 1) + 
  geom_smooth(aes(group = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("abundance of\nBacteroidetes in implanted donor")+
  ylim(c(0, max(realValues$donorBacteroidetesAbundance) * 1.05))
options(repr.plot.width=14, repr.plot.height=8) 
# This figure isn't part of the final manuscript   
ggsave(plot = pp, filename = "figures/BestDonorsPlot1_BacteroidetesAbundance.pdf", width = 4.25, height = 2.1)





############################################
### cumulative abundance Oral bacteria ###
############################################
cumAbundanceOralDummyDonor[which(map_lgl(cumAbundanceOralDummyDonor, function(x) length(x) == 0))] <- NULL
cumAbundanceOralDummyDonor <- do.call('rbind', map2(cumAbundanceOralDummyDonor, names(cumAbundanceOralDummyDonor), function(x, nb) {
  return(data.frame(acceptor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][2]),
                    donor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][3]),
                    postFMTOralCumAb = unlist(x),
                    dataset = nb))
}))
cumAbundanceOralDummyDonor[which(map_lgl(cumAbundanceOralDummyDonor, function(x) length(x) == 0))] <- NULL
cumAbundanceOralDummyDonor <- cumAbundanceOralDummyDonor %>% relocate(postFMTOralCumAb, .after = donor_fmt_id)   
rownames(cumAbundanceOralDummyDonor) <- NULL                                      
cumAbundanceOralDummyDonor <- as_tibble(cumAbundanceOralDummyDonor)
cumAbundanceOralDummyDonor <- inner_join(cumAbundanceOralDummyDonor %>% 
                                  mutate(predictedPostFMTOralCumAb = postFMTOralCumAb) %>%
                                  mutate(acceptor_fmt_id = as.double(as.vector(acceptor_fmt_id))) %>%
                                  select(acceptor_fmt_id,
                                         donor_fmt_id, 
                                         dataset, 
                                         predictedPostFMTOralCumAb), 
                                abundancePredictedUnalteredInstance %>% 
                                  mutate(predictedCumAbundanceOral = cumAbundanceOral) %>%
                                  select(fmt_id, predictedCumAbundanceOral),
                                by = c('acceptor_fmt_id' = 'fmt_id'))
cumAbundanceOralDummyDonor <- cumAbundanceOralDummyDonor %>% 
  mutate(differenceOralAbundance = predictedPostFMTOralCumAb - predictedCumAbundanceOral)
cumAbundanceOralDummyDonor$donor_fmt_id <- factor(as.vector(cumAbundanceOralDummyDonor$donor_fmt_id), 
                                                  levels = cumAbundanceOralDummyDonor %>% 
  group_by(donor_fmt_id) %>% 
  summarize(differenceOralAbundance = mean(differenceOralAbundance)) %>% 
  arrange(desc(differenceOralAbundance)) %>%
  pull(donor_fmt_id))               
pp <- ggplot(cumAbundanceOralDummyDonor, aes(x = donor_fmt_id, y = differenceOralAbundance)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ylab("Difference in\npredicted post-FMT\nOral bacterial abundance upon donor exchange") +
  xlab("Donor individual")
pp3 <- ggplot(realValues %>% 
              filter(fmt_id %in% cumAbundanceOralDummyDonor$donor_fmt_id) %>%
              mutate(fmt_id = factor(fmt_id, levels = levels(cumAbundanceOralDummyDonor$donor_fmt_id))), 
              aes(x = fmt_id, 
                  y = donorOralAbundance)) + 
  geom_point(size = 1) + 
  geom_smooth(aes(group = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("abundance of\nOral bacteria in implanted donor")+
  ylim(c(0, max(realValues$oralAbundance) * 1.05))
options(repr.plot.width=14, repr.plot.height=8)    
# This figure isn't part of the final manuscript
ggsave(plot = pp, filename = "figures/BestDonorsPlot1_OralAbundance.pdf", width = 4.25, height = 2.1)





############################################
### cumulative abundance Predict1 bacteria ###
############################################
cumAbundanceTopPredict1DummyDonor[which(map_lgl(cumAbundanceTopPredict1DummyDonor, function(x) length(x) == 0))] <- NULL
cumAbundanceTopPredict1DummyDonor <- do.call('rbind', map2(cumAbundanceTopPredict1DummyDonor, names(cumAbundanceTopPredict1DummyDonor), function(x, nb) {
  return(data.frame(acceptor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][2]),
                    donor_fmt_id = map_chr(names(x), function(z) str_split(z, "___")[[1]][3]),
                    postFMTPredict1CumAb = unlist(x),
                    dataset = nb))
}))
cumAbundanceTopPredict1DummyDonor[which(map_lgl(cumAbundanceTopPredict1DummyDonor, function(x) length(x) == 0))] <- NULL
cumAbundanceTopPredict1DummyDonor <- cumAbundanceTopPredict1DummyDonor %>% relocate(postFMTPredict1CumAb, .after = donor_fmt_id)   
rownames(cumAbundanceTopPredict1DummyDonor) <- NULL                                      
cumAbundanceTopPredict1DummyDonor <- as_tibble(cumAbundanceTopPredict1DummyDonor)
cumAbundanceTopPredict1DummyDonor <- inner_join(cumAbundanceTopPredict1DummyDonor %>% 
                                  mutate(predictedPostFMTPredict1CumAb = postFMTPredict1CumAb) %>%
                                  mutate(acceptor_fmt_id = as.double(as.vector(acceptor_fmt_id))) %>%
                                  select(acceptor_fmt_id,
                                         donor_fmt_id, 
                                         dataset, 
                                         predictedPostFMTPredict1CumAb), 
                                abundancePredictedUnalteredInstance %>% 
                                  mutate(predictedCumAbundancePredict1 = cumAbundancePredict1) %>%
                                  select(fmt_id, predictedCumAbundancePredict1),
                                by = c('acceptor_fmt_id' = 'fmt_id'))
cumAbundanceTopPredict1DummyDonor <- cumAbundanceTopPredict1DummyDonor %>% 
  mutate(differencePredict1Abundance = predictedPostFMTPredict1CumAb - predictedCumAbundancePredict1)
cumAbundanceTopPredict1DummyDonor$donor_fmt_id <- factor(as.vector(cumAbundanceTopPredict1DummyDonor$donor_fmt_id), 
                                                  levels = cumAbundanceTopPredict1DummyDonor %>% 
  group_by(donor_fmt_id) %>% 
  summarize(differencePredict1Abundance = mean(differencePredict1Abundance)) %>% 
  arrange(desc(differencePredict1Abundance)) %>%
  pull(donor_fmt_id))               
pp <- ggplot(cumAbundanceTopPredict1DummyDonor, aes(x = donor_fmt_id, y = differencePredict1Abundance)) + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ylab("Difference in\npredicted post-FMT\nPredict1 abundance upon donor exchange") +
  xlab("Donor individual")
pp3 <- ggplot(realValues %>% 
              filter(fmt_id %in% cumAbundanceTopPredict1DummyDonor$donor_fmt_id) %>%
              mutate(fmt_id = factor(fmt_id, levels = levels(cumAbundanceTopPredict1DummyDonor$donor_fmt_id))), 
              aes(x = fmt_id, 
                  y = predict1Abundance)) + 
  geom_point(size = 1) + 
  geom_smooth(aes(group = 1)) +
  theme_embl() + 
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  ylab("abundance of\nPredict1 bacteria in implanted donor")+
  ylim(c(0, max(realValues$predict1Abundance) * 1.05))
options(repr.plot.width=14, repr.plot.height=8)    
# This figure isn't part of the final manuscript
ggsave(plot = pp, filename = "figures/BestDonorsPlot1_Predict1Abundance.pdf", width = 4.25, height = 2.1)

richnessEvalPlots <- list()
abundanceEvalPlots <- list()
options(warn = 0)

plotEval <- function(o, var1, var2, var3, x1, y1, xLab1, yLab1, xLab2, yLab2, lims, logScale = NULL) {
    dat <- rbind(realValues %>% 
               ungroup() %>%
               mutate(type = "real",
                      acceptor_fmt_id = fmt_id) %>% 
               select(all_of(c("acceptor_fmt_id", var1, "type"))), 
               #select(acceptor_fmt_id, postFMTRichness, type ),
             o %>% 
               ungroup() %>% 
               mutate(!!var1 := .data[[var2]],
                      acceptor_fmt_id = fmt_id) %>%
               #select(acceptor_fmt_id, postFMTRichness) %>%
               select(all_of(c("acceptor_fmt_id", var1))) %>%
               mutate(type = "predicted")) %>% 
  pivot_wider(id_cols = acceptor_fmt_id, names_from = type, values_from = .data[[var1]])
if (!is.null(logScale)) {
    p <- ggplot(dat, aes(x = predicted+logScale, y = real+logScale))
} else {
    p <- ggplot(dat, aes(x = predicted, y = real))
}
p <- p + 
  theme_embl() +
  geom_point(size = 1, alpha = 0.5) + 
  #xlab("Richness post-FMT (Predicted)") +
  #ylab("Richness post-FMT (Real)") +
  xlab(xLab1) +
  ylab(yLab1) +
  xlim(lims) +
  ylim(lims) +
  annotate(geom = 'text', x = x1, y = y1, label = str_c("Spearman cor.: ", round(cor(dat$real, dat$predicted, method = 'spearman'), 2)), size = 4.25) +
  NULL
if (!is.null(logScale)) {
    p <- p +
    scale_x_log10(limits = lims + logScale) +
    scale_y_log10(limits = lims + logScale)
#     xlim(lims + logScale) + 
#     ylim(lims + logScale)
}    
    
#ggsave("figures/Richness_real_vs_pred_scatter.pdf", plot = p, width = 3.25, height = 3.25)
if (!is.null(logScale)) {
    p2 <- ggplot(realValues, aes(x = .data[[var3]]+logScale, y = .data[[var1]]+logScale)) + 
    geom_point(size = 1, alpha = 0.5)
} else {
    p2 <- ggplot(realValues, aes(x = .data[[var3]], y = .data[[var1]])) + geom_point(size = 1, alpha = 0.5)
}
p2 <- p2 + 
theme_embl() +
#xlab("Richness (Donor)") +
#ylab("Richness (post-FMT)") + 
  xlab(xLab2) +
  ylab(yLab2) +
  xlim(lims) +
  ylim(lims) +
annotate(geom = 'text', x = x1, y = y1, label = str_c("Spearman cor: ", round(cor(realValues %>% pull(.data[[var3]]), 
                                                                                  realValues %>% pull(.data[[var1]]), 
                                                                                                      method = 'spearman'), 2)), size = 4.25)
if (!is.null(logScale)) {
    p2 <- p2 + 
    scale_x_log10(limits = lims + logScale) +
    scale_y_log10(limits = lims + logScale) 
    #xlim(lims + logScale) + 
    #ylim(lims + logScale)
}    
#ggsave("figures/Total_richness_eval.pdf", plot = p+p2, width = 6.5, height = 3.25)
#ggsave("figures/Total_richness_eval.pdf", plot = p+p2, width = 6.5, height = 3.25)
#options(repr.plot.width=6.5, repr.plot.height=3.25)
return(list(p, p2))
}

l <- plotEval(o = richnessPredictedUnalteredInstance, 
              var1 = "postFMTRichness",
        var2 = 'richness',
        var3 = "DonorRichness",
        x1 = 210,
        y1 = 10, 
        xLab1 = "Richness post-FMT (predicted)",
        yLab1 = "Richness post-FMT (real)",
        xLab2 = "Richness (donor)",
        yLab2 = "Richness (post-FMT)",
        lims = c(0, 450))
options(repr.plot.width=6.5, repr.plot.height=3.25)
richnessEvalPlots[[length(richnessEvalPlots)  + 1]] <- l[[1]]
richnessEvalPlots[[length(richnessEvalPlots)  + 1]] <- l[[2]]

l <- plotEval(o = richnessPredictedUnalteredInstance, 
              var1 = "proteoPrevalence",
        var2 = 'cumPrevalenceProteobacteria',
        var3 = "donorProteoPrevalence",
        x1 = 12,
        y1 = 0.75, 
        xLab1 = "Proteobacterial\nrichness post-FMT (redicted)",
        yLab1 = "Proteobacterial\nrichness post-FMT (real)",
        xLab2 = "Proteobacterial\nrichness (donor)",
        yLab2 = "Proteobacterial\nrichness (post-FMT)",
        lims = c(0, 20))
richnessEvalPlots[[length(richnessEvalPlots)  + 1]] <- l[[1]]
richnessEvalPlots[[length(richnessEvalPlots)  + 1]] <- l[[2]]

l <- plotEval(o = richnessPredictedUnalteredInstance, 
              var1 = "oralPrevalence",
        var2 = 'cumPrevalenceOral',
        var3 = "donorOralPrevalence",
        x1 = 23,
        y1 = 1, 
        xLab1 = "Oral bacteria\nrichness post-FMT (predicted)",
        yLab1 = "Oral bacteria\nrichness post-FMT (real)",
        xLab2 = "Oral bacteria\nrichness (donor)",
        yLab2 = "Oral bacteria\nrichness (post-FMT)",
        lims = c(0, 40))
richnessEvalPlots[[length(richnessEvalPlots)  + 1]] <- l[[1]]
richnessEvalPlots[[length(richnessEvalPlots)  + 1]] <- l[[2]]

l <- plotEval(o = richnessPredictedUnalteredInstance, 
              var1 = "predict1Prevalence",
var2 ='cumPrevalencePredict1',
var3 ='donorPredict1Prevalence',
x1 = 20,
y1 = 70,
xLab1 = "Top Predict 1 bacteria\nrichness post-FMT (predicted)",
yLab1 = "Top Predict 1 bacteria\nrichness post-FMT (real)",
xLab2 = "Top Predict 1 bacteria\nrichness (donor)",
yLab2 = "Top Predict 1 bacteria\nrichness (post-FMT)",
lims = c(0,70))
richnessEvalPlots[[length(richnessEvalPlots)  + 1]] <- l[[1]]
richnessEvalPlots[[length(richnessEvalPlots)  + 1]] <- l[[2]]

l <- plotEval(o = richnessPredictedUnalteredInstance, 
              var1 = "FirmicutePrevalence",
var2 ='cumPrevalenceFirmicutes',
var3 ='donorFirmicutePrevalence',
x1 = 250,
y1 = 5,
xLab1 = "Firmicutes bacteria\nrichness post-FMT (predicted)",
yLab1 = "Firmicutes bacteria\nrichness post-FMT (real)",
xLab2 = "Firmicutes bacteria\nrichness (donor)",
yLab2 = "Firmicutes bacteria\nrichness (post-FMT)",
lims = c(0,350))
richnessEvalPlots[[length(richnessEvalPlots)  + 1]] <- l[[1]]
richnessEvalPlots[[length(richnessEvalPlots)  + 1]] <- l[[2]]

l <- plotEval(o = richnessPredictedUnalteredInstance, 
              var1 = "BacteroidetesPrevalence",
var2 ='cumPrevalenceBacteroidetes',
var3 ='donorBacteroidetesPrevalence',
x1 = 40,
y1 = 5,
xLab1 = "Bacteroidetes bacteria\nrichness post-FMT (predicted)",
yLab1 = "Bacteroidetes bacteria\nrichness post-FMT(real)",
xLab2 = "Bacteroidetes bacteria\nrichness (donor)",
yLab2 = "Bacteroidetes bacteria\nrichness (post-FMT)",
lims = c(0,60))
richnessEvalPlots[[length(richnessEvalPlots)  + 1]] <- l[[1]]
richnessEvalPlots[[length(richnessEvalPlots)  + 1]] <- l[[2]]

test <- tibble(a = richnessEvalPlots) %>%
    mutate(group = unlist(map(c("Total richness",
                             "Proteobacterial richness",
                             "Oral bacterial richness",
                             "Predict 1 richness",
                             "Firmicutes richness",
                             "Bacteroidetes richness"), function(x) rep(x, 2)))) %>%
    group_by(group) %>%
                              nest() %>%
                              mutate(plots = map2(data, group, function(x, b) {
                                  #return((x$a[[1]] + x$a[[2]]) + plot_annotation(
                                  return((x$a[[1]] + 
                                         ggtitle(b) + 
                                         theme(plot.title = element_text(hjust = 1.75))) +
                                         x$a[[2]])
    #title = b,
    #caption = 'made with patchwork',
    #theme = theme(plot.title = element_text(size = 16,
     #                                      hjust = 0.5))
                              }))
wrap_plots(test$plots, nrow = 3)                              
ggsave(plot = wrap_plots(test$plots, nrow = 3), 
       filename = "figures/ED_Figure_9.pdf", 
       width = 16, height = 12)

l <- plotEval(o = abundancePredictedUnalteredInstance,
              var1 = "proteoAbundance",
var2 ='cumAbundanceProteobacteria',
var3 ='donorProteoAbundance',
# x1 = 26,
# y1 = 1,
x1 = 0.001,
y1 = 0.0001,                  
xLab1 = "Proteobacterial\nrelative abundance post-FMT (predicted, %)",
yLab1 = "Proteobacterial\nrelative abundance post-FMT (real, %)",
xLab2 = "Proteobacterial\nrelative abundance (donor, %)",
yLab2 = "Proteobacterial\nrelative abundance (post-FMT, %)",
logScale = 1E-5,
lims = c(0,40))
options(repr.plot.width=8, repr.plot.height=4)
abundanceEvalPlots[[length(abundanceEvalPlots)  + 1]] <- l[[1]]
abundanceEvalPlots[[length(abundanceEvalPlots)  + 1]] <- l[[2]]


l <- plotEval(o = abundancePredictedUnalteredInstance,
              var1 = "FirmicuteAbundance",
var2 ='cumAbundanceFirmicutes',
var3 ='donorFirmicuteAbundance',
# x1 = 65,
# y1 = 2,
x1 = 1,
y1 = 1, 
xLab1 = "Firmicutes\nabundance post-FMT (predicted, %)",
yLab1 = "Firmicutes\nabundance post-FMT (real, %)",
xLab2 = "Firmicutes\nabundance (donor, %)",
yLab2 = "Firmicutes\nabundance (post-FMT, %)",
logScale = 1E-5,   
lims = c(0.1,110))
abundanceEvalPlots[[length(abundanceEvalPlots)  + 1]] <- l[[1]]
abundanceEvalPlots[[length(abundanceEvalPlots)  + 1]] <- l[[2]]


l <- plotEval(o = abundancePredictedUnalteredInstance,
              var1 = "BacteroidetesAbundance",
var2 ='cumAbundanceBacteroidetes',
var3 ='donorBacteroidetesAbundance',
# x1 = 15,
# y1 = 79,
x1 = 0.001,
y1 = 0.0001,               
xLab1 = "Bacteroidetes\nabundance post-FMT (predicted, %)",
yLab1 = "Bacteroidetes\nabundance post-FMT (real, %)",
xLab2 = "Bacteroidetes\nabundance (donor, %)",
yLab2 = "Bacteroidetes\nabundance (post-FMT, %)",
logScale = 1E-5,   
lims = c(0,100))
abundanceEvalPlots[[length(abundanceEvalPlots)  + 1]] <- l[[1]]
abundanceEvalPlots[[length(abundanceEvalPlots)  + 1]] <- l[[2]]


l <- plotEval(o = abundancePredictedUnalteredInstance,
              var1 = "predict1Abundance",
var2 ='cumAbundancePredict1',
var3 ='donorPredict1Abundance',
# x1 = 30,
# y1 = 5,
x1 = 0.001,
y1 = 0.0001,                  
xLab1 = "post-FMT relative abundance of\ntop Predict 1 bacteria (predicted, %)",
yLab1 = "post-FMT relative abundance of\ntop Predict 1 bacteria (real, %)",
xLab2 = "Cumulative relative abundance of\ntop Predict 1 bacteria (donor, %)",
yLab2 = "Cumulative relative abundance of\ntop Predict 1 bacteria (post-FMT, %)",
logScale = 1E-5,              
lims = c(0,41))
abundanceEvalPlots[[length(abundanceEvalPlots)  + 1]] <- l[[1]]
abundanceEvalPlots[[length(abundanceEvalPlots)  + 1]] <- l[[2]]
ggsave(plot = l[[1]] + l[[2]],
      filename = "figures/cumAbundancePredict1.pdf", width = 7, height = 3.2)

l <- plotEval(o = abundancePredictedUnalteredInstance,
              var1 = "oralAbundance",
var2 ='cumAbundanceOral',
var3 ='donorOralAbundance',
#x1 = 5,
#y1 = 30,
x1 = 0.001,
y1 = 0.0001,              
xLab1 = "Oral Bacteria\nrelative abundance post-FMT (predicted)",
yLab1 = "Oral Bacteria\nrelative abundance post-FMT (real)",
xLab2 = "Oral Bacteria\nabundance (donor)",
yLab2 = "Oral Bacteria\nabundance (post-FMT)",
logScale = 1E-5,              
lims = c(0,30))
abundanceEvalPlots[[length(abundanceEvalPlots)  + 1]] <- l[[1]]
abundanceEvalPlots[[length(abundanceEvalPlots)  + 1]] <- l[[2]]

test <- tibble(a = abundanceEvalPlots) %>%
    mutate(group = unlist(map(c("Proteobacterial cum. abundance",
                             "Firmicutes cum. abundance",
                             "Bacteroidetes cum. abundance",
                             "Predict 1 cum. abundance",
                             "Oral bacterial cum. abundance"), function(x) rep(x, 2)))) %>%
    group_by(group) %>%
                              nest() %>%
                              mutate(plots = map2(data, group, function(x, b) {
                                  #return((x$a[[1]] + x$a[[2]]) + plot_annotation(
                                  return((x$a[[1]] + 
                                         ggtitle(b) + 
                                         theme(plot.title = element_text(hjust = 1.75))) +
                                         x$a[[2]])
    #title = b,
    #caption = 'made with patchwork',
    #theme = theme(plot.title = element_text(size = 16,
     #                                      hjust = 0.5))
                              }))
ggsave(plot = wrap_plots(test$plots, nrow = 3), 
       filename = "figures/ED_Figure_10.pdf", 
       width = 16, height = 12)

