#DAG3_bact <- read.delim("/Users/a.p.petrenko/Desktop/awk/filtering/raw_data/DAG3_ECs_sorted_bact.txt", sep="\t") 

#loading DAG3 dataset, it's alredy filtered with cazymes from cazy.org and they are bacterial
DAG3_bact <- read.delim("/Users/a.p.petrenko/research/CAZymes/output/dag3_filtered_by_ECs_and_bact.txt", sep="\t") 

ECs_bact <- read.delim("//Users/a.p.petrenko/back_up/CAZy_ecs_bacterial.tsv", sep="\t") 
library(dplyr)
ECs_bact <- ECs_bact[,2:3]

ECs_bact <- ECs_bact %>% distinct(ECs_bact$EC, .keep_all = TRUE)
ECs_bact <- ECs_bact[,1:2]

#dag3_filtered_by_ECs_and_bact <- dag3_filtered_by_ECs[dag3_filtered_by_ECs$ECs %in% EC_bact, ]
library(readr)
#loading metadata for DAG3 8229 samples, 35 variables 
md_dag3<-as.data.frame(read_tsv("/Users/a.p.petrenko/back_up/metadata/dag3_phenotypes.txt"))
#use row names (sample ID) as the first column
rownames(md_dag3) <- md_dag3$DAG3_sampleID
md_dag3 <- md_dag3[,-1]
#deleting the column represents the birth length 
sum <-  md_dag3[, -5]
#transposing the dataframe
DAG3_bact_transp <- t(DAG3_bact)
#save it as dataframe
DAG3_bact_transp <- as.data.frame(DAG3_bact_transp)
# Assign the first row as column names
colnames(DAG3_bact_transp) <- as.character(unlist(DAG3_bact_transp[1, ]))
# Remove the first row
DAG3_bact_transp <- DAG3_bact_transp[-1, ]
#rename
df <- DAG3_bact_transp
#check for duplicates
duplicate_rows <- df[duplicated(df), ]
df <- unique(df)
# Convert all values to numeric
df[] <- lapply(df, as.numeric)
#checkiing if the min value is correct (should be 0)
min_value <- min(df[3])

#shortnening sample ID in df to merge it with metadata
# Extract the new row names by removing everything after the third underline
new_row_names <- gsub("^([^_]+_[^_]+_[^_]+)_.*$", "\\1", rownames(df))
# Remove any duplicate row names
new_row_names <- make.unique(new_row_names, sep = "_")
rownames(df) <- new_row_names
# Remove duplicate rows
df <- df[!duplicated(new_row_names), ]
df2 <- df
df1 <- md_dag3
# Merge the dataframes based on row names
merged_df <- merge(df1, df2, by = "row.names")
#save sample names
DAG3_md_df_sample_names <- merged_df$Row.names
# Remove the extra row names column
merged_df$row.names <- NULL
# Set the merged row names
rownames(merged_df) <- merged_df$Row.names
# Remove the Row.names column
merged_df$Row.names <- NULL
#rename ready df as DAG3
DAG3 <- merged_df
DAG3values <- DAG3[,35:205]
DAG3md <- DAG3[,1:34]
#Calculate the count of non-zero values for each column
non_zero_counts <- colSums(DAG3values != 0)
# Create a dataframe with the column names and non-zero counts
result_df <- data.frame(Column = colnames(DAG3values), NonZeroCount = non_zero_counts)
sum <- length(DAG3md[,1])
result_df$NonZeroPercentage <- (result_df$NonZeroCount / sum) * 100
sorted_df <- result_df[order(result_df[, 3]), ]
colnames(result_df)[1] <- "Ecs"

#barplot(sorted_df$NonZeroPercentage)
#percentages <- sorted_df$NonZeroPercentage
#
library(ggplot2)
# Sort the data frame by the percentage column
sorted_df <- sorted_df[order(sorted_df$NonZeroPercentage), ]

sorted_df$Column <- factor(sorted_df$Column, levels = sorted_df$Column)

ggplot(sorted_df, aes(x = sorted_df$Column, y = sorted_df$NonZeroPercentage)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "CAZymes", y = "Percentage of population") +
  ggtitle("Percentage of Population presents specific CAZyme") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18)
  ) +
  scale_x_discrete(labels = function(x) ifelse(x == "CAZymes", "CAZymes", x)) +
geom_hline(yintercept = 50, linetype = "solid", color = "red")


#sorting by thresholds
# Define the threshold values
thresholds <- c(99.9, 99.5, 99, 98, 95, 90, 80, 70, 60, 50, 40, 30, 20, 10, 5, 0)
# Initialize a vector to store the counts
counts <- numeric(length(thresholds))
# Calculate the counts for each threshold
for (i in 1:length(thresholds)) {
  counts[i] <- sum(result_df$NonZeroPercentage > thresholds[i])
}
# Print the counts for each threshold
for (i in 1:length(thresholds)) {
  cat("Number of CAZymes in population percentage", thresholds[i], ":", counts[i], "\n")
}

# Calculate the count of values equal to 100 and greater than or equal to 50
count100  <- result_df$NonZeroPercentage == 100
count50  <- result_df$NonZeroPercentage >= 50
EC_filter_100p <- sum(count100)
EC_filter_100p <- result_df$Ecs[count100]
EC_filter_50p <- sum(count50)
EC_filter_50p <- result_df$Ecs[count50]





library(dplyr)
DAG3_values_100 <- subset(DAG3values, select = all_of(EC_filter_100p))
DAG3_values_50 <- subset(DAG3values, select = all_of(EC_filter_50p))

#compositional data problem

#calculation the pseudocount for DAG3_values_50 
non_zero_values_DAG_values_50 <- DAG3_values_50[DAG3_values_50 != 0]
min_value <- min(non_zero_values_DAG_values_50)
pseudocount <- min_value/2
#add pseudocount to all values in matrix
DAG3_values_50_plus_pseudo <- DAG3_values_50 + pseudocount


#CLR
install.packages("compositions")
library(compositions)
clr_DAG3_50 <- clr(as.matrix(DAG3_values_50_plus_pseudo))
# Apply the CLR transformation to the dataframe
clr_DAG3_50 <- clr(as.matrix(DAG3_values_50_plus_pseudo))
# Convert the transformed matrix back to a dataframe
clr_DAG3_50 <- as.data.frame(clr_DAG3_50)

#visualising the change of value distribution before presudocount ass, after and after CLR

ggplot(DAG3_values_50, aes(x = DAG3_values_50[,1])) +
  geom_histogram(fill = "steelblue", color = "black", bins = 100) + theme_bw() +
  labs(x = "RPKM", y = "number of samples", title = "Distribution of RPKM for CAZyme 2.4.4.1 before adding the pseudocount")

ggplot(DAG3_values_50_plus_pseudo, aes(x = DAG3_values_50_plus_pseudo[,1])) +
  geom_histogram(fill = "steelblue", color = "black", bins = 100) + theme_bw() +
  labs(x = "RPKM", y = "number of samples", title = "Distribution of RPKM for CAZyme 2.4.4.1 with pseudocount")

ggplot(clr_DAG3_50, aes(x = clr_DAG3_50[,1])) +
  geom_histogram(fill = "steelblue", color = "black", bins = 100) + theme_bw() +
  labs(x = "CLR (RPKM)", y = "number of samples", title = "Distribution of CLR (RPKM) after CLR for CAZyme 2.4.4.1")



ggplot(DAG3_values_50, aes(x = DAG3_values_50[,52])) +
  geom_histogram(fill = "steelblue", color = "black", bins = 100) +
  labs(x = "RPKM", y = "number of samples", title = "Distribution of RPKM for CAZyme 3.2.1.23 before adding the pseudocount")

ggplot(DAG3_values_50_plus_pseudo, aes(x = DAG3_values_50_plus_pseudo[,52])) +
  geom_histogram(fill = "steelblue", color = "black", bins = 100) +
  labs(x = "RPKM", y = "number of samples", title = "Distribution of RPKM for CAZyme 3.2.1.23 with pseudocount")

ggplot(clr_DAG3_50, aes(x = clr_DAG3_50[,52])) +
  geom_histogram(fill = "steelblue", color = "black", bins = 100) +
  labs(x = "log((RPKM+pseudocount)/GM)", y = "number of samples", title = "Distribution of RPKM after CLR for CAZyme 3.2.1.23")

ad_results <- data.frame(Column = character(),
                         A2 = double(),
                         p_value = double(),
                         stringsAsFactors = FALSE)

for (col in 1:ncol(clr_DAG3_50)) {
  result <- ad.test(clr_DAG3_50[, col])
  column_name <- colnames(clr_DAG3_50)[col]
  ad_results <- rbind(ad_results, 
                      data.frame(Column = column_name,
                                 A2 = result$statistic,
                                 p_value = result$p.value))
}



#sorting md for those who has IBS-C
IBS_C_md <- DAG3md[DAG3md$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeC == "Y", ]
IBS_C_md <- IBS_C_md[complete.cases(IBS_C_md$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeC), ]

#sorting md for those who has IBS-D
IBS_D_md <- DAG3md[DAG3md$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeD == "Y", ]
IBS_D_md <- IBS_D_md[complete.cases(IBS_D_md$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeD), ]

#sorting md for those who has IBS-M
IBS_M_md <- DAG3md[DAG3md$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeM == "Y", ]
IBS_M_md <- IBS_M_md[complete.cases(IBS_M_md$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeM), ]

#combining these sorted md with CAZymes values from clr_DAG3_50

IBS_C <- merge(IBS_C_md, clr_DAG3_50, by = "row.names")
IBS_D <- merge(IBS_D_md, clr_DAG3_50, by = "row.names")
IBS_M <- merge(IBS_M_md, clr_DAG3_50, by = "row.names")


#IBS_C vs IBS_M
df1_W <- IBS_C
df2_W <- IBS_M
result_df <- data.frame()
for (col in 37:126) {
  test_result <- wilcox.test(df1_W[, col], df2_W[, col], paired = FALSE)
  p_value <- test_result$p.value
  row <- data.frame(Column = colnames(df1_W)[col], P_Value = p_value)
  result_df <- rbind(result_df, row)
}
IBS_C_vs_M <- result_df


#IBS_C vs IBS_D
df1_W <- IBS_C
df2_W <- IBS_D
result_df <- data.frame()
for (col in 37:126) {
  test_result <- wilcox.test(df1_W[, col], df2_W[, col], paired = FALSE)
  p_value <- test_result$p.value
  row <- data.frame(Column = colnames(df1_W)[col], P_Value = p_value)
  result_df <- rbind(result_df, row)
}
IBS_C_vs_D <- result_df


#IBS_M vs IBS_D
df1_W <- IBS_M
df2_W <- IBS_D
result_df <- data.frame()
for (col in 37:126) {
  test_result <- wilcox.test(df1_W[, col], df2_W[, col], paired = FALSE)
  p_value <- test_result$p.value
  row <- data.frame(Column = colnames(df1_W)[col], P_Value = p_value)
  result_df <- rbind(result_df, row)
}
IBS_M_vs_D <- result_df

IBS_M_vs_D <- IBS_M_vs_D[order(IBS_M_vs_D$P_Value), ]
IBS_C_vs_M <- IBS_C_vs_M[order(IBS_C_vs_M$P_Value), ]
IBS_C_vs_D <- IBS_C_vs_D[order(IBS_C_vs_D$P_Value), ]


significant_IBS_M_vs_D <- IBS_M_vs_D$Column[IBS_M_vs_D$P_Value < 0.05]
significant_IBS_C_vs_D <- IBS_C_vs_D$Column[IBS_C_vs_D$P_Value < 0.05]
significant_IBS_C_vs_M <- IBS_C_vs_M$Column[IBS_C_vs_M$P_Value < 0.05]


significant_IBS_M_vs_D_ad <- IBS_M_vs_D$Column[IBS_M_vs_D$P_Value < 0.05/270]
significant_IBS_C_vs_D_ad <- IBS_C_vs_D$Column[IBS_C_vs_D$P_Value < 0.05/270]
significant_IBS_C_vs_M_ad <- IBS_C_vs_M$Column[IBS_C_vs_M$P_Value < 0.05/270]



#the same thing but with BMI

table(DAG3md$BMI_Group)

#sorting md for those who has IBS-C
BMI_U_md <- DAG3md[DAG3md$BMI_Group == "Underweight", ]
BMI_H_md <- DAG3md[DAG3md$BMI_Group == "Healthy Weight", ]
BMI_Ov_md <- DAG3md[DAG3md$BMI_Group == "Overweight", ]
BMI_Ob_md <- DAG3md[DAG3md$BMI_Group == "Obesity", ]


#combining these sorted md with CAZymes values from clr_DAG3_50

BMI_U <- merge(BMI_U_md, clr_DAG3_50, by = "row.names")
BMI_H <- merge(BMI_H_md, clr_DAG3_50, by = "row.names")
BMI_Ov <- merge(BMI_Ov_md, clr_DAG3_50, by = "row.names")
BMI_Ob <- merge(BMI_Ob_md, clr_DAG3_50, by = "row.names")


#BMI_U vs BMI_H
df1_W <- BMI_U
df2_W <- BMI_H
result_df <- data.frame()
for (col in 37:126) {
  test_result <- wilcox.test(df1_W[, col], df2_W[, col], paired = FALSE)
  p_value <- test_result$p.value
  row <- data.frame(Column = colnames(df1_W)[col], P_Value = p_value)
  result_df <- rbind(result_df, row)
}
BMI_U_vs_H <- result_df
BMI_U_vs_H <- BMI_U_vs_H[order(BMI_U_vs_H$P_Value),]

#BMI_U vs BMI_Ov
df1_W <- BMI_U
df2_W <- BMI_Ov
result_df <- data.frame()
for (col in 37:126) {
  test_result <- wilcox.test(df1_W[, col], df2_W[, col], paired = FALSE)
  p_value <- test_result$p.value
  row <- data.frame(Column = colnames(df1_W)[col], P_Value = p_value)
  result_df <- rbind(result_df, row)
}
BMI_U_vs_Ov <- result_df
BMI_U_vs_Ov <- BMI_U_vs_Ov[order(BMI_U_vs_Ov$P_Value),]

#BMI_U vs BMI_Ob
df1_W <- BMI_U
df2_W <- BMI_Ob
result_df <- data.frame()
for (col in 37:126) {
  test_result <- wilcox.test(df1_W[, col], df2_W[, col], paired = FALSE)
  p_value <- test_result$p.value
  row <- data.frame(Column = colnames(df1_W)[col], P_Value = p_value)
  result_df <- rbind(result_df, row)
}
BMI_U_vs_Ob <- result_df
BMI_U_vs_Ob <- BMI_U_vs_Ob[order(BMI_U_vs_Ob$P_Value),]


#BMI_H vs BMI_Ov
df1_W <- BMI_H
df2_W <- BMI_Ov
result_df <- data.frame()
for (col in 37:126) {
  test_result <- wilcox.test(df1_W[, col], df2_W[, col], paired = FALSE)
  p_value <- test_result$p.value
  row <- data.frame(Column = colnames(df1_W)[col], P_Value = p_value)
  result_df <- rbind(result_df, row)
}
BMI_H_vs_Ov <- result_df
BMI_H_vs_Ov <- BMI_H_vs_Ov[order(BMI_H_vs_Ov$P_Value),]


#BMI_H vs BMI_Ob
df1_W <- BMI_H
df2_W <- BMI_Ob
result_df <- data.frame()
for (col in 37:126) {
  test_result <- wilcox.test(df1_W[, col], df2_W[, col], paired = FALSE)
  p_value <- test_result$p.value
  row <- data.frame(Column = colnames(df1_W)[col], P_Value = p_value)
  result_df <- rbind(result_df, row)
}
BMI_H_vs_Ob <- result_df
BMI_H_vs_Ob <- BMI_H_vs_Ob[order(BMI_H_vs_Ob$P_Value),]

#BMI_Ov vs BMI_Ob
df1_W <- BMI_Ov
df2_W <- BMI_Ob
result_df <- data.frame()
for (col in 37:126) {
  test_result <- wilcox.test(df1_W[, col], df2_W[, col], paired = FALSE)
  p_value <- test_result$p.value
  row <- data.frame(Column = colnames(df1_W)[col], P_Value = p_value)
  result_df <- rbind(result_df, row)
}
BMI_Ov_vs_Ob <- result_df
BMI_Ov_vs_Ob <- BMI_Ov_vs_Ob[order(BMI_Ov_vs_Ob$P_Value),]

Wilc_result <- BMI_Ov_vs_Ob

significant_Wilc <- list()

significant_Wilc$BMI_Ov_vs_Ob <- Wilc_result$Column[Wilc_result$P_Value < 0.05]

Wilc_result <- BMI_H_vs_Ob
significant_Wilc$BMI_H_vs_Ob <- Wilc_result$Column[Wilc_result$P_Value < 0.05]

Wilc_result <- BMI_U_vs_Ob
significant_Wilc$BMI_U_vs_Ob <- Wilc_result$Column[Wilc_result$P_Value < 0.05]

Wilc_result <- BMI_U_vs_H
significant_Wilc$BMI_U_vs_H <- Wilc_result$Column[Wilc_result$P_Value < 0.05]

Wilc_result <- BMI_U_vs_Ov
significant_Wilc$BMI_U_vs_Ov <- Wilc_result$Column[Wilc_result$P_Value < 0.05]

Wilc_result <- BMI_H_vs_Ov
significant_Wilc$BMI_H_vs_Ov <- Wilc_result$Column[Wilc_result$P_Value < 0.05]

significant_Wilc_before_ad <- significant_Wilc
significant_Wilc <- list()

significant_Wilc$BMI_Ov_vs_Ob <- Wilc_result$Column[Wilc_result$P_Value < 0.05/540]

Wilc_result <- BMI_H_vs_Ob
significant_Wilc$BMI_Ov_vs_Ob <- Wilc_result$Column[Wilc_result$P_Value < 0.05/540]

Wilc_result <- BMI_U_vs_Ob
significant_Wilc$BMI_U_vs_Ob <- Wilc_result$Column[Wilc_result$P_Value < 0.05/540]

Wilc_result <- BMI_U_vs_H
significant_Wilc$BMI_U_vs_H <- Wilc_result$Column[Wilc_result$P_Value < 0.05/540]

Wilc_result <- BMI_U_vs_Ov
significant_Wilc$BMI_U_vs_Ov <- Wilc_result$Column[Wilc_result$P_Value < 0.05/540]

Wilc_result <- BMI_H_vs_Ov
significant_Wilc$BMI_H_vs_Ov <- Wilc_result$Column[Wilc_result$P_Value < 0.05/540]


significant_Wilc_after_ad <- significant_Wilc
View(significant_Wilc)

BMI_U_vs_Ob
BMI_U_vs_H
BMI_U_vs_Ov
BMI_H_vs_Ov
BMI_H_vs_Ob
BMI_Ov_vs_Ob

# List of dataframes
dataframes <- list(BMI_U_vs_Ob, BMI_U_vs_H, BMI_U_vs_Ov, BMI_H_vs_Ov, BMI_H_vs_Ob, BMI_Ov_vs_Ob)
View(dataframes)
# Loop through each dataframe
for (i in seq_along(dataframes)) {
  dataframes[[i]] <- dataframes[[i]][dataframes[[i]][, 2] < 0.05/540, ]
}
View(dataframes)




# Function to change column names
change_colnames <- function(df) {
  colnames(df) <- c("EC", "p_value")
  return(df)
}

# Apply the function to all dataframes
dataframes <- lapply(dataframes, change_colnames)

# Function to merge dataframes by surname
merge_dataframes <- function(df) {
  merged_df <- merge(df, ECs_bact, by = "EC")
  return(merged_df)
}

# Apply the function to each dataframe in the list
merged_dataframes <- lapply(dataframes, merge_dataframes)

View(merged_dataframes)
View(BMI_U_vs_Ov)
EC_names_BMI_U_vs_Ob <- merged_dataframes[[1]]
EC_names_BMI_U_vs_H <- merged_dataframes[[2]]
EC_names_BMI_U_vs_Ov <- merged_dataframes[[3]]
EC_names_BMI_H_vs_Ov <- merged_dataframes[[4]]
EC_names_BMI_H_vs_Ob <- merged_dataframes[[5]]
EC_names_BMI_Ov_vs_Ob<- merged_dataframes[[6]]




BMI_U_vs_Ob
BMI_U_vs_H
BMI_U_vs_Ov
BMI_H_vs_Ov
BMI_H_vs_Ob
BMI_Ov_vs_Ob






#summmary statistics of population:
median(DAG3_50_CLR$ANTHRO.BMI)
median(DAG3_50_CLR$ANTHRO.AGE)
fivenum(DAG3_50_CLR$ANTHRO.BMI)

mean(DAG3_50_CLR$ANTHRO.BMI)
table(DAG3_50_CLR$ANTHRO.Sex)
mean(as.numeric(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Discomfort.Unwell))
mean_value <- mean(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Discomfort.Unwell, na.rm = TRUE)
symptoms_mean <- colMeans(DAG3_50_CLR[, 29:34], na.rm = TRUE)

column_range <- 29:34
install.packages("matrixStats")
library(matrixStats)
# Assuming you have a data frame called DAG3_50_CLR and a column range called column_range

# Calculate the median, minimum, maximum, and standard deviation for columns 29 to 34, ignoring missing values
column_medians <- colMedians(as.matrix(DAG3_50_CLR[, column_range]), na.rm = TRUE)
column_mins <- apply(DAG3_50_CLR[, column_range], 2, min, na.rm = TRUE)
column_maxs <- apply(DAG3_50_CLR[, column_range], 2, max, na.rm = TRUE)
column_stddevs <- apply(DAG3_50_CLR[, column_range], 2, sd, na.rm = TRUE)

# Get the column names
column_names <- names(DAG3_50_CLR[column_range])

# Combine column names, medians, minimums, maximums, and standard deviations
result <- data.frame(Column = column_names, Median = column_medians, Minimum = column_mins, Maximum = column_maxs, StdDev = column_stddevs)


column_range <- 29:34

# Calculate the number and percentage of values greater than 1 for each column
count_greater_than_1 <- apply(DAG3_50_CLR[, column_range], 2, function(x) sum(x >= 1, na.rm = TRUE))
percentage_greater_than_1 <- count_greater_than_1 / nrow(DAG3_50_CLR) * 100

# Get the column names
column_names <- names(DAG3_50_CLR[column_range])

# Combine column names, count, and percentage
result <- data.frame(
  Column = column_names,
  CountGreaterThan1 = count_greater_than_1,
  PercentageGreaterThan1 = percentage_greater_than_1
)

# Print the result
View(result)


#the same thing but with gender

gender_F_md <- DAG3md[DAG3md$ANTHRO.Sex == "F", ]
gender_F_md <- gender_F_md[complete.cases(gender_F_md$ANTHRO.Sex), ]

gender_M_md <- DAG3md[DAG3md$ANTHRO.Sex == "M", ]
gender_M_md <- gender_M_md[complete.cases(gender_M_md$ANTHRO.Sex), ]

#combining these sorted md with CAZymes values from clr_DAG3_50

gender_F <- merge(gender_F_md, clr_DAG3_50, by = "row.names")
gender_M <- merge(gender_M_md , clr_DAG3_50, by = "row.names")



#gender_F vs gender_M
df1_W <- gender_F
df2_W <- gender_M
result_df <- data.frame()
for (col in 37:126) {
  test_result <- wilcox.test(df1_W[, col], df2_W[, col], paired = FALSE)
  p_value <- test_result$p.value
  row <- data.frame(Column = colnames(df1_W)[col], P_Value = p_value)
  result_df <- rbind(result_df, row)
}
Gender_F_vs_M <- result_df
Gender_F_vs_M <- Gender_F_vs_M[order(Gender_F_vs_M$P_Value), ]
significant_Wilc_gender <- Gender_F_vs_M$Column[Gender_F_vs_M$P_Value < 0.05]
significant_Wilc_gender_adj <- Gender_F_vs_M$Column[Gender_F_vs_M$P_Value < 0.05/90]

significant_Wilc_gender_adj
sort(significant_Wilc_gender_adj)
View(Gender_F_vs_M )
significant_Wilc_gender
significant_Wilc_gender_adj <- c(significant_Wilc_gender_adj)

gender_wilc_EC <- ECs_bact[ECs_bact$EC %in% significant_Wilc_gender_adj,]


colnames(Gender_F_vs_M) <- c("EC", "p_value")
Gender_EC_p_value_ad <- merge(Gender_F_vs_M, gender_wilc_EC, by = "EC")
#final
View(Gender_EC_p_value_ad)






















#keep only selected metadata from DAG3md

DAG3md_selected_names <- colnames(DAG3md)
DAG3md_selected_names
DAG3md_selected <- DAG3md[, c(1, 2, 4, 11, 12, 13, 14, 17, 18, 21, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35)]

df <- DAG3md[, c(23, 24, 25)]
colnames(df) <- c("IBS_C", "IBS_M","IBS_D")

#Function to check symptom overlaps
check_symptom_overlap <- function(df) {
  # Check if any samples have all three symptoms
  all_three_symptoms <- df[df$IBS_C == "Y" & df$IBS_M == "Y" & df$IBS_D == "Y", "sample"]
  
  # Check if any samples have the first and third symptoms
  first_third_symptoms <- df[df$IBS_C == "Y" & df$IBS_D == "Y" & df$IBS_M != "Y", "sample"]
  
  # Check if any samples have the second and third symptoms
  second_third_symptoms <- df[df$IBS_M == "Y" & df$IBS_D == "Y" & df$IBS_C != "Y", "sample"]
  
  # Check if any samples have the first and second symptoms
  first_second_symptoms <- df[df$IBS_C == "Y" & df$IBS_M == "Y" & df$IBS_D != "Y", "sample"]
  
  # Return the overlapping samples for each scenario
  return(
    list(
      all_three_symptoms = all_three_symptoms,
      first_third_symptoms = first_third_symptoms,
      second_third_symptoms = second_third_symptoms,
      first_second_symptoms = first_second_symptoms
    )
  )
}

# Call the function
overlap_samples <- check_symptom_overlap(df)

# Print the overlapping samples
print(overlap_samples)


table(DAG3md_selected[, c(10)])









#"BMI_U_vs_H","BMI_U_vs_Ov","BMI_U_vs_Ob","BMI_H_vs_Ov", "BMI_H_vs_Ob","BMI_Ov_vs_Ob"












table(DAG3md$MED.DISEASES.Gastrointestinal.Rome3_IBS.Any)
table(DAG3md$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeD)




#fivenum for filter 100 and plot
five_num_DAG3_values_100 <- as.data.frame(sapply(DAG3_values_100, fivenum))
row_stat_names <- c('The minimum', 'The first quartile', 'The median', 'The third quartile', 'The maximum')
rownames(five_num_DAG3_values_100) <- row_stat_names
require(reshape2)
ggplot(data = melt(DAG3_values_100), aes(x=variable, y=value)) +
geom_boxplot(aes(fill=variable)) + theme_bw() +
theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
labs(title = "Summary of selected CAZymes (filter 100)", x = "CAZymes", y = "Count") +
theme(legend.position="none")

#fivenum for filter 50 and plot
five_num_DAG3_values_50 <- as.data.frame(sapply(DAG3_values_50, fivenum))
row_stat_names <- c('The minimum', 'The first quartile', 'The median', 'The third quartile', 'The maximum')
rownames(five_num_DAG3_values_50) <- row_stat_names
ggplot(data = melt(DAG3_values_50), aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=variable)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(title = "Summary of selected CAZymes (filter 50)", x = "CAZymes", y = "Count") +
  theme(legend.position="none")

#DAG3_values_100 20 cazymes for 8204 samples
#DAG3_values_50 90 cazymes for 8204 samples

#checking md for normality

install.packages("nortest")
library(nortest)
# Perform Anderson-Darling and Q-Q plot
ad.test(DAG3md$ANTHRO.AGE)
qqnorm(DAG3md$ANTHRO.AGE)
qqline(DAG3md$ANTHRO.AGE)
fivenum(DAG3md$ANTHRO.AGE)

ggplot(DAG3md, aes(x = DAG3md$ANTHRO.AGE)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  labs(x = "Age", y = "People", title = "Age Distribution") + theme_bw()


ad.test(DAG3md$ANTHRO.BMI)
qqnorm(DAG3md$ANTHRO.BMI)
qqline(DAG3md$ANTHRO.BMI)
fivenum(DAG3md$ANTHRO.BMI)

ggplot(DAG3md, aes(x = DAG3md$ANTHRO.BMI)) +
  geom_histogram(binwidth = 0.5, fill = "steelblue", color = "white") +
  labs(x = "Age", y = "People", title = "BMI Distribution") + theme_bw()



DAG3md$BMI_Group <- cut(DAG3md$ANTHRO.BMI, 
                        breaks = c(0, 18.5, 25, 30, Inf),
                        labels = c("Underweight", "Healthy Weight", "Overweight", "Obesity"))


ggplot(DAG3md, aes(x = BMI_Group, fill = BMI_Group)) +
  geom_bar() +
  labs(x = "BMI Group", y = "People", title = "BMI Distribution by Group in DAG3") + theme_bw() +
  scale_fill_manual(values = c("Underweight" = "blue",
                               "Healthy Weight" = "green",
                               "Overweight" = "yellow",
                               "Obesity" = "red"))


subset_df <- DAG3_values_50[1:15, 1:50]
write_delim(subset_df, "subset_Df.txt", sep="\t", quote = F)

write.table(subset_df, "checking_phenotypes_NEXT.txt", sep = "\t", quote = F, row.names = T)
#clustering Ecs by first 3 digits:
df <- DAG3_values_50

df <- data.frame(df, check.names = FALSE)

# Sample dataframe
df <- data.frame(
  SampleID = c("Sample1", "Sample2", "Sample3"),
  `2.3.10.10` = c(10, 5, 8),
  `2.4.4.108` = c(6, 3, 2),
  `2.4.40.10` = c(4, 7, 9),
  `2.4.40.11` = c(4, 7, 9),
  check.names = FALSE
)


DAG3_50_CLR <- merge( DAG3md, clr_DAG3_50, by = "row.names")
row.names(DAG3_50_CLR) <- DAG3_50_CLR$Row.names
DAG3_50_CLR <- DAG3_50_CLR[,-1]
# 1- age, 2- bmi, 4- gender, 
# 87 - 3.2.1.23 
# 38 - 2.4.1.11

test_CLR <- DAG3_50_CLR[,c(1,2,4,38,87)]

#AGE vs CLR(RPKM)
# Create the scatter plot for 3.2.1.23  cazyme and age
ggplot(test_CLR, aes(x = test_CLR[,1], y = test_CLR[,5])) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Age", y = "CLR(RPKM)", title = "Age and CLR(RPKM) of CAZyme 3.2.1.23 ") +
  theme_minimal()
x <- test_CLR[,1]
y <- test_CLR[,5]
model <- lm(y ~ x)
correlation <- cor.test(x, y)$estimate
p_value <- summary(model)$coefficients[2, "Pr(>|t|)"]
correlation
p_value

#BMI vs CLR(RPKM)
# Create the scatter plot for 3.2.1.23  cazyme and age
ggplot(test_CLR, aes(x = test_CLR[,2], y = test_CLR[,5])) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "BMI", y = "CLR(RPKM)", title = "BMI and CLR(RPKM) of CAZyme 3.2.1.23 ") +
  theme_minimal()
x <- test_CLR[,2]
y <- test_CLR[,5]
model <- lm(y ~ x)
correlation <- cor.test(x, y)$estimate
p_value <- summary(model)$coefficients[2, "Pr(>|t|)"]
correlation
p_value

#gender vs CLR(RPLM)
#wilcoxon gender final for CLR(RPLM)

View(DAG3_50_CLR)
#gender_F vs gender_M
df1_W <- DAG3_50_CLR[which(DAG3_50_CLR$ANTHRO.Sex == "F"), ]
df2_W <- DAG3_50_CLR[which(DAG3_50_CLR$ANTHRO.Sex == "M"), ]
result_df <- data.frame()
for (col in 35:124) {
  test_result <- wilcox.test(DAG3_50_CLR[,col]~ DAG3_50_CLR$ANTHRO.Sex, paired = FALSE)
  #test_result <- wilcox.test(df1_W[, col], df2_W[, col], paired = FALSE)
  p_value <- test_result$p.value
  direction <- ifelse(median(df1_W[, col]) > median(df2_W[, col]), "Decrease", "Increase")
  row <- data.frame(EC = colnames(df1_W)[col], P_Value = p_value, direction = direction)
  result_df <- rbind(result_df, row)
}
Gender_F_vs_M <- result_df
Gender_F_vs_M <- Gender_F_vs_M[order(Gender_F_vs_M$P_Value), ] #no LM correction yet
View(Gender_F_vs_M)



#wilcoxon FGIDs final for CLR(RPLM)
#FGID Bloating
df1_W <- DAG3_50_CLR[which(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.FGID.Bloating == "Y"), ]
df2_W <- DAG3_50_CLR[which(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.FGID.Bloating == "N"), ]
result_df <- data.frame()
for (col in 35:124) {
  test_result <- wilcox.test(df1_W[, col], df2_W[, col], paired = FALSE)
  p_value <- test_result$p.value
  direction <- ifelse(median(df1_W[, col]) > median(df2_W[, col]), "Decrease", "Increase")
  row <- data.frame(EC = colnames(df1_W)[col], p_value = p_value, direction = direction)
  result_df <- rbind(result_df, row)
}
FGID_Bloating_wilc <- result_df
FGID_Bloating_wilc <- FGID_Bloating_wilc[order(FGID_Bloating_wilc$p_value), ] #no LM correction yet
View(FGID_Bloating_wilc)

#wilcoxon FGIDs final for CLR(RPLM)
#FGID Constipation
df1_W <- DAG3_50_CLR[which(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.FGID.Constipation == "Y"), ]
df2_W <- DAG3_50_CLR[which(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.FGID.Constipation == "N"), ]
result_df <- data.frame()
for (col in 35:124) {
  test_result <- wilcox.test(df1_W[, col], df2_W[, col], paired = FALSE)
  p_value <- test_result$p.value
  direction <- ifelse(median(df1_W[, col]) > median(df2_W[, col]), "Decrease", "Increase")
  row <- data.frame(EC = colnames(df1_W)[col], p_value = p_value, direction = direction)
  result_df <- rbind(result_df, row)
}
FGID_Constipation_wilc <- result_df
FGID_Constipation_wilc <- FGID_Constipation_wilc[order(FGID_Constipation_wilc$p_value), ] #no LM correction yet
View(FGID_Constipation_wilc)

#wilcoxon FGIDs final for CLR(RPLM)
#FGID Diarrhea
df1_W <- DAG3_50_CLR[which(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.FGID.Diarrhea == "Y"), ]
df2_W <- DAG3_50_CLR[which(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.FGID.Diarrhea == "N"), ]
result_df <- data.frame()
for (col in 35:124) {
  test_result <- wilcox.test(df1_W[, col], df2_W[, col], paired = FALSE)
  p_value <- test_result$p.value
  direction <- ifelse(median(df1_W[, col]) > median(df2_W[, col]), "Decrease", "Increase")
  row <- data.frame(EC = colnames(df1_W)[col], p_value = p_value, direction = direction)
  result_df <- rbind(result_df, row)
}
FGID_Diarrhea_wilc <- result_df
FGID_Diarrhea_wilc <- FGID_Diarrhea_wilc[order(FGID_Diarrhea_wilc$p_value), ] #no LM correction yet
View(FGID_Diarrhea_wilc)















#significant_Wilc_gender <- Gender_F_vs_M$Column[Gender_F_vs_M$P_Value < 0.05]
#significant_Wilc_gender_adj <- Gender_F_vs_M$Column[Gender_F_vs_M$P_Value < 0.05/90]
#significant_Wilc_gender_adj
#sort(significant_Wilc_gender_adj)
#View(Gender_F_vs_M )
#significant_Wilc_gender
#significant_Wilc_gender_adj <- c(significant_Wilc_gender_adj)
#gender_wilc_EC <- ECs_bact[ECs_bact$EC %in% significant_Wilc_gender_adj,]


#SPEARMAN correlation

#final Spearman age on DAG3_50_CLR
# Initialize an empty data frame to store the results
spearman_age_result <- data.frame(CAZyme = character(),
                                  p_value = numeric(),
                                  rho = numeric(),
                                  stringsAsFactors = FALSE)

# Iterate over each CAZyme column
for (i in 35:124) {
  cazyme <- colnames(DAG3_50_CLR)[i]  # Get the CAZyme name
  x <- DAG3_50_CLR[, 1]  # Age column
  y <- DAG3_50_CLR[, i]  # CLR(RPKM) column for the current CAZyme
  
  # Calculate Spearman correlation
  result <- cor.test(x, y, method = "spearman")
  
  # Store the results in the data frame
  spearman_age_result <- rbind(spearman_age_result,
                               data.frame(CAZyme = cazyme,
                                          p_value = result$p.value,
                                          rho = result$estimate,
                                          stringsAsFactors = FALSE))
}

# Print the resulting data frame
print(spearman_age_result)
colnames(spearman_age_result) <- c("EC", "p_value", "rho")
#View(spearman_age_result)
spearman_age_result #no correction for multiple testing yet
#spearman_age_result_filt <- spearman_age_result[spearman_age_result$p_value < 0.05/90, ]
#spearman_age_result_filt <- merge(spearman_age_result_filt, ECs_bact, by = "EC")
#spearman_age_result_filt <- spearman_age_result_filt[order(spearman_age_result_filt$p_value),]


#final Spearman BMI on DAG3_50_CLR
# Initialize an empty data frame to store the results
spearman_BMI_result <- data.frame(CAZyme = character(),
                                  p_value = numeric(),
                                  rho = numeric(),
                                  stringsAsFactors = FALSE)

# Iterate over each CAZyme column
for (i in 35:124) {
  cazyme <- colnames(DAG3_50_CLR)[i]  # Get the CAZyme name
  x <- DAG3_50_CLR[, 2]  # Age column
  y <- DAG3_50_CLR[, i]  # CLR(RPKM) column for the current CAZyme
  
  # Calculate Spearman correlation
  result <- cor.test(x, y, method = "spearman")
  
  # Store the results in the data frame
  spearman_BMI_result <- rbind(spearman_BMI_result,
                               data.frame(CAZyme = cazyme,
                                          p_value = result$p.value,
                                          rho = result$estimate,
                                          stringsAsFactors = FALSE))
}

# Print the resulting data frame
print(spearman_BMI_result)
colnames(spearman_BMI_result) <- c("EC", "p_value", "rho")
spearman_BMI_result #no LM correction yet


#SPEARMAN ON SYMPTOMS

#checking md for normality
# Perform Anderson-Darling and Q-Q plot
ad.test(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Discomfort.Unwell)
qqnorm(DAG3md$ANTHRO.AGE)
qqline(DAG3md$ANTHRO.AGE)
fivenum(DAG3md$ANTHRO.AGE)
fivenum(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Discomfort.Bloating)

ggplot(DAG3_50_CLR, aes(x = DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Discomfort.Unwell)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "white") +
  labs(x = "Discomfort Unwell rate", y = "People", title = "Discomfort Unwell Distribution") + theme_bw()

ggplot(DAG3_50_CLR, aes(x = DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Discomfort.Pain)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "white") +
  labs(x = "Discomfort Pain rate", y = "People", title = "Discomfort Pain Distribution") + theme_bw()

ggplot(DAG3_50_CLR, aes(x = DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Discomfort.Bloating)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "white") +
  labs(x = "Discomfort Bloating rate", y = "People", title = "Discomfort Bloating Distribution") + theme_bw()

ggplot(DAG3_50_CLR, aes(x = DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Discomfort.Flatulence)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue", color = "white") +
  labs(x = "Discomfort Flatulence rate", y = "People", title = "Discomfort Flatulence Distribution") + theme_bw()

ggplot(DAG3_50_CLR, aes(x = DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Discomfort.Nausea)) +
  geom_histogram(binwidth = 0.09, fill = "steelblue", color = "white") +
  labs(x = "Discomfort Nausea rate", y = "People", title = "Discomfort Nausea Distribution") + theme_bw()

fivenum(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Discomfort.Unwell)
fivenum(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Discomfort.Pain)

fivenum(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Discomfort.Bloating)
fivenum(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Discomfort.Flatulence)

fivenum(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Discomfort.Nausea)



#Sprearman on symptoms Unwell (29)
#final Spearman Unwell on DAG3_50_CLR
# Initialize an empty data frame to store the results
spearman_unwell_result <- data.frame(CAZyme = character(),
                                  p_value = numeric(),
                                  rho = numeric(),
                                  stringsAsFactors = FALSE)

# Iterate over each CAZyme column
for (i in 35:124) {
  cazyme <- colnames(DAG3_50_CLR)[i]  # Get the CAZyme name
  x <- DAG3_50_CLR[, 29]  # Unwell column
  y <- DAG3_50_CLR[, i]  # CLR(RPKM) column for the current CAZyme
  
  # Calculate Spearman correlation
  result <- cor.test(x, y, method = "spearman")
  
  # Store the results in the data frame
  spearman_unwell_result <- rbind(spearman_unwell_result,
                               data.frame(CAZyme = cazyme,
                                          p_value = result$p.value,
                                          rho = result$estimate,
                                          stringsAsFactors = FALSE))
}

# Print the resulting data frame
#print(spearman_unwell_result)
colnames(spearman_unwell_result) <- c("EC", "p_value", "rho")
View(spearman_unwell_result)
#spearman_unwell_result #no correction for multiple testing yet

#Sprearman on symptoms Pain (30)
#final Spearman Pain on DAG3_50_CLR
# Initialize an empty data frame to store the results
spearman_pain_result <- data.frame(CAZyme = character(),
                                     p_value = numeric(),
                                     rho = numeric(),
                                     stringsAsFactors = FALSE)

# Iterate over each CAZyme column
for (i in 35:124) {
  cazyme <- colnames(DAG3_50_CLR)[i]  # Get the CAZyme name
  x <- DAG3_50_CLR[, 30]  # Pain column
  y <- DAG3_50_CLR[, i]  # CLR(RPKM) column for the current CAZyme
  
  # Calculate Spearman correlation
  result <- cor.test(x, y, method = "spearman")
  
  # Store the results in the data frame
  spearman_pain_result <- rbind(spearman_pain_result,
                                  data.frame(CAZyme = cazyme,
                                             p_value = result$p.value,
                                             rho = result$estimate,
                                             stringsAsFactors = FALSE))
}

# Print the resulting data frame
colnames(spearman_pain_result) <- c("EC", "p_value", "rho")
#View(spearman_pain_result)
#spearman_pain_result #no correction for multiple testing yet





#Sprearman on symptoms Bloating (31)
#final Spearman Bloating on DAG3_50_CLR
# Initialize an empty data frame to store the results
spearman_bloating_result <- data.frame(CAZyme = character(),
                                   p_value = numeric(),
                                   rho = numeric(),
                                   stringsAsFactors = FALSE)
# Iterate over each CAZyme column
for (i in 35:124) {
  cazyme <- colnames(DAG3_50_CLR)[i]  # Get the CAZyme name
  x <- DAG3_50_CLR[, 31]  # Age column
  y <- DAG3_50_CLR[, i]  # CLR(RPKM) column for the current CAZyme
  
  # Calculate Spearman correlation
  result <- cor.test(x, y, method = "spearman")
  
  # Store the results in the data frame
  spearman_bloating_result <- rbind(spearman_bloating_result,
                                data.frame(CAZyme = cazyme,
                                           p_value = result$p.value,
                                           rho = result$estimate,
                                           stringsAsFactors = FALSE))
}

# Print the resulting data frame
#print(spearman_bloating_result)
colnames(spearman_bloating_result) <- c("EC", "p_value", "rho")
#View(spearman_bloating_result)
#spearman_pain_result #no correction for multiple testing yet


#Sprearman on symptoms Flatulence (32)
#final Spearman Bloating on DAG3_50_CLR
# Initialize an empty data frame to store the results
spearman_flatulence_result <- data.frame(CAZyme = character(),
                                       p_value = numeric(),
                                       rho = numeric(),
                                       stringsAsFactors = FALSE)
# Iterate over each CAZyme column
for (i in 35:124) {
  cazyme <- colnames(DAG3_50_CLR)[i]  # Get the CAZyme name
  x <- DAG3_50_CLR[, 32]  # Age column
  y <- DAG3_50_CLR[, i]  # CLR(RPKM) column for the current CAZyme
  
  # Calculate Spearman correlation
  result <- cor.test(x, y, method = "spearman")
  
  # Store the results in the data frame
  spearman_flatulence_result <- rbind(spearman_flatulence_result,
                                    data.frame(CAZyme = cazyme,
                                               p_value = result$p.value,
                                               rho = result$estimate,
                                               stringsAsFactors = FALSE))
}

# Print the resulting data frame
colnames(spearman_flatulence_result) <- c("EC", "p_value", "rho")
View(spearman_flatulence_result)
#spearman_pain_result #no correction for multiple testing yet




#Sprearman on symptoms Burping (33)
#final Spearman Bloating on DAG3_50_CLR
# Initialize an empty data frame to store the results
spearman_burping_result <- data.frame(CAZyme = character(),
                                         p_value = numeric(),
                                         rho = numeric(),
                                         stringsAsFactors = FALSE)
# Iterate over each CAZyme column
for (i in 35:124) {
  cazyme <- colnames(DAG3_50_CLR)[i]  # Get the CAZyme name
  x <- DAG3_50_CLR[,33]  # Age column
  y <- DAG3_50_CLR[, i]  # CLR(RPKM) column for the current CAZyme
  
  # Calculate Spearman correlation
  result <- cor.test(x, y, method = "spearman")
  
  # Store the results in the data frame
  spearman_burping_result <- rbind(spearman_burping_result,
                                      data.frame(CAZyme = cazyme,
                                                 p_value = result$p.value,
                                                 rho = result$estimate,
                                                 stringsAsFactors = FALSE))
}

# Print the resulting data frame
colnames(spearman_burping_result) <- c("EC", "p_value", "rho")
#View(spearman_burping_result)


#Sprearman on symptoms Nausea (34)
#final Spearman Bloating on DAG3_50_CLR
# Initialize an empty data frame to store the results
spearman_nausea_result <- data.frame(CAZyme = character(),
                                      p_value = numeric(),
                                      rho = numeric(),
                                      stringsAsFactors = FALSE)
# Iterate over each CAZyme column
for (i in 35:124) {
  cazyme <- colnames(DAG3_50_CLR)[i]  # Get the CAZyme name
  x <- DAG3_50_CLR[,34]  # Age column
  y <- DAG3_50_CLR[, i]  # CLR(RPKM) column for the current CAZyme
  
  # Calculate Spearman correlation
  result <- cor.test(x, y, method = "spearman")
  
  # Store the results in the data frame
  spearman_nausea_result <- rbind(spearman_nausea_result,
                                   data.frame(CAZyme = cazyme,
                                              p_value = result$p.value,
                                              rho = result$estimate,
                                              stringsAsFactors = FALSE))
}

# Print the resulting data frame
colnames(spearman_nausea_result) <- c("EC", "p_value", "rho")
#View(spearman_nausea_result)






#spearman_BMI_result_filt <- spearman_BMI_result[spearman_BMI_result$p_value < 0.05/90, ]
#spearman_BMI_result_filt <- merge(spearman_BMI_result_filt, ECs_bact, by = "EC")
#spearman_BMI_result_filt <- spearman_BMI_result_filt[order(spearman_BMI_result_filt$p_value),]
#BMI <- spearman_BMI_result_filt$EC
#age <- spearman_age_result_filt$EC
#install.packages("VennDiagram")

#LM correction

#Gender_F_vs_M is results from wilcoxon


colnames(Gender_F_vs_M) <- c("EC", "p_value", "direction")
Gender_F_vs_M_wilc <- Gender_F_vs_M


FGID_Bloating_wilc$rho <- c("")
FGID_Constipation_wilc$rho <- c("")
FGID_Diarrhea_wilc$rho <- c("")
Gender_F_vs_M_wilc$rho <- c("")
spearman_BMI_result$direction <- c("")
spearman_age_result$direction <- c("")
spearman_unwell_result$direction <- c("")
spearman_pain_result$direction <- c("")
spearman_bloating_result$direction <- c("")
spearman_flatulence_result$direction <- c("")
spearman_burping_result$direction <- c("")
spearman_nausea_result$direction <- c("")

FGID_Bloating_wilc$test <- c("FGID_Bloating_wilc")
FGID_Constipation_wilc$test <- c("FGID_Constipation_wilc")
FGID_Diarrhea_wilc$test <- c("FGID_Diarrhea_wilc")

Gender_F_vs_M_wilc$test <- c("Sex_wilc")
spearman_BMI_result$test <- c("spearman_BMI")
spearman_age_result$test <- c("spearman_age")

spearman_unwell_result$test <- c("spearman_unwell")
spearman_pain_result$test <- c("spearman_pain")
spearman_bloating_result$test <- c("spearman_bloating")

spearman_flatulence_result$test <- c("spearman_flatulence")
spearman_burping_result$test <- c("spearman_burping")
spearman_nausea_result$test <- c("spearman_nausea")


#FGID_Bloating_wilc, FGID_Constipation_wilc, FGID_Diarrhea_wilc, 
#Gender_F_vs_M_wilc, spearman_BMI_result, spearman_age_result, 
#spearman_unwell_result, spearman_pain_result, spearman_bloating_result, 
#spearman_flatulence_result, spearman_burping_result, spearman_nausea_result

# Create an empty data frame to store the combined results
phenotype_results <- data.frame(test = character(),
                                EC = character(),
                                p_value = numeric(),
                                direction = character(),
                                rho = numeric(),
                                stringsAsFactors = FALSE)

# List of the 12 data frames
data_frames_list <- list(FGID_Bloating_wilc, FGID_Constipation_wilc, FGID_Diarrhea_wilc,
                         Gender_F_vs_M_wilc, spearman_BMI_result, spearman_age_result,
                         spearman_unwell_result, spearman_pain_result, spearman_bloating_result,
                         spearman_flatulence_result, spearman_burping_result, spearman_nausea_result)

# Verify column names in each data frame and combine the data frames
for (df in data_frames_list) {
  required_columns <- c("test", "EC", "p_value", "direction", "rho")
  missing_columns <- setdiff(required_columns, colnames(df))
  
  if (length(missing_columns) > 0) {
    stop(paste("Missing column(s) in the data frame:", paste(missing_columns, collapse = ", ")))
  }
  
  phenotype_results <- rbind(phenotype_results, df[, required_columns])
}

View(phenotype_results)


phenotype_results$BH_FDR_p_value <- p.adjust(phenotype_results$p_value, method = "BH")
phenotype_results$BF_p_value <- p.adjust(phenotype_results$p_value, method="bonferroni")

phenotype_results_prot_names <- merge(phenotype_results, ECs_bact_clean, by = "EC")
write.csv(phenotype_results_prot_names, file = "phenotype_results_prot_names.csv")
View(phenotype_results_prot_names)
table(phenotype_results_prot_names$test)

wilc_blot <- phenotype_results_prot_names[which(phenotype_results_prot_names$test == "FGID_Bloating_wilc"), ]
wilc_const <- phenotype_results_prot_names[which(phenotype_results_prot_names$test == "FGID_Constipation_wilc"), ]
wilc_diar <- phenotype_results_prot_names[which(phenotype_results_prot_names$test == "FGID_Diarrhea_wilc"), ]
wilc_sex <- phenotype_results_prot_names[which(phenotype_results_prot_names$test == "Sex_wilc"), ]

spear_age <- phenotype_results_prot_names[which(phenotype_results_prot_names$test == "spearman_age"), ]
spear_blot <- phenotype_results_prot_names[which(phenotype_results_prot_names$test == "spearman_bloating"), ]
spear_BMI <- phenotype_results_prot_names[which(phenotype_results_prot_names$test == "spearman_BMI"), ]
spear_burp <- phenotype_results_prot_names[which(phenotype_results_prot_names$test == "spearman_burping"), ]

spear_flat <- phenotype_results_prot_names[which(phenotype_results_prot_names$test == "spearman_flatulence"), ]
spear_naus <- phenotype_results_prot_names[which(phenotype_results_prot_names$test == "spearman_nausea"), ]
spear_pain <- phenotype_results_prot_names[which(phenotype_results_prot_names$test == "spearman_pain"), ]
spear_unwel <- phenotype_results_prot_names[which(phenotype_results_prot_names$test == "spearman_unwell"), ]

write.csv(wilc_blot, file = "wilc_blot.csv")
write.csv(wilc_const, file = "wilc_const.csv")
write.csv(wilc_diar, file = "wilc_diar.csv")
write.csv(wilc_sex, file = "wilc_sex.csv")

write.csv(spear_age, file = "spear_age.csv")
write.csv(spear_blot, file = "spear_blot.csv")
write.csv(spear_BMI, file = "spear_BMI.csv")
write.csv(spear_burp, file = "spear_burp.csv")

write.csv(spear_flat, file = "spear_flat.csv")
write.csv(spear_naus, file = "spear_naus.csv")
write.csv(spear_pain, file = "spear_pain.csv")
write.csv(spear_unwel, file = "spear_unwel.csv")



#heatmap of phenotype_results_prot_names

library(corrplot)


# Reshape the data
spear_all_reshaped <- reshape(spear_all, idvar = "EC", timevar = "test", direction = "wide")

# Extract the rho-values and p-values
rho_values <- as.matrix(spear_all_reshaped[, grep("rho", names(spear_all_reshaped))])
p_values <- as.matrix(spear_all_reshaped[, grep("BH_FDR_p_value", names(spear_all_reshaped))])

# Convert rho-values to numeric
rho_values <- apply(rho_values, 2, as.numeric)

# Set the color scale for p-values
color_scale <- colorRampPalette(c("white", "red"))(100)

# Create the heatmap using corrplot
corrplot(rho_values, method = "color", type = "upper", tl.cex = 0.8, col = color_scale, p.mat = p_values, pch.col = "black")


library(corrplot)

# Reshape the data
spear_all_reshaped <- reshape(spear_all, idvar = "EC", timevar = "test", direction = "wide")

# Extract the rho-values and p-values
rho_values <- as.matrix(spear_all_reshaped[, grep("rho", names(spear_all_reshaped))])
p_values <- as.matrix(spear_all_reshaped[, grep("BH_FDR_p_value", names(spear_all_reshaped))])

# Convert rho-values to numeric
rho_values <- apply(rho_values, 2, as.numeric)

# Set the color scale for p-values
color_scale <- colorRampPalette(c("white", "red"))(100)

# Create the heatmap using corrplot
corrplot(rho_values, method = "color", type = "upper", tl.cex = 0.8, col = color_scale, p.mat = p_values, pch.col = "black",
         addvertlabels = TRUE, addtext = TRUE)



library(corrplot)



# Reshape the data
spear_all_reshaped <- reshape(spear_all, idvar = "EC", timevar = "test", direction = "wide")

# Extract the rho-values and p-values
rho_values <- as.matrix(spear_all_reshaped[, grep("rho", names(spear_all_reshaped))])
p_values <- as.matrix(spear_all_reshaped[, grep("BH_FDR_p_value", names(spear_all_reshaped))])

# Convert rho-values to numeric
rho_values <- apply(rho_values, 2, as.numeric)

# Set the color scale for p-values
color_scale <- colorRampPalette(c("white", "red"))(100)

# Add row and column names to the correlation matrix
rownames(rho_values) <- spear_all_reshaped$EC
colnames(rho_values) <- unique(spear_all_reshaped$test)

# Create the heatmap using corrplot
corrplot(rho_values, method = "color", type = "upper", tl.cex = 0.8, col = color_scale, p.mat = p_values, pch.col = "black", addvertlabels = TRUE, addcolnames = TRUE, number.cex = 0.6)




# Assuming you have p-values and rho-values stored in variables
p_values <- spear_age$BH_FDR_p_value
rho_values <- spear_age$rho

p_values <- as.numeric(p_values)
rho_values <- as.numeric(rho_values)
# Calculate correlation matrix
corr_matrix <- cor(matrix(rho_values, nrow = 2))

# Set up the correlation matrix with p-values
corr_matrix_p <- corr_matrix
p_index <- lower.tri(corr_matrix_p, diag = FALSE)
corr_matrix_p[p_index] <- p_values[1:sum(p_index)]

# Load the corrplot library
library(corrplot)

# Create the heatmap
corrplot(corr_matrix_p, method = "color", type = "lower", tl.cex = 0.7, tl.col = "black")




#multivariable linear regression model 
#IBS_M vs healthy
#exclude all IBS_C and IBS_D to keep IBS_C and healthy
df <- DAG3_50_CLR[which(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeC != "Y"), ]
df_upd <- df[which(df$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeD != "Y"), ]
#delete all columns that are not Age + Sex + BMI + Bristool (Mean) + IBS-M
#IBS_M is 24th in df
c <- c(1,2,4,17,24, 35:124)
df_upd <- df_upd[,c]

library(dplyr)

# Rename CAZyme columns
cazyme_columns <- colnames(df_upd)[6:95]
new_cazyme_columns <- paste0("C", gsub("\\.", "_", cazyme_columns))
colnames(df_upd)[6:95] <- new_cazyme_columns

# Create an empty list to store the regression results
regression_results <- list()

# Loop over the modified CAZyme columns
for (cazyme_col in new_cazyme_columns) {
  # Construct the formula for the regression model
  formula <- paste(cazyme_col, "~ ANTHRO.AGE + ANTHRO.Sex + ANTHRO.BMI + META.POOP.BristolMean + MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeM")
  
  # Fit the multivariable linear regression model
  model <- lm(formula, data = df_upd)
  
  # Store the regression summary in the list
  regression_summary <- summary(model)
  regression_results[[cazyme_col]] <- as.data.frame(t(regression_summary$coefficients))
}


# Combine the regression results into a single data frame
summary_df <- bind_rows(regression_results, .id = "CAZyme")



filtered_summary_df <- summary_df[seq(nrow(summary_df)) %% 4 == 0, ]
filtered_summary_df$CAZyme <- gsub("_", ".", gsub("^C", "", filtered_summary_df$CAZyme))
colnames(filtered_summary_df) <- c("EC", "intercept_p_value", "age_p_value", "sex_p_value", "BMI_p_value", "Bristol_p_value", "IBS_M")
filtered_summary_df$BH_p_value <- p.adjust(filtered_summary_df$IBS_M, method = "BH")
filtered_summary_df$BF_p_value <- p.adjust(filtered_summary_df$IBS_M, method="bonferroni")


filtered_summary_df <- merge(filtered_summary_df, ECs_bact_clean, by = "EC")
mult_model_IBS_M <- filtered_summary_df



#multivariable linear regression model 
#IBS_C vs healthy
#exclude all IBS_M and IBS_D to keep IBS_C and healthy
df <- DAG3_50_CLR[which(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeM != "Y"), ]
df_upd <- df[which(df$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeD != "Y"), ]
#delete all columns that are not Age + Sex + BMI + Bristool (Mean) + IBS-C
#IBS_C is 23rd in df
c <- c(1,2,4,17,23, 35:124)
df_upd <- df_upd[,c]

library(dplyr)

# Rename CAZyme columns
cazyme_columns <- colnames(df_upd)[6:95]
new_cazyme_columns <- paste0("C", gsub("\\.", "_", cazyme_columns))
colnames(df_upd)[6:95] <- new_cazyme_columns

# Create an empty list to store the regression results
regression_results <- list()

# Loop over the modified CAZyme columns
for (cazyme_col in new_cazyme_columns) {
  # Construct the formula for the regression model
  formula <- paste(cazyme_col, "~ ANTHRO.AGE + ANTHRO.Sex + ANTHRO.BMI + META.POOP.BristolMean + MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeC")
  
  # Fit the multivariable linear regression model
  model <- lm(formula, data = df_upd)
  
  # Store the regression summary in the list
  regression_summary <- summary(model)
  regression_results[[cazyme_col]] <- as.data.frame(t(regression_summary$coefficients))
}


# Combine the regression results into a single data frame
summary_df <- bind_rows(regression_results, .id = "CAZyme")

# Print the summary data frame
#View(summary_df)

filtered_summary_df <- summary_df[seq(nrow(summary_df)) %% 4 == 0, ]
filtered_summary_df$CAZyme <- gsub("_", ".", gsub("^C", "", filtered_summary_df$CAZyme))
colnames(filtered_summary_df) <- c("EC", "intercept_p_value", "age_p_value", "sex_p_value", "BMI_p_value", "Bristol_p_value", "IBS_C")
filtered_summary_df$BH_p_value <- p.adjust(filtered_summary_df$IBS_C, method = "BH")
filtered_summary_df$BF_p_value <- p.adjust(filtered_summary_df$IBS_C, method="bonferroni")


filtered_summary_df <- merge(filtered_summary_df, ECs_bact_clean, by = "EC")
#View(filtered_summary_df)
mult_model_IBS_C <- filtered_summary_df





#multivariable linear regression model 
#IBS_C vs healthy
#exclude all IBS_M and IBS_D to keep IBS_C and healthy
df <- DAG3_50_CLR[which(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeM != "Y"), ]
df_upd <- df[which(df$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeD != "Y"), ]
#delete all columns that are not Age + Sex + BMI + Bristool (Mean) + IBS-C
#IBS_C is 23rd in df
c <- c(1,2,4,17,23, 35:124)
df_upd <- df_upd[,c]

library(dplyr)

# Rename CAZyme columns
cazyme_columns <- colnames(df_upd)[6:95]
new_cazyme_columns <- paste0("C", gsub("\\.", "_", cazyme_columns))
colnames(df_upd)[6:95] <- new_cazyme_columns

# Create an empty list to store the regression results
regression_results <- list()

# Loop over the modified CAZyme columns
for (cazyme_col in new_cazyme_columns) {
  # Construct the formula for the regression model
  formula <- paste(cazyme_col, "~ ANTHRO.AGE + ANTHRO.Sex + ANTHRO.BMI + META.POOP.BristolMean + MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeC")
  
  # Fit the multivariable linear regression model
  model <- lm(formula, data = df_upd)
  
  # Store the regression summary in the list
  regression_summary <- summary(model)
  regression_results[[cazyme_col]] <- as.data.frame(t(regression_summary$coefficients))
}
# Combine the regression results into a single data frame
summary_df <- bind_rows(regression_results, .id = "CAZyme")
# Print the summary data frame
#View(summary_df)
filtered_summary_df <- summary_df[seq(nrow(summary_df)) %% 4 == 0, ]
filtered_summary_df$CAZyme <- gsub("_", ".", gsub("^C", "", filtered_summary_df$CAZyme))
colnames(filtered_summary_df) <- c("EC", "intercept_p_value", "age_p_value", "sex_p_value", "BMI_p_value", "Bristol_p_value", "IBS_C")
filtered_summary_df$BH_p_value <- p.adjust(filtered_summary_df$IBS_C, method = "BH")
filtered_summary_df$BF_p_value <- p.adjust(filtered_summary_df$IBS_C, method="bonferroni")

filtered_summary_df <- merge(filtered_summary_df, ECs_bact_clean, by = "EC")
#View(filtered_summary_df)
mult_model_IBS_C <- filtered_summary_df



#multivariable linear regression model 
#IBS_D vs healthy
#exclude all IBS_M and IBS_C to keep IBS_D and healthy
df <- DAG3_50_CLR[which(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeM != "Y"), ]
df_upd <- df[which(df$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeC != "Y"), ]
#delete all columns that are not Age + Sex + BMI + Bristool (Mean) + IBS-C
#IBS_D is 23rd in df
c <- c(1,2,4,17,25, 35:124)
df_upd <- df_upd[,c]
# Rename CAZyme columns
cazyme_columns <- colnames(df_upd)[6:95]
new_cazyme_columns <- paste0("C", gsub("\\.", "_", cazyme_columns))
colnames(df_upd)[6:95] <- new_cazyme_columns

# Create an empty list to store the regression results
regression_results <- list()

# Loop over the modified CAZyme columns
for (cazyme_col in new_cazyme_columns) {
  # Construct the formula for the regression model
  formula <- paste(cazyme_col, "~ ANTHRO.AGE + ANTHRO.Sex + ANTHRO.BMI + META.POOP.BristolMean + MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeD")
  
  # Fit the multivariable linear regression model
  model <- lm(formula, data = df_upd)
  
  # Store the regression summary in the list
  regression_summary <- summary(model)
  regression_results[[cazyme_col]] <- as.data.frame(t(regression_summary$coefficients))
}

# Combine the regression results into a single data frame
summary_df <- bind_rows(regression_results, .id = "CAZyme")

# Print the summary data frame
#View(summary_df)
filtered_summary_df <- summary_df[seq(nrow(summary_df)) %% 4 == 0, ]
filtered_summary_df$CAZyme <- gsub("_", ".", gsub("^C", "", filtered_summary_df$CAZyme))
colnames(filtered_summary_df) <- c("EC", "intercept_p_value", "age_p_value", "sex_p_value", "BMI_p_value", "Bristol_p_value", "IBS_D")
filtered_summary_df$BH_p_value <- p.adjust(filtered_summary_df$IBS_D, method = "BH")
filtered_summary_df$BF_p_value <- p.adjust(filtered_summary_df$IBS_D, method="bonferroni")

filtered_summary_df <- merge(filtered_summary_df, ECs_bact_clean, by = "EC")
#View(filtered_summary_df)
mult_model_IBS_D <- filtered_summary_df

count <- sum(mult_model_IBS_D$BH_p_value < 0.05)







df <- df_upd[,46]

residuals(lm(‘cazyme ~ cofactors’ , data=df, na.action = na.omit))


#correcting values for residuals 
#IBS_D vs healthy

df <- DAG3_50_CLR[which(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeM != "Y"), ]
df_upd <- df[which(df$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeC != "Y"), ]
#delete all columns that are not Age + Sex + BMI + Bristool (Mean) + IBS-C
#IBS_D is 25rd in df
c <- c(1,2,4,17,25, 35:124)
df_upd <- df_upd[,c]
# Rename CAZyme columns
cazyme_columns <- colnames(df_upd)[6:95]
new_cazyme_columns <- paste0("C", gsub("\\.", "_", cazyme_columns))
colnames(df_upd)[6:95] <- new_cazyme_columns

na_rows <- apply(df_upd, 1, function(row) any(is.na(row)))
df_upd <- df_upd[!na_rows, ]

formula <- paste("C2_4_1_211 ~ ANTHRO.AGE + ANTHRO.Sex + ANTHRO.BMI + META.POOP.BristolMean")
residuals_cazyme_D <- as.data.frame(residuals(lm(formula, data = df_upd)))

residuals_cazyme_D$IBS <- df_upd$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeD



df <- DAG3_50_CLR[which(DAG3_50_CLR$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeD != "Y"), ]
df_upd <- df[which(df$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeC != "Y"), ]
#delete all columns that are not Age + Sex + BMI + Bristool (Mean) + IBS-M
#IBS_M is 24th in df 
c <- c(1,2,4,17,24, 35:124)
df_upd <- df_upd[,c]
# Rename CAZyme columns
cazyme_columns <- colnames(df_upd)[6:95]
new_cazyme_columns <- paste0("C", gsub("\\.", "_", cazyme_columns))
colnames(df_upd)[6:95] <- new_cazyme_columns

na_rows <- apply(df_upd, 1, function(row) any(is.na(row)))
df_upd <- df_upd[!na_rows, ]

formula <- paste("C2_4_1_211 ~ ANTHRO.AGE + ANTHRO.Sex + ANTHRO.BMI + META.POOP.BristolMean")
residuals_cazyme_M <- as.data.frame(residuals(lm(formula, data = df_upd)))

residuals_cazyme_M$IBS <- df_upd$MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeM

residuals_cazyme_M$IBS[residuals_cazyme_M$IBS == "Y"] <- "IBS_M"
residuals_cazyme_M$IBS[residuals_cazyme_M$IBS == "N"] <- "No_IBS"

residuals_cazyme_D$IBS[residuals_cazyme_D$IBS == "Y"] <- "IBS_D"
residuals_cazyme_D$IBS[residuals_cazyme_D$IBS == "N"] <- "No_IBS"

residuals_M_only <- residuals_cazyme_M[residuals_cazyme_M$IBS == "IBS_M", ]

residuals_cazyme_MD_Healthy <- rbind(residuals_M_only, residuals_cazyme_D)
colnames(residuals_cazyme_MD_Healthy) <- c("RPKP_CLR_resid_cor", "IBS_status")

# Create a box plot and violin plot using ggplot2
plot <- ggplot(residuals_cazyme_MD_Healthy, aes(x = IBS_status, y = RPKP_CLR_resid_cor)) +
  geom_boxplot(width = 0.3, fill = "lightblue", outlier.shape = NA) +
  geom_violin(fill = "lightgray", alpha = 0.5) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 15, face = "bold"),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15)
  ) +
  labs(x = "CAZyme galacto-N-biose phosphorylase (EC 2.4.1.211)", y = "CLR (RPKM) corrected for residuals", title = "The distribution of the most significant CAZyme for IBS-M, IBS-D and non-IBS")

print(plot)

library(ggsignif)

# Create a box plot and violin plot using ggplot2
plot <- ggplot(residuals_cazyme_MD_Healthy, aes(x = IBS_status, y = RPKP_CLR_resid_cor)) +
  geom_boxplot(width = 0.3, fill = "lightblue", outlier.shape = NA) +
  geom_violin(fill = "lightgray", alpha = 0.5) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 17, face = "bold"),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15)
  ) +
  labs(x = "CAZyme galacto-N-biose phosphorylase (EC 2.4.1.211)", y = "CLR (RPKM) corrected for residuals", title = "The distribution of the most significant CAZyme for IBS-M, IBS-D and non-IBS")

# Add significance stars
plot_with_significance <- plot +
  geom_signif(comparisons = list(c("IBS-M", "non-IBS"), c("IBS-D", "non-IBS")), 
              map_signif_level = FALSE, 
              textsize = 3,
              y_position = c(25, 20),
              tip_length = 0.01)

print(plot_with_significance)


residuals_cazyme_MD_Healthy <- residuals_cazyme_MD_Healthy[complete.cases(residuals_cazyme_MD_Healthy), ]

# Convert IBS_status to a factor variable
residuals_cazyme_MD_Healthy$IBS_status <- factor(residuals_cazyme_MD_Healthy$IBS_status)

# Create a box plot and violin plot using ggplot2
plot <- ggplot(residuals_cazyme_MD_Healthy, aes(x = IBS_status, y = RPKP_CLR_resid_cor)) +
  geom_boxplot(width = 0.3, fill = "lightblue", outlier.shape = NA) +
  geom_violin(fill = "lightgray", alpha = 0.5) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 17, face = "bold"),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_text(size = 15)
  ) +
  labs(x = "CAZyme galacto-N-biose phosphorylase (EC 2.4.1.211)", y = "CLR (RPKM) corrected for residuals", title = "The distribution of the most significant CAZyme for IBS-M, IBS-D and non-IBS")

# Add significance stars
plot_with_significance <- plot +
  geom_signif(comparisons = list(c("IBS-M", "non-IBS"), c("IBS-D", "non-IBS")), 
              map_signif_level = TRUE, 
              textsize = 3,
              y_position = c(25, 20),
              tip_length = 0.01)

print(plot_with_significance)


#heatmap

heatmap <- rbind(spear_age, spear_BMI, spear_blot, spear_burp, spear_flat, spear_pain, spear_unwel)
                 
columm_num <- c(1,2,5,6)
heatmap <- heatmap[,columm_num]

colnames(heatmap) <- c("EC","test","rho","p_BH")



filtered_heatmap <- heatmap[heatmap$p_BH < 0.05, ]
filtered_heatmap$rho <- as.numeric(as.character(filtered_heatmap$rho))

test_mapping <- c("spearman_age" = "Age", "spearman_bloating" = "Bloating", "spearman_BMI" = "BMI", "spearman_burping" = "Burping", "spearman_flatulence" = "Flatulence", "spearman_pain" = "Pain", "spearman_unwell" = "Unwell")
filtered_heatmap <- filtered_heatmap %>%
  mutate(test = recode(test, !!!test_mapping, .default = test))

star_labels <- c("***" = "p-value <= 0.005", "**" = "0.005 < p-value <= 0.01", "*" = "0.01 < p-value < 0.05")

ggplot(filtered_heatmap, aes(x = test, y = EC)) +
  geom_tile(aes(fill = rho), color = "white") +
  geom_text(aes(label = ifelse(p_BH < 0.05, ifelse(p_BH <= 0.005, "***", "**"), "*")),
            color = "black", size = 5, angle = 90, vjust = 0.7) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = "Phenotypes and symptoms", y = "EC numbers of CAZymes") +
  coord_flip()+
  ggtitle("Spearman correlation for phenotypes and symptoms") + theme_bw()+
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )



ggplot(filtered_heatmap, aes(x = test, y = EC)) +
  geom_tile(aes(fill = rho), color = "white") +
  geom_text(aes(label = ifelse(p_BH < 0.05, ifelse(p_BH < 0.005, "***", ifelse(p_BH < 0.0005, "**", "*")), "")),
            color = "black", size = 5, angle = 90, vjust = 0.7) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = "Phenotypes and symptoms", y = "EC numbers of CAZymes") +
  coord_flip()+
  ggtitle("Spearman correlation for phenotypes and symptoms") + theme_bw()+
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )


ggplot(filtered_heatmap, aes(x = test, y = EC)) +
  geom_tile(aes(fill = rho), color = "white") +
  geom_text(aes(label = ifelse(p_BH < 0.05, ifelse(p_BH < 0.005, "***", ifelse(p_BH < 0.0005, "**", "*")), "")),
            color = "black", size = 5, angle = 90, vjust = 0.7) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = "Phenotypes and symptoms", y = "EC numbers of CAZymes") +
  coord_flip() +
  ggtitle("Spearman correlation for phenotypes and symptoms") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  annotate("text", x = Inf, y = Inf, label = "Number of Stars:\n***: p-value < 0.0005\n** : 0.0005 <= p-value < 0.005\n*   : 0.005 <= p-value < 0.05",
           hjust = 1, vjust = 1, size = 5, color = "black")


ggplot(filtered_heatmap, aes(x = test, y = EC)) +
  geom_tile(aes(fill = rho), color = "white") +
  geom_text(aes(label = ifelse(p_BH < 0.05, ifelse(p_BH < 0.005, "***", ifelse(p_BH < 0.0005, "**", "*")), "")),
            color = "black", size = 5, angle = 90, vjust = 0.7) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = "Phenotypes and symptoms", y = "EC numbers of CAZymes") +
  coord_flip() +
  ggtitle("Spearman correlation for phenotypes and symptoms") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    plot.margin = margin(t = 0, r = 1.5, b = 0, l = 0)
  ) +
  annotate("text", x = 1.02, y = 0.5, label = "Number of Stars:\n***: p-value < 0.0005\n** : 0.0005 <= p-value < 0.005\n*   : 0.005 <= p-value < 0.05",
           hjust = 0, vjust = 0.5, size = 5, color = "black")


spear_age$BH_FDR_p_value

count_age <- sum(spear_age$BH_FDR_p_value < 0.05)
count_blot <- sum(spear_blot$BH_FDR_p_value < 0.05)
count_burp <- sum(spear_burp$BH_FDR_p_value < 0.05)
count_flat <- sum(spear_flat$BH_FDR_p_value < 0.05)
count_pain <- sum(spear_pain$BH_FDR_p_value < 0.05)
count_unwell <- sum(spear_unwel$BH_FDR_p_value < 0.05)
count_sex <- sum(wilc_sex$BH_FDR_p_value < 0.05)
subset_df_wilc_sex <- wilc_sex[wilc_sex$BH_FDR_p_value < 0.05, ]

subset_df_spear_BMI <- spear_BMI[spear_BMI$BH_FDR_p_value < 0.05, ]
subset_df_spear_age <- spear_age[spear_age$BH_FDR_p_value < 0.05, ]




mult_model_IBS_D <- mult_model_IBS_D[order(mult_model_IBS_D$BH_p_value), ]
mult_model_IBS_C <- mult_model_IBS_C[order(mult_model_IBS_C$BH_p_value), ]
mult_model_IBS_M <- mult_model_IBS_M[order(mult_model_IBS_M$BH_p_value), ]

row.names(mult_model_IBS_D) <- NULL
row.names(mult_model_IBS_C) <- NULL
row.names(mult_model_IBS_M) <- NULL

write.csv(mult_model_IBS_C, file = "mult_model_IBS_C.csv")
write.csv(mult_model_IBS_M, file = "mult_model_IBS_M.csv")
write.csv(mult_model_IBS_D, file = "mult_model_IBS_D.csv")

c <- c(1,8,10)

mult_model_IBS_C_top_10 <- mult_model_IBS_C[1:10,c]
mult_model_IBS_D_top_10 <- mult_model_IBS_D[1:10,c]
mult_model_IBS_M_top_10 <- mult_model_IBS_M[1:10,c]

mult_model_IBS_D_top_24 <- mult_model_IBS_D[1:24,c]
mult_model_IBS_M_top_24 <- mult_model_IBS_M[1:24,c]

overlap1 <- intersect(mult_model_IBS_M_top_24$EC, mult_model_IBS_D_top_24$EC)
View(mult_model_IBS_D_top_10)
mult_model_IBS_C_top_10
mult_model_IBS_M_top_10

mult_model_IBS_C_top_10$EC

mult_model_IBS_M$BH_p_value

sign_mult_model_IBS_M <- mult_model_IBS_M$EC[mult_model_IBS_M$BH_p_value < 0.05]
sign_mult_model_IBS_D <- mult_model_IBS_M$EC[mult_model_IBS_D$BH_p_value < 0.05]

overlap <- intersect(mult_model_IBS_M_top_10$EC, mult_model_IBS_D_top_10$EC)


#boxplot of top1 IBS_D and IBS_M

values_M <- IBS_M[, 47]
# Create a corresponding dataframe with "name" and "values" columns
df_M <- data.frame(name = "IBS_M", values = values_M)

# Extract values from 47th column of IBS_D
values_D <- IBS_D[, 47]
# Create a corresponding dataframe with "name" and "values" columns
df_D <- data.frame(name = "IBS_D", values = values_D)



# Combine the two dataframes
output_df <- rbind(df_M, df_D)
# Create a box plot using ggplot2
boxplot <- ggplot(output_df, aes(x = name, y = values)) +
  geom_boxplot() +
  theme_bw() +
  xlab("2.4.1.211") +
  ylab("CLR (RPKM)") 
print(boxplot)

# Create a box plot and violin plot using ggplot2
plot <- ggplot(output_df, aes(x = name, y = values)) +
  geom_boxplot(width = 0.3, fill = "lightblue", outlier.shape = NA) +
  geom_violin(fill = "lightgray", alpha = 0.5) +
  theme_bw() +
  labs(x = "CAZyme galacto-N-biose phosphorylase (EC 2.4.1.211)", y = "CLR (RPKM)", title = "The distribution of the most significant CAZyme for IBS-M and IBS-D")


# Display the combined plot
print(plot)


# Create the box plot
ggplot(test, aes(x = gender, y = CAZy_2_4_1_10)) +
  geom_boxplot() +
  labs(x = "Gender", y = "CAZy 2.4.1.10 rpkm", title = "rpkm value of 3.2.1.23 CAZyme by Gender") +
  theme_minimal()




# Display the diagram in RStudio
venn



#list_correction <- list(Gender_F_vs_M, spearman_BMI_result, spearman_age_result, 
#                        spearman_unwell_result, spearman_pain_result, spearman_bloating_result, 
#                        spearman_flatulence_result, spearman_burping_result, spearman_nausea_result)
#names(list_correction)[1] <- "Gender_F_vs_M_wilc"
#names(list_correction)[2] <- "spearman_BMI_result"
#names(list_correction)[3] <- "spearman_age_result"

#names(list_correction)[4] <- "spearman_unwell_result"
#names(list_correction)[5] <- "spearman_pain_result"
#names(list_correction)[6] <- "spearman_bloating_result"

#names(list_correction)[7] <- "spearman_flatulence_result"
#names(list_correction)[8] <- "spearman_burping_result"
#names(list_correction)[9] <- "spearman_nausea_result"


#list_correction_BH <- lapply(list_correction, function(df) {
#  df$BH_FDR_p_value <- p.adjust(df$p_value, method = "fdr")
#  return(df)
#})

#correlation_coefficient <- 0.5
#number_of_tests <- sum(sapply(list_correction, function(df) length(df$p_value)))
#correction_factor <- sum(1 / (1 - correlation_coefficient)) / number_of_tests

#list_correction_BH_LM <- lapply(list_correction_BH, function(df) {
#  df$LM_p_value <- df$p_value * correction_factor
#  return(df)
#})


write.csv(list_correction_BH_LM, file = "list_correction_BH_LM.csv")


#View(list_correction[[1]])
# Define the correlation coefficient
#correlation_coefficient <- 0.5

# Calculate the correction factor
#number_of_tests <- sum(sapply(list_correction, function(df) length(df$p_value)))
#correction_factor <- sum(1 / (1 - correlation_coefficient)) / number_of_tests

# Apply LM correction to each data frame in the list
#list_correction <- lapply(list_correction, function(df) {
#  df$corrected_p_values <- df$p_value * correction_factor
#  return(df)
#})


#ECs_bact[12:30,1]
ECs_bact$EC <- trimws(ECs_bact$EC)
ECs_bact <- ECs_bact[grepl('\\d+\\.\\d+\\.\\d+\\.\\d+', ECs_bact$EC), ]
ECs_bact_clean <- ECs_bact[!duplicated(ECs_bact$EC), ]
list_correction_BH_LM_name <- lapply(list_correction_BH_LM, function(df) merge(df, ECs_bact_clean, by = "EC"))
write.csv(list_correction_BH_LM_name, file = "list_correction_BH_LM_name.csv")




library(VennDiagram)

venn_result <- venn.diagram(
  x = list(bmi, age),
  category.names = c("EC_spearman_BMI_result_filt", "EC_spearman_age_result_filt"),
  filename = NULL
)

# Plot the Venn diagram
plot(venn_result)

install.packages("UpSetR")

# Load the UpSetR library
library(UpSetR)

# Create the UpSet plot
upset_result <- upset(fromList(list(BMI = BMI, age = age)),
                      sets = c("BMI", "age"),
                      order.by = "freq",
                      main.bar.color = "steelblue",
                      matrix.color = "steelblue",
                      text.scale = c(1.5, 1.2, 1.2),
                      point.size = 4,
                      show.numbers = TRUE)

# Display the UpSet plot
print(upset_result)



# 137 - 3.2.1.23 
# 36 - 2.4.1.10
test <- DAG3[, c(1, 2, 4, 36, 137)]

#test[,2:3] <- round(test[,2:3], 0)

colnames(test) <- c("age", "bmi","gender", "CAZy_2_4_1_10", "CAZy_3_2_1_23")


# Add a small constant to the RPKM counts
small_constant <- 0.1
test$CAZy_2_4_1_10_counts_adjusted <- test$CAZy_2_4_1_10+ small_constant
test$CAZy_3_2_1_23_counts_adjusted <- test$CAZy_3_2_1_23+ small_constant


#AGE vs RPKM
# Create the scatter plot for first cazyme and age
ggplot(test, aes(x = age, y = CAZy_2_4_1_10_counts_adjusted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Age", y = "RPKM Counts", title = "Age and rpkm of CAZyme 2.4.1.10") +
  theme_minimal()
x <- test$age
y <- test$CAZy_2_4_1_10
model <- lm(y ~ x)
correlation <- cor.test(x, y)$estimate
p_value <- summary(model)$coefficients[2, "Pr(>|t|)"]
correlation
p_value

# Create the scatter plot for second cazyme and age
ggplot(test, aes(x = age, y = CAZy_3_2_1_23_counts_adjusted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "Age", y = "RPKM Counts", title = "Age and rpkm of CAZyme 3.2.1.23") +
  theme_minimal()
x <- test$age
y <- test$CAZy_3_2_1_23
model <- lm(y ~ x)
correlation <- cor.test(x, y)$estimate
p_value <- summary(model)$coefficients[2, "Pr(>|t|)"]
correlation
p_value



#BMI vs RPKM
# Create the scatter plot for first cazyme and BMI
ggplot(test, aes(x = bmi, y = CAZy_2_4_1_10_counts_adjusted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "BMI", y = "RPKM Counts", title = "BMI and rpkm of CAZyme 2.4.1.10") +
  theme_minimal()
x <- test$bmi
y <- test$CAZy_2_4_1_10
model <- lm(y ~ x)
correlation <- cor.test(x, y)$estimate
p_value <- summary(model)$coefficients[2, "Pr(>|t|)"]
correlation
p_value

# Create the scatter plot for second cazyme and BMI
ggplot(test, aes(x = bmi, y = CAZy_3_2_1_23_counts_adjusted)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  labs(x = "BMI", y = "RPKM Counts", title = "BMI and rpkm of CAZyme 3.2.1.23") +
  theme_minimal()
x = test$bmi
y <- test$CAZy_3_2_1_23
model <- lm(y ~ x)
correlation <- cor.test(x, y)$estimate
p_value <- summary(model)$coefficients[2, "Pr(>|t|)"]
correlation
p_value

# Create the box plot
ggplot(test, aes(x = gender, y = CAZy_2_4_1_10)) +
  geom_boxplot() +
  labs(x = "Gender", y = "CAZy 2.4.1.10 rpkm", title = "rpkm value of 3.2.1.23 CAZyme by Gender") +
  theme_minimal()

# Create the box plot
ggplot(test, aes(x = gender, y = CAZy_3_2_1_23)) +
  geom_boxplot() +
  labs(x = "Gender", y = "CAZy 3.2.1.23 rpkm", title = "rpkm value of 3.2.1.23 CAZyme by Gender") +
  theme_minimal()









age <- c(25, 30, 35, 40, 45, 45, 45)
rpkm <- c(1.2, 3.4, 2.1, 0.8, 4.05, 4.0, 3.1)
df <- data.frame(age, rpkm)

# Create heatmap
ggplot(df, aes(x = age, y = rpkm, fill = rpkm)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal()







#3.2.1.23 galactosidase

beta_galactosidase <- df100[,14]
ggplot(data = melt(beta_galactosidase), aes(x=variable, y=value)) + geom_violin()

# Create a violin plot of the third column (z)
ggplot(df, aes(x = "", y = df100[,14])) +
  geom_violin(fill = "lightblue", color = "black") +
  theme_classic() +
  xlab("3.2.1.23 beta-galactosidase") +
  ylab("counts")

#gender and CAZymes


d <- md_dag3
names_dag3_md <- rownames(md_dag3)
rownames(d) <- NULL
md_dag3_1col_names <- cbind(names_dag3_md,d)
female_md  <- md_dag3_1col_names$ANTHRO.Sex == 'F'
female_names <- md_dag3_1col_names$names_dag3_md[female_md]

male_md  <- md_dag3_1col_names$ANTHRO.Sex == 'M'
male_names <- md_dag3_1col_names$names_dag3_md[male_md]

#df100 male/female

# extracting data frame rows
df100_f <- df100[rownames(df100) %in% female_names, ]
df100_m <- df100[rownames(df100) %in% male_names, ]

df_f <- df[rownames(df) %in% female_names, ]
df_m <- df[rownames(df) %in% male_names, ]

#zero counts by gender 

zero_count_f <- sum(df_f == 0, na.rm = TRUE)
zero_count_m <- sum(df_m == 0, na.rm = TRUE)

ggplot(data = melt(df100_f), aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=variable)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(title = "Summary of selected CAZymes for females", x = "CAZymes", y = "Count") +
  theme(legend.position="none")

ggplot(data = melt(df100_m), aes(x=variable, y=value)) +
  geom_boxplot(aes(fill=variable)) + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(title = "Summary of selected CAZymes for males", x = "CAZymes", y = "Count") +
  theme(legend.position="none")
















# Calculate the count of values equal to 100
count_100 <- sum(count100)

# Subset the names where values in the third column are equal to 100
names_100 <- result_df$Ecs[count100]
print(names_100)













# Print the result
if (has_duplicates) {
  print("Duplicate row names found.")
} else {
  print("No duplicate row names.")
}










library(ggplot2)

ggplot(DAG3, aes(x = ANTHRO.AGE, y = ANTHRO.Sex, fill = ANTHRO.Sex)) + theme_bw() +
  geom_violin(color = "black", alpha = 0.5) +
  geom_boxplot(width = 0.2, color = "black", fill = "white", alpha = 0.5) +
  labs(title = "Age Distribution by Sex", x = "Sex", y = "Age") +
  scale_fill_manual(values = c("lightblue", "lightpink"), guide = FALSE) +
  facet_wrap(~ ANTHRO.Sex, nrow = 1)


ggplot(DAG3, aes(x = DAG3$ANTHRO.AGE, fill = DAG3$ANTHRO.Sex)) +
  geom_histogram(color = "black", bins = 30) +
  labs(title = "Age Distribution by Sex", x = "Age", y = "Count", fill = "Gender") + theme_bw() +
  scale_fill_manual(values = c("lightblue", "lightpink")) +
  facet_wrap(~ DAG3$ANTHRO.Sex, nrow = 1)

#prevalent/rarest via min/max of colsum 
# Convert all columns to numeric
cazyme_data <- sapply(df, as.numeric)
# Calculate the sum of reads for each CAZyme
cazyme_counts <- colSums(cazyme_data, na.rm = TRUE)
most_prevalent_cazyme <- names(cazyme_counts)[which.max(cazyme_counts)]
rarest_cazyme <- names(cazyme_counts)[which.min(cazyme_counts)]


#3.2.1.157 [85]





# Create an empty dataframe to store the results
output_df <- data.frame(Column = character(),
                        NonZeroPercentage = numeric(),
                        ZeroPercentage = numeric(),
                        stringsAsFactors = FALSE)

# Assuming your dataframe is named 'df'

non_zero_count <- sum(df[, 7] != 0)

print(paste("Number of non-zero values in the first column:", non_zero_count))














# Assuming your data frame is named 'merged_data'

# Selecting the columns for bloating symptoms and CAZymes
bloating_symptoms <- merged_df100$MED.DISEASES.Gastrointestinal.Discomfort.Bloating
cazyme_columns <- colnames(merged_df100)[36:55]

# Creating an empty data frame to store regression results
regression_results <- data.frame(CAZyme = character(), Coefficient = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Performing linear regression for each CAZyme
for (cazyme in cazyme_columns) {
  regression_model <- lm(bloating_symptoms ~ merged_df100[[cazyme]], data = DAG3)
  coefficient <- coef(regression_model)[[2]]
  p_value <- summary(regression_model)$coefficients[2, 4]
  regression_results <- rbind(regression_results, data.frame(CAZyme = cazyme, Coefficient = coefficient, p_value = p_value))
}

# Sorting the results by p-value in ascending order
regression_results <- regression_results[order(regression_results$p_value), ]

# Printing the top CAZymes associated with bloating symptoms
print(regression_results)










# Selecting the columns for bloating symptoms and CAZymes
bloating_symptoms <- merged_df100$MED.DISEASES.Gastrointestinal.Discomfort.Bloating
cazyme_columns <- colnames(merged_df100)[36:55]

# Creating an empty data frame to store regression results
regression_results <- data.frame(CAZyme = character(), Coefficient = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Performing linear regression for each CAZyme
for (cazyme in cazyme_columns) {
  regression_model <- lm(bloating_symptoms ~ merged_df100[[cazyme]], data = DAG3)
  coefficient <- coef(regression_model)[[2]]
  p_value <- summary(regression_model)$coefficients[2, 4]
  regression_results <- rbind(regression_results, data.frame(CAZyme = cazyme, Coefficient = coefficient, p_value = p_value))
}

# Sorting the results by p-value in ascending order
regression_results <- regression_results[order(regression_results$p_value), ]

# Printing the top CAZymes associated with bloating symptoms
print(regression_results)