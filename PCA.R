# 加载必要的包
if (!require(readxl)) install.packages("readxl")
if (!require(factoextra)) install.packages("factoextra")
if (!require(cluster)) install.packages("cluster")
library(readxl)
library(factoextra)
library(cluster)

# 读取数据
file_path <- "D:/OC/PCA_eva.xls"
data_raw <- read_excel(file_path, col_names = TRUE)

if (is.character(data_raw[[1]]) || is.factor(data_raw[[1]])) {
  rownames(data_raw) <- data_raw[[1]]
  data_raw <- data_raw[, -1]
}

scores <- as.data.frame(data_raw[, 1:2])
colnames(scores) <- c("PC1", "PC2")
scores[] <- lapply(scores, as.numeric)

set.seed(123) 
gap_result <- clusGap(scores, FUN = kmeans, nstart = 25, K.max = 10, B = 50)


print(gap_result, method = "firstmax")

fviz_gap_stat(gap_result) + 
  labs(title = "Gap Statistic for Number of Clusters (PC1 & PC2)")