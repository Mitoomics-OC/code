# 加载包
library(readxl)
library(pROC)

# 文件路径
file_path <- "D:/BMD.xlsx"
output_dir <- dirname(file_path)

# 读取数据
data <- read_excel(file_path)

# 预处理
data$label <- factor(data$label, levels = c("BAM", "OC"))
data$ln_CA125 <- log(pmax(data$CA125, 0.01))
data$ln_HE4   <- log(pmax(data$HE4, 0.01))

# 划分训练集和验证集
train_data <- data[data$Group == "Training", ]
valid_data <- data[data$Group == "Validation", ]

# 拟合逻辑回归模型
lr_model <- glm(label ~ ln_CA125 + ln_HE4 + BMDscore, 
                data = train_data, family = binomial)

# 预测概率
train_prob <- predict(lr_model, type = "response")
valid_prob <- predict(lr_model, newdata = valid_data, type = "response")

# 固定训练集灵敏度为95%，确定阈值
roc_obj <- roc(train_data$label, train_prob, levels = c("BAM", "OC"), direction = "<")
target_sens <- 0.95
idx <- which.min(abs(roc_obj$sensitivities - target_sens))
threshold <- roc_obj$thresholds[idx]

# 根据阈值得到预测分类
train_pred <- ifelse(train_prob >= threshold, "OC", "BAM")
valid_pred <- ifelse(valid_prob >= threshold, "OC", "BAM")

# 将预测值添加回原数据框
data$pred_prob <- NA
data$pred_prob[data$Group == "Training"] <- train_prob
data$pred_prob[data$Group == "Validation"] <- valid_prob
data$pred_class <- NA
data$pred_class[data$Group == "Training"] <- train_pred
data$pred_class[data$Group == "Validation"] <- valid_pred

# 保存完整预测结果（CSV格式，便于查看）
output_csv <- file.path(output_dir, "predictions_CA125_HE4_BMDscore.csv")
write.csv(data[, c("state", "CA125", "HE4", "BMDscore", "label", "Group", "pred_prob", "pred_class")], 
          output_csv, row.names = FALSE)

# 控制台输出前10行
cat("预测值示例（前10行）：\n")
print(head(data[, c("label", "Group", "pred_prob", "pred_class")], 10))

cat("\n阈值（训练集灵敏度95%）：", threshold)
cat("\n预测结果已保存至：", output_csv)