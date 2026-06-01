# 加载包
library(readxl)
library(dplyr)
library(pROC)
library(logistf)   # Firth 逻辑回归（处理完全分离）

# ---------- 1. 读取数据 ----------
file_path <- r"(D:\OC)"
df_raw <- read_excel(file_path)

print(colnames(df_raw))

label_col <- "label"         
group_col <- "Group"          
ca125_col <- "CA125"          
he4_col <- "HE4"              
# ------------------------------------------------

# 数据预处理
df <- df_raw %>%
  select(all_of(c(label_col, group_col, ca125_col, he4_col))) %>%
  rename(
    Y = all_of(label_col),
    Group = all_of(group_col),
    CA125 = all_of(ca125_col),
    HE4 = all_of(he4_col)
  ) %>%
  filter(Y %in% c(0, 1)) %>%
  mutate(
    logCA125 = log(CA125 + 1e-6),
    logHE4 = log(HE4 + 1e-6)
  ) %>%
  na.omit()


train_df <- df %>% filter(Group == "Training")
valid_df <- df %>% filter(Group == "Validation")


fit <- tryCatch({
  glm(Y ~ logCA125 + logHE4, data = train_df, family = binomial)
}, warning = function(w) {
  if (grepl("拟合概率算出来是数值零或一", w$message)) {
    message("检测到完全分离，改用 Firth 逻辑回归")
    return(logistf(Y ~ logCA125 + logHE4, data = train_df))
  } else {
    warning(w)
    return(NULL)
  }
})

if (is.null(fit)) {
  fit <- logistf(Y ~ logCA125 + logHE4, data = train_df)
}

# 输出模型摘要
cat("\n========== 模型系数 ==========\n")
if (inherits(fit, "logistf")) {
  coefs <- coef(fit)
  ors <- exp(coefs)
  ci <- exp(confint(fit))
  for (i in 1:length(coefs)) {
    cat(names(coefs)[i], ": 系数 =", round(coefs[i], 4), 
        "  OR =", round(ors[i], 3), 
        "  95% CI =", round(ci[i,1], 3), "-", round(ci[i,2], 3), "\n")
  }
} else {
  summary(fit)
  cat("OR 及 95% CI:\n")
  print(exp(cbind(OR = coef(fit), confint(fit))))
}

# ---------- 3. 训练集预测概率 ----------
if (inherits(fit, "logistf")) {
  train_prob <- predict(fit, type = "response")
} else {
  train_prob <- predict(fit, type = "response")
}
roc_train <- roc(train_df$Y, train_prob, quiet = TRUE)
auc_train <- auc(roc_train)


target_sens <- 0.95
all_coords <- coords(roc_train, "all", ret = c("threshold", "sensitivity", "specificity"))
valid_coords <- all_coords[all_coords$sensitivity >= target_sens, ]

if (nrow(valid_coords) == 0) {
  
  best_idx <- which.max(all_coords$sensitivity)
  cutoff_threshold <- all_coords$threshold[best_idx]
  train_sens_at_cutoff <- all_coords$sensitivity[best_idx]
  train_spec_at_cutoff <- all_coords$specificity[best_idx]
  cat("警告：无法达到目标灵敏度", target_sens, "，取实际最高灵敏度", round(train_sens_at_cutoff, 4), "\n")
} else {
  best_idx <- which.max(valid_coords$specificity)
  cutoff_threshold <- valid_coords$threshold[best_idx]
  train_sens_at_cutoff <- valid_coords$sensitivity[best_idx]
  train_spec_at_cutoff <- valid_coords$specificity[best_idx]
}

cat("\n基于训练集 (灵敏度 ≥ 95%) 选定的阈值:", round(cutoff_threshold, 4), "\n")
cat("实际训练集灵敏度:", round(train_sens_at_cutoff, 4), 
    " 特异度:", round(train_spec_at_cutoff, 4), "\n")

if (inherits(fit, "logistf")) {
  valid_prob <- predict(fit, newdata = valid_df, type = "response")
} else {
  valid_prob <- predict(fit, newdata = valid_df, type = "response")
}
valid_pred_class <- ifelse(valid_prob >= cutoff_threshold, 1, 0)

# 验证集四格表
TP <- sum(valid_pred_class == 1 & valid_df$Y == 1)
FN <- sum(valid_pred_class == 0 & valid_df$Y == 1)
TN <- sum(valid_pred_class == 0 & valid_df$Y == 0)
FP <- sum(valid_pred_class == 1 & valid_df$Y == 0)

valid_sens <- TP / (TP + FN)
valid_spec <- TN / (TN + FP)

# 95% CI (Wilson)
sens_ci <- prop.test(TP, TP+FN, correct = FALSE)$conf.int
spec_ci <- prop.test(TN, TN+FP, correct = FALSE)$conf.int

cat("\n========== 验证集性能 (使用训练集阈值) ==========\n")
cat("阈值:", round(cutoff_threshold, 4), "\n")
cat("灵敏度:", round(valid_sens, 3), " (95% CI: ", round(sens_ci[1], 3), "-", round(sens_ci[2], 3), ")\n")
cat("特异度:", round(valid_spec, 3), " (95% CI: ", round(spec_ci[1], 3), "-", round(spec_ci[2], 3), ")\n")
cat("\n四格表:\n")
cat("           真实阳性  真实阴性\n")
cat("预测阳性      ", TP, "      ", FP, "\n")
cat("预测阴性      ", FN, "      ", TN, "\n")

# ---------- 6. 输出 ROC 曲线坐标及预测结果 ----------
out_dir <- dirname(file_path)

# ROC 坐标
roc_coords <- all_coords
roc_coords$FPR <- 1 - roc_coords$specificity
roc_coords <- roc_coords[order(roc_coords$FPR), ]
write.csv(roc_coords, file.path(out_dir, "ROC_coordinates_Training_CA125_HE4_sens95.csv"), row.names = FALSE)

# 预测概率
train_pred_df <- data.frame(True_Label = train_df$Y, Pred_Prob = train_prob)
write.csv(train_pred_df, file.path(out_dir, "Training_predictions_sens95.csv"), row.names = FALSE)

valid_pred_df <- data.frame(True_Label = valid_df$Y, Pred_Prob = valid_prob, Pred_Class = valid_pred_class)
write.csv(valid_pred_df, file.path(out_dir, "Validation_predictions_sens95.csv"), row.names = FALSE)

cat("\n输出文件已保存至:", out_dir, "\n")

# ---------- 7. 绘制训练集 ROC 曲线并标出 cutoff 点 ----------
png(file.path(out_dir, "ROC_Training_CA125_HE4_sens95.png"), width = 800, height = 600)
plot(roc_train, col = "blue", lwd = 2, 
     main = paste("训练集 ROC (AUC =", round(auc_train, 3), ")"))
points(1 - train_spec_at_cutoff, train_sens_at_cutoff, col = "red", pch = 19, cex = 1.5)
legend("bottomright", legend = c(paste("阈值 =", round(cutoff_threshold, 3)),
                                 paste("灵敏度 =", round(train_sens_at_cutoff, 3)),
                                 paste("特异度 =", round(train_spec_at_cutoff, 3))),
       col = "red", pch = 19, bty = "n")
dev.off()
cat("ROC 曲线图已保存\n")


# 读取 BOT 数据
bot_path <- r"(D:\OC\BOT.xlsx)"
if (!file.exists(bot_path)) {
  stop("BOT.xlsx 文件不存在，请检查路径")
}
bot_df_raw <- read_excel(bot_path)


print(colnames(bot_df_raw))
bot_df <- bot_df_raw %>%
  select(label, CA125, HE4) %>%
  rename(Y = label) %>%
  filter(Y %in% c(0, 1)) %>%
  mutate(
    logCA125 = log(CA125 + 1e-6),
    logHE4 = log(HE4 + 1e-6)
  ) %>%
  na.omit()

cat("BOT 数据集样本量:", nrow(bot_df), 
    "  阳性(1):", sum(bot_df$Y), 
    " 阴性(0):", nrow(bot_df) - sum(bot_df$Y), "\n")

# 使用已训练的模型进行预测
if (inherits(fit, "logistf")) {
  bot_prob <- predict(fit, newdata = bot_df, type = "response")
} else {
  bot_prob <- predict(fit, newdata = bot_df, type = "response")
}


bot_pred_class <- ifelse(bot_prob >= cutoff_threshold, 1, 0)


bot_pred_df <- data.frame(
  True_Label = bot_df$Y,
  Pred_Prob = bot_prob,
  Pred_Class = bot_pred_class
)
write.csv(bot_pred_df, file.path(out_dir, "BOT_predictions_sens95.csv"), row.names = FALSE)
cat("\nBOT 预测结果已保存至:", file.path(out_dir, "BOT_predictions_sens95.csv"), "\n")