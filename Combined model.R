# ============================================================
# 逻辑回归 + ROC 分析 (CA125, HE4, BMDscore 预测 label)
# 输入文件: combined_model.xls (列: CA125, HE4, BMDscore, label, Group)
# label: BAM / OC  (BAM -> 0, OC -> 1)
# Group: Training / Validation
# 输出文件保存至输入文件所在目录
# ============================================================

# 加载必要的包
library(readxl)
library(dplyr)
library(pROC)
library(logistf)   # 处理完全分离

# ---------- 1. 参数设置 ----------
# 输入文件路径（请根据实际情况修改，注意使用正斜杠或双反斜杠）
file_path <- "D:/daydayup/卵巢癌早诊/11课题/样本及测序数据/0922数据整理/model2/combined_model.xls"

# 列名（与 Excel 完全一致）
ca125_col <- "CA125"
he4_col   <- "HE4"
bmd_col   <- "BMDscore"
label_col <- "label"
group_col <- "Group"

# 标签映射: 设定阳性类别（OC 为 1，BAM 为 0）
positive_label <- "OC"
negative_label <- "BAM"

# 可选：若需要灵敏度不低于某个值的阈值，请设置 target_sens，否则使用最大约登指数
target_sens <- 0.95   # 固定灵敏度为 95%
# target_sens <- NULL

# ---------- 2. 读取并清洗数据 ----------
df_raw <- read_excel(file_path)

# 查看原始列名（调试用）
cat("原始列名:\n")
print(colnames(df_raw))

# 选择所需列，重命名，并转换 label 为 0/1
df <- df_raw %>%
  select(all_of(c(ca125_col, he4_col, bmd_col, label_col, group_col))) %>%
  rename(
    CA125 = all_of(ca125_col),
    HE4   = all_of(he4_col),
    BMDscore = all_of(bmd_col),
    Label = all_of(label_col),
    Group = all_of(group_col)
  ) %>%
  mutate(
    # 将 label 转为 0/1（BAM=0, OC=1）
    Y = case_when(
      Label == positive_label ~ 1,
      Label == negative_label ~ 0,
      TRUE ~ NA_real_
    ),
    # 确保 BMDscore 为数值（若读入为字符则转换）
    BMDscore = as.numeric(as.character(BMDscore)),
    # 对 CA125 和 HE4 取对数（加微小常数避免 log(0)）
    logCA125 = log(CA125 + 1e-6),
    logHE4   = log(HE4 + 1e-6)
  ) %>%
  filter(!is.na(Y)) %>%          # 剔除无法识别 label 的行
  filter(Group %in% c("Training", "Validation")) %>%  # 只保留指定分组
  na.omit()                      # 删除任何含有 NA 的行

cat("\n清洗后数据行数:", nrow(df), "\n")
cat("label 分布:\n")
print(table(df$Y, df$Label, dnn = c("Y", "Label")))

# 分割训练集和验证集
train_df <- df %>% filter(Group == "Training")
valid_df <- df %>% filter(Group == "Validation")

cat("\n训练集行数:", nrow(train_df), " 验证集行数:", nrow(valid_df), "\n")

# 检查训练集是否有足够变异
if (length(unique(train_df$Y)) < 2) {
  stop("训练集中 Y 只有一类，无法拟合逻辑回归！请检查数据划分。")
}
if (length(unique(train_df$BMDscore)) == 1) {
  warning("训练集中 BMDscore 只有一个唯一值，将从模型中剔除。")
  use_bmd <- FALSE
} else {
  use_bmd <- TRUE
}

# ---------- 3. 拟合逻辑回归 ----------
# 构造公式
if (use_bmd) {
  formula_lr <- Y ~ logCA125 + logHE4 + BMDscore
} else {
  formula_lr <- Y ~ logCA125 + logHE4
}

# 尝试标准 glm，若出现完全分离则自动使用 Firth 回归
fit <- tryCatch({
  glm(formula_lr, data = train_df, family = binomial)
}, warning = function(w) {
  if (grepl("拟合概率算出来是数值零或一|glm.fit: fitted probabilities numerically 0 or 1", w$message)) {
    message("检测到完全分离，改用 Firth 逻辑回归 (logistf)")
    return(logistf(formula_lr, data = train_df))
  } else {
    warning(w)
    return(NULL)
  }
})

if (is.null(fit)) {
  fit <- logistf(formula_lr, data = train_df)
}

# 输出模型系数和 OR
cat("\n========== 模型系数 ==========\n")
if (inherits(fit, "logistf")) {
  coefs <- coef(fit)
  ors <- exp(coefs)
  ci <- exp(confint(fit))
  results_coef <- data.frame(
    Variable = names(coefs),
    Coefficient = round(coefs, 4),
    OR = round(ors, 3),
    CI_lower = round(ci[,1], 3),
    CI_upper = round(ci[,2], 3)
  )
  print(results_coef)
} else {
  summary_fit <- summary(fit)
  print(summary_fit)
  cat("\nOR 及 95% CI:\n")
  or_ci <- exp(cbind(OR = coef(fit), confint(fit)))
  print(round(or_ci, 3))
}

# ---------- 4. 训练集 ROC 分析 ----------
# 预测概率
train_prob <- predict(fit, type = "response")
roc_train <- roc(train_df$Y, train_prob, quiet = TRUE)
auc_train <- auc(roc_train)

# 选择阈值：固定灵敏度 ≥ target_sens 且特异度最高
all_coords <- coords(roc_train, "all", ret = c("threshold", "sensitivity", "specificity"))
if (!is.null(target_sens)) {
  valid_coords <- all_coords[all_coords$sensitivity >= target_sens, ]
  if (nrow(valid_coords) == 0) {
    best_idx <- which.max(all_coords$sensitivity)
    cat("\n警告：无法达到目标灵敏度", target_sens, "，取实际最高灵敏度\n")
  } else {
    best_idx <- which.max(valid_coords$specificity)
  }
} else {
  # 最大约登指数 (Youden)
  best_idx <- which.max(all_coords$sensitivity + all_coords$specificity)
}
best_threshold <- all_coords$threshold[best_idx]
train_sens <- all_coords$sensitivity[best_idx]
train_spec <- all_coords$specificity[best_idx]

cat("\n========== 训练集最佳阈值 ==========\n")
cat("阈值:", round(best_threshold, 4), "\n")
cat("灵敏度 (固定≥95%):", round(train_sens, 4), "\n")
cat("特异度:", round(train_spec, 4), "\n")
cat("AUC:", round(auc_train, 4), "\n")

# ---------- 5. 验证集评估 ----------
valid_prob <- predict(fit, newdata = valid_df, type = "response")
valid_pred_class <- ifelse(valid_prob >= best_threshold, 1, 0)

# 四格表指标
TP <- sum(valid_pred_class == 1 & valid_df$Y == 1)
FN <- sum(valid_pred_class == 0 & valid_df$Y == 1)
TN <- sum(valid_pred_class == 0 & valid_df$Y == 0)
FP <- sum(valid_pred_class == 1 & valid_df$Y == 0)

valid_sens <- TP / (TP + FN)
valid_spec <- TN / (TN + FP)

# Wilson 置信区间
sens_ci <- prop.test(TP, TP+FN, correct = FALSE)$conf.int
spec_ci <- prop.test(TN, TN+FP, correct = FALSE)$conf.int

cat("\n========== 验证集性能 (使用训练集阈值) ==========\n")
cat("阈值:", round(best_threshold, 4), "\n")
cat("灵敏度:", round(valid_sens, 3), " (95% CI: ", round(sens_ci[1], 3), "-", round(sens_ci[2], 3), ")\n")
cat("特异度:", round(valid_spec, 3), " (95% CI: ", round(spec_ci[1], 3), "-", round(spec_ci[2], 3), ")\n")
cat("\n四格表:\n")
cat("           真实阳性  真实阴性\n")
cat("预测阳性      ", TP, "      ", FP, "\n")
cat("预测阴性      ", FN, "      ", TN, "\n")

# ---------- 6. 输出文件 ----------
out_dir <- dirname(file_path)

# 6.1 模型摘要 (文本)
sink(file.path(out_dir, "model_summary.txt"))
cat("===== 逻辑回归模型摘要 =====\n\n")
if (inherits(fit, "logistf")) {
  cat("Firth 逻辑回归结果 (处理完全分离)\n")
  print(results_coef)
} else {
  print(summary_fit)
  cat("\nOR 及 95% CI:\n")
  print(or_ci)
}
cat("\n===== 训练集 ROC 信息 =====\n")
cat("AUC:", round(auc_train, 4), "\n")
cat("固定灵敏度 ≥ 0.95 时选择的阈值:", round(best_threshold, 4), "\n")
cat("训练集灵敏度:", round(train_sens, 4), "\n")
cat("训练集特异度:", round(train_spec, 4), "\n")
cat("\n===== 验证集性能 =====\n")
cat("灵敏度:", round(valid_sens, 3), " (95% CI: ", round(sens_ci[1], 3), "-", round(sens_ci[2], 3), ")\n")
cat("特异度:", round(valid_spec, 3), " (95% CI: ", round(spec_ci[1], 3), "-", round(spec_ci[2], 3), ")\n")
cat("\n===== 特异度总结 =====\n")
cat("固定灵敏度为 0.95 时，训练集特异度为", round(train_spec, 4), "，验证集特异度为", round(valid_spec, 4), "\n")
sink()

# 6.2 ROC 曲线坐标
roc_coords <- all_coords
roc_coords$FPR <- 1 - roc_coords$specificity
roc_coords <- roc_coords[order(roc_coords$FPR), ]
write.csv(roc_coords, file.path(out_dir, "ROC_coordinates.csv"), row.names = FALSE)

# 6.3 训练集预测概率
train_out <- data.frame(True_Label = train_df$Y, Pred_Prob = train_prob)
write.csv(train_out, file.path(out_dir, "Training_predictions.csv"), row.names = FALSE)

# 6.4 验证集预测概率及类别
valid_out <- data.frame(True_Label = valid_df$Y, Pred_Prob = valid_prob, Pred_Class = valid_pred_class)
write.csv(valid_out, file.path(out_dir, "Validation_predictions.csv"), row.names = FALSE)

# 6.5 ROC 曲线图 (训练集)
png(file.path(out_dir, "ROC_Training.png"), width = 800, height = 600)
plot(roc_train, col = "blue", lwd = 2, 
     main = paste("训练集 ROC 曲线 (AUC =", round(auc_train, 3), ")"))
points(1 - train_spec, train_sens, col = "red", pch = 19, cex = 1.5)
legend("bottomright", 
       legend = c(paste("阈值 =", round(best_threshold, 3)),
                  paste("灵敏度 =", round(train_sens, 3)),
                  paste("特异度 =", round(train_spec, 3))),
       col = "red", pch = 19, bty = "n")
dev.off()

cat("\n所有输出文件已保存至:", out_dir, "\n")
cat("生成的文件:\n")
cat("  - model_summary.txt\n  - ROC_coordinates.csv\n  - Training_predictions.csv\n  - Validation_predictions.csv\n  - ROC_Training.png\n")

# ---------- 7. 补充特异度输出（直接显示在控制台） ----------
cat("\n========================================\n")
cat("固定灵敏度为 95% 时，最终结果如下：\n")
cat("  所选阈值 =", round(best_threshold, 4), "\n")
cat("  训练集特异度 =", round(train_spec, 4), "\n")
cat("  验证集特异度 =", round(valid_spec, 4), "\n")
cat("========================================\n")

# ---------- 8. 计算训练集预测类别及 95% CI ----------
# 使用与验证集相同的阈值 best_threshold 对训练集进行分类
train_pred_class <- ifelse(train_prob >= best_threshold, 1, 0)

# 训练集四格表
TP_train <- sum(train_pred_class == 1 & train_df$Y == 1)
FN_train <- sum(train_pred_class == 0 & train_df$Y == 1)
TN_train <- sum(train_pred_class == 0 & train_df$Y == 0)
FP_train <- sum(train_pred_class == 1 & train_df$Y == 0)

# 训练集灵敏度和特异度的 Wilson 95% 置信区间
sens_ci_train <- prop.test(TP_train, TP_train + FN_train, correct = FALSE)$conf.int
spec_ci_train <- prop.test(TN_train, TN_train + FP_train, correct = FALSE)$conf.int

cat("\n========== 训练集性能 (使用选定阈值) ==========\n")
cat("阈值:", round(best_threshold, 4), "\n")
cat("灵敏度:", round(train_sens, 4), 
    " (95% CI: ", round(sens_ci_train[1], 4), " - ", round(sens_ci_train[2], 4), ")\n")
cat("特异度:", round(train_spec, 4), 
    " (95% CI: ", round(spec_ci_train[1], 4), " - ", round(spec_ci_train[2], 4), ")\n")
cat("四格表: TP=", TP_train, " FN=", FN_train, " TN=", TN_train, " FP=", FP_train, "\n")

# 可选：将训练集 CI 追加到 model_summary.txt 中
sink(file.path(out_dir, "model_summary.txt"), append = TRUE)
cat("\n===== 训练集性能 (使用选定阈值) =====\n")
cat("阈值:", round(best_threshold, 4), "\n")
cat("灵敏度:", round(train_sens, 4), 
    " (95% CI: ", round(sens_ci_train[1], 4), " - ", round(sens_ci_train[2], 4), ")\n")
cat("特异度:", round(train_spec, 4), 
    " (95% CI: ", round(spec_ci_train[1], 4), " - ", round(spec_ci_train[2], 4), ")\n")
cat("四格表: TP=", TP_train, " FN=", FN_train, " TN=", TN_train, " FP=", FP_train, "\n")
sink()

cat("\n训练集 95% 置信区间已追加至 model_summary.txt\n")