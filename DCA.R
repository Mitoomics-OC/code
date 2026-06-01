if (!require("rmda")) install.packages("rmda")
library(rmda)
library(readxl)

setwd("D:\\OC")
wd <- "D:\\OC"
file_path <- paste0(wd, "/DCA1.xlsx")  # 替换为实际文件名
if (!file.exists(file_path)) stop("文件不存在！")
data <- readxl::read_excel(file_path)

prevalence <- 0.5

## 决策曲线分析（病例对照设计）
dca_result <- decision_curve(
  outcome ~ pred2,
  data = data,
  thresholds = seq(0, 1, by = 0.01),
  study.design = "case-control",
  population.prevalence = prevalence
)


#绘图
while (!is.null(dev.list())) dev.off()

# 绘制决策曲线
plot_decision_curve(
  dca_result,
  curve.names = c('Model2'),         # 曲线名字
  cost.benefit.axis = FALSE,        # 是否显示cost-benefit轴
  col = c('red','blue','green'),   # 曲线颜色
  confidence.intervals = FALSE,     # 是否画置信区间
  standardize = FALSE
)


# 导出净获益数据（每个阈值下的净获益，适合绘图和进一步分析）
write.csv(dca_result[["derived.data"]], file = "dca_derived_data_0527_0.5_model_combined.csv", row.names = FALSE)
