#!/usr/bin/env Rscript
# =============================================================================
# generate_figures.R
# 마이크로바이옴 CRC 예측 논문용 그래프 생성 스크립트
#
# 전제 조건:
#   run_analysis.R 이 완료되어 C:/Users/minta/ICBBS/r_results/ 에 결과 저장됨
#
# 출력:
#   r_results/figures/ 폴더에 PNG 파일 저장
#   - fig1_roc_dev.png        Figure 1: ROC (개발 코호트)
#   - fig2_roc_ext.png        Figure 2: ROC (외부 검증)
#   - fig3_perf_barplot.png   Figure 3: Performance barplot
#   - fig4_feature_importance Figure 4: XGBoost feature importance
#   - fig5_generalizability   Figure 5: AUC generalizability
#   - fig_cal.png             Figure 6: Calibration
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(jsonlite)
})

RESULTS_DIR <- "C:/Users/minta/ICBBS/r_results"
FIG_DIR     <- file.path(RESULTS_DIR, "figures")
dir.create(FIG_DIR, showWarnings=FALSE, recursive=TRUE)

DPI    <- 300
WIDTH  <- 7
HEIGHT <- 5.5
MODEL_COLORS <- c("LR"="#4472C4", "RF"="#ED7D31", "XGB"="#70AD47")
MODEL_LABELS <- c("LR"="Logistic Regression", "RF"="Random Forest", "XGB"="XGBoost")

cat("=== 그래프 생성 시작 ===\n")

dev_path <- file.path(RESULTS_DIR, "metrics_dev_test.tsv")
ext_path <- file.path(RESULTS_DIR, "metrics_external.tsv")
ece_path <- file.path(RESULTS_DIR, "ece_summary.tsv")
roc_path <- file.path(RESULTS_DIR, "roc_curves.rds")
if (!file.exists(dev_path)) stop("[ERROR] metrics_dev_test.tsv 없음. run_analysis.R 먼저 실행.")
res_dev <- read.delim(dev_path, stringsAsFactors=FALSE)
res_ext <- read.delim(ext_path, stringsAsFactors=FALSE)
cat("[OK] 결과 파일 로드 완료\n")

# Figure 1 & 2: ROC Curves
if (file.exists(roc_path)) {
  roc_list <- readRDS(roc_path)
  make_roc_df <- function(roc_obj, model_name) {
    data.frame(FPR=1-roc_obj$specificities, TPR=roc_obj$sensitivities, Model=model_name)
  }
  roc_df_dev <- bind_rows(make_roc_df(roc_list$dev$LR,"LR"), make_roc_df(roc_list$dev$RF,"RF"), make_roc_df(roc_list$dev$XGB,"XGB"))
  roc_df_ext <- bind_rows(make_roc_df(roc_list$ext$LR,"LR"), make_roc_df(roc_list$ext$RF,"RF"), make_roc_df(roc_list$ext$XGB,"XGB"))

  auc_dev <- setNames(round(res_dev$AUC,3), res_dev$Model)
  dev_labs <- setNames(sprintf("%s (AUC=%.3f, 95%%CI:%.3f–%.3f)", MODEL_LABELS[res_dev$Model], res_dev$AUC, res_dev$AUC_CI_L, res_dev$AUC_CI_U), res_dev$Model)
  g1 <- ggplot(roc_df_dev, aes(x=FPR,y=TPR,color=Model,linetype=Model)) +
    geom_line(linewidth=1) + geom_abline(intercept=0,slope=1,linetype="dashed",color="grey50") +
    scale_color_manual(values=MODEL_COLORS, labels=dev_labs) +
    scale_linetype_manual(values=c(LR="solid",RF="dashed",XGB="solid"), labels=dev_labs) +
    coord_equal(xlim=c(0,1),ylim=c(0,1)) +
    labs(title="Figure 1. ROC Curves – Development Cohort (ZellerG_2014)",
         subtitle="95% CI by 2,000 bootstrap iterations",
         x="False Positive Rate", y="True Positive Rate", color="Model", linetype="Model") +
    theme_classic(base_size=12) +
    theme(legend.position=c(0.65,0.25), legend.background=element_rect(fill="white",color="grey80"), legend.text=element_text(size=8.5))
  ggsave(file.path(FIG_DIR,"fig1_roc_dev.png"), g1, width=WIDTH, height=HEIGHT, dpi=DPI)
  cat("[OK] fig1_roc_dev.png\n")

  ext_labs <- setNames(sprintf("%s (AUC=%.3f, 95%%CI:%.3f–%.3f)", MODEL_LABELS[res_ext$Model], res_ext$AUC, res_ext$AUC_CI_L, res_ext$AUC_CI_U), res_ext$Model)
  g2 <- ggplot(roc_df_ext, aes(x=FPR,y=TPR,color=Model,linetype=Model)) +
    geom_line(linewidth=1) + geom_abline(intercept=0,slope=1,linetype="dashed",color="grey50") +
    scale_color_manual(values=MODEL_COLORS, labels=ext_labs) +
    scale_linetype_manual(values=c(LR="solid",RF="dashed",XGB="solid"), labels=ext_labs) +
    coord_equal(xlim=c(0,1),ylim=c(0,1)) +
    labs(title="Figure 2. ROC Curves – External Validation (YachidaS_2019)",
         subtitle="Models trained on ZellerG_2014, applied without retraining",
         x="False Positive Rate", y="True Positive Rate", color="Model", linetype="Model") +
    theme_classic(base_size=12) +
    theme(legend.position=c(0.65,0.25), legend.background=element_rect(fill="white",color="grey80"), legend.text=element_text(size=8.5))
  ggsave(file.path(FIG_DIR,"fig2_roc_ext.png"), g2, width=WIDTH, height=HEIGHT, dpi=DPI)
  cat("[OK] fig2_roc_ext.png\n")
} else {
  cat("[WARN] roc_curves.rds 없음 → AUC bar chart 대체\n")
  auc_df <- bind_rows(res_dev%>%mutate(Cohort="Development"), res_ext%>%mutate(Cohort="External")) %>%
    select(Model,Cohort,AUC,AUC_CI_L,AUC_CI_U)
  for (coh in c("Development","External")) {
    fn <- if(coh=="Development") "fig1_roc_dev.png" else "fig2_roc_ext.png"
    g <- ggplot(auc_df%>%filter(Cohort==coh), aes(x=Model,y=AUC,fill=Model)) +
      geom_col(width=0.6) + geom_errorbar(aes(ymin=AUC_CI_L,ymax=AUC_CI_U),width=0.15,linewidth=0.8) +
      geom_text(aes(label=sprintf("%.3f",AUC)),vjust=-1.0,size=3.5) +
      scale_fill_manual(values=MODEL_COLORS) + scale_y_continuous(limits=c(0,1.05)) +
      labs(title=sprintf("AUC – %s Cohort",coh), x="Model", y="AUC") +
      theme_classic(base_size=12) + theme(legend.position="none")
    ggsave(file.path(FIG_DIR,fn), g, width=6, height=5, dpi=DPI)
    cat("[OK]",fn,"(AUC bar)\n")
  }
}

# Figure 3: Performance Barplot
perf_df <- bind_rows(res_dev%>%mutate(Cohort="Development (n=87)"), res_ext%>%mutate(Cohort="External (n=212)")) %>%
  select(Model,Cohort,AUC,Recall,F1) %>%
  pivot_longer(cols=c(AUC,Recall,F1), names_to="Metric", values_to="Value") %>%
  mutate(Model=factor(Model,levels=c("LR","RF","XGB")), Metric=factor(Metric,levels=c("AUC","Recall","F1")),
         Cohort=factor(Cohort,levels=c("Development (n=87)","External (n=212)")))
g3 <- ggplot(perf_df, aes(x=Model,y=Value,fill=Model)) +
  geom_col(position="dodge",width=0.7) +
  geom_text(aes(label=sprintf("%.3f",Value)),position=position_dodge(0.7),vjust=-0.5,size=3) +
  scale_fill_manual(values=MODEL_COLORS, labels=MODEL_LABELS) +
  scale_y_continuous(limits=c(0,1.10),breaks=seq(0,1,0.2)) +
  facet_grid(Cohort~Metric) +
  labs(title="Figure 3. Performance Comparison", x="Model", y="Score", fill="Model") +
  theme_classic(base_size=11) +
  theme(strip.background=element_rect(fill="#1F4973",color=NA), strip.text=element_text(color="white",face="bold"), legend.position="bottom")
ggsave(file.path(FIG_DIR,"fig3_perf_barplot.png"), g3, width=9, height=6.5, dpi=DPI)
cat("[OK] fig3_perf_barplot.png\n")

# Figure 4: Feature Importance
imp_path <- file.path(RESULTS_DIR,"xgb_importance.rds")
if (file.exists(imp_path)) {
  imp_df <- readRDS(imp_path)
  top20  <- imp_df %>% arrange(desc(Gain)) %>% head(20) %>%
    mutate(Feature=factor(Feature,levels=rev(Feature)),
           Highlight=grepl("Fusobacterium|Peptostreptococcus|Parvimonas|Porphyromonas",Feature,ignore.case=TRUE))
  g4 <- ggplot(top20, aes(x=Gain,y=Feature,fill=Highlight)) +
    geom_col() + scale_fill_manual(values=c("FALSE"="#4472C4","TRUE"="#C0504D"), labels=c("FALSE"="Other","TRUE"="CRC-associated")) +
    labs(title="Figure 4. XGBoost Feature Importance (Top-20)", subtitle="Gain = reduction in loss", x="Importance (Gain)", y="Taxa", fill="Category") +
    theme_classic(base_size=11) + theme(axis.text.y=element_text(size=9), legend.position="bottom")
  ggsave(file.path(FIG_DIR,"fig4_feature_importance.png"), g4, width=WIDTH, height=7, dpi=DPI)
  cat("[OK] fig4_feature_importance.png\n")
} else cat("[WARN] xgb_importance.rds 없음 → Figure 4 건너뜀\n")

# Figure 5: Generalizability
gen_df <- bind_rows(res_dev%>%mutate(Cohort="Development"), res_ext%>%mutate(Cohort="External")) %>%
  select(Model,Cohort,AUC,AUC_CI_L,AUC_CI_U) %>%
  mutate(Model=factor(Model,levels=c("LR","RF","XGB")), Cohort=factor(Cohort,levels=c("Development","External")))
g5 <- ggplot(gen_df, aes(x=Model,y=AUC,color=Cohort,group=Cohort)) +
  geom_point(size=4,position=position_dodge(0.4)) + geom_line(position=position_dodge(0.4),linewidth=1) +
  geom_errorbar(aes(ymin=AUC_CI_L,ymax=AUC_CI_U),width=0.15,linewidth=0.8,position=position_dodge(0.4)) +
  geom_text(aes(label=sprintf("%.3f",AUC)),position=position_dodge(0.4),vjust=-1.2,size=3.2) +
  scale_color_manual(values=c("Development"="#1F4973","External"="#C0504D")) +
  scale_y_continuous(limits=c(0.60,1.05),breaks=seq(0.6,1.0,0.1)) +
  labs(title="Figure 5. AUC Generalizability",subtitle="Error bars: 95% bootstrap CI",x="Model",y="AUC",color="Cohort") +
  theme_classic(base_size=12) + theme(legend.position="top")
ggsave(file.path(FIG_DIR,"fig5_generalizability.png"), g5, width=6.5, height=5.5, dpi=DPI)
cat("[OK] fig5_generalizability.png\n")

# Figure 6: Calibration
cal_rds <- file.path(RESULTS_DIR,"cal_data.rds")
if (file.exists(cal_rds)) {
  cat("[OK] cal_data.rds 존재 → reliability curve 생성 (implementation 필요시 확장)\n")
} else if (file.exists(ece_path)) {
  ece_df <- read.delim(ece_path, stringsAsFactors=FALSE) %>%
    mutate(Model=factor(Model,levels=c("LR","RF","XGB")),
           Cohort=factor(Cohort,levels=c("DEV_Test","External"),labels=c("Development","External")))
  g_cal <- ggplot(ece_df, aes(x=Model,y=ECE,fill=Model)) +
    geom_col(width=0.6) + geom_text(aes(label=sprintf("%.3f",ECE)),vjust=-0.4,size=3.5) +
    scale_fill_manual(values=MODEL_COLORS) + scale_y_continuous(limits=c(0,max(ece_df$ECE)*1.3)) +
    facet_wrap(~Cohort) +
    labs(title="Figure 6. Expected Calibration Error (ECE)",subtitle="Lower ECE = better calibration",x="Model",y="ECE") +
    theme_classic(base_size=11) +
    theme(strip.background=element_rect(fill="#1F4973"), strip.text=element_text(color="white",face="bold"), legend.position="none")
  ggsave(file.path(FIG_DIR,"fig_cal.png"), g_cal, width=8, height=5, dpi=DPI)
  cat("[OK] fig_cal.png\n")
}

cat("\n=== 그래프 생성 완료 ===\n")
cat("저장 경로:", FIG_DIR, "\n")
figs <- list.files(FIG_DIR, pattern="\\.png$")
for (f in figs) cat(" -", f, "\n")
cat("\n다음: python add_figures.py\n")
