# ============================================================
# Microbiome CRC Prediction — 완전 재현 가능 R 파이프라인
# ============================================================
# 실행 전 준비사항:
#   1. R 4.1 이상 설치 확인
#   2. 이 스크립트를 RStudio에서 열고 [Source] 버튼 클릭
#   3. 최초 실행 시 패키지 자동 설치 (인터넷 필요, 약 10-20분)
#   4. 결과는 C:/Users/minta/ICBBS/r_results/ 에 자동 저장
# ============================================================

cat("========================================\n")
cat("  Microbiome CRC Prediction Pipeline\n")
cat("========================================\n\n")

# STEP 0: 패키지 자동 설치
cran_pkgs <- c("data.table","dplyr","tidyr","vegan","pROC","PRROC",
               "randomForest","xgboost","ggplot2","patchwork","scales","jsonlite")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) install.packages(pkg, repos="https://cloud.r-project.org")
}
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager", repos="https://cloud.r-project.org")
bioc_pkgs <- c("curatedMetagenomicData","TreeSummarizedExperiment")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly=TRUE)) BiocManager::install(pkg, ask=FALSE, update=FALSE)
}
suppressPackageStartupMessages({
  library(curatedMetagenomicData); library(TreeSummarizedExperiment)
  library(data.table); library(dplyr); library(tidyr); library(vegan)
  library(pROC); library(PRROC); library(randomForest); library(xgboost)
  library(ggplot2); library(patchwork); library(scales); library(jsonlite)
})
cat("  All packages loaded.\n\n")

# STEP 1: 설정값
RARE_DEPTH  <- 10000
PREV_CUT    <- 0.05
N_BINS_CAL  <- 10
BOOT_N      <- 2000
TRAIN_RATIO <- 0.70
POS_LABEL   <- "CRC"
SEED        <- 42
OUTDIR <- "C:/Users/minta/ICBBS/r_results"
dir.create(OUTDIR, recursive=TRUE, showWarnings=FALSE)
cat(sprintf("  Output directory: %s\n\n", OUTDIR))

# STEP 2: 데이터 로드
cat("[2] Loading public datasets from curatedMetagenomicData...\n")
dev_se_all <- curatedMetagenomicData("ZellerG_2014.relative_abundance", dryrun=FALSE, counts=TRUE)[[1]]
dev_keep <- colData(dev_se_all)$study_condition %in% c("CRC","control")
dev_se   <- dev_se_all[, dev_keep]
dev_counts <- t(assay(dev_se))
dev_meta   <- as.data.frame(colData(dev_se))
dev_y      <- ifelse(dev_meta$study_condition == POS_LABEL, 1L, 0L)
names(dev_y) <- rownames(dev_meta)
cat(sprintf("  DEV raw samples: %d  →  CRC+control 필터 후: %d\n", ncol(dev_se_all), nrow(dev_counts)))
cat(sprintf("  DEV (ZellerG_2014): %d samples (%d CRC, %d Control), %d taxa\n",
            nrow(dev_counts), sum(dev_y==1), sum(dev_y==0), ncol(dev_counts)))

ext_se_all <- curatedMetagenomicData("YachidaS_2019.relative_abundance", dryrun=FALSE, counts=TRUE)[[1]]
ext_keep <- colData(ext_se_all)$study_condition %in% c("CRC","control")
ext_se   <- ext_se_all[, ext_keep]
ext_counts <- t(assay(ext_se))
ext_meta   <- as.data.frame(colData(ext_se))
ext_y      <- ifelse(ext_meta$study_condition == POS_LABEL, 1L, 0L)
names(ext_y) <- rownames(ext_meta)
cat(sprintf("  EXT raw samples: %d  →  CRC+control 필터 후: %d\n", ncol(ext_se_all), nrow(ext_counts)))
cat(sprintf("  EXT (YachidaS_2019): %d samples (%d CRC, %d Control), %d taxa\n\n",
            nrow(ext_counts), sum(ext_y==1), sum(ext_y==0), ncol(ext_counts)))

cohort_info <- list(
  dev=list(n_total=nrow(dev_counts), n_crc=sum(dev_y==1), n_ctrl=sum(dev_y==0), n_taxa=ncol(dev_counts)),
  ext=list(n_total=nrow(ext_counts), n_crc=sum(ext_y==1), n_ctrl=sum(ext_y==0), n_taxa=ncol(ext_counts))
)

# STEP 3: 전처리 함수
cat("[3] Preprocessing...\n")
safe_rarefy <- function(counts_mat, depth) {
  rs <- rowSums(counts_mat); keep <- rs >= depth
  if (!all(keep)) message(sprintf("  [WARN] %d samples dropped (reads < %d)", sum(!keep), depth))
  vegan::rrarefy(counts_mat[keep,,drop=FALSE], sample=depth)
}
prevalence_filter <- function(counts_mat, cut) {
  prev <- colSums(counts_mat > 0) / nrow(counts_mat)
  keep <- names(prev)[prev >= cut]
  list(filtered=counts_mat[,keep,drop=FALSE], kept_taxa=keep)
}
clr_transform <- function(counts_mat, pseudocount=1) {
  X <- counts_mat + pseudocount; log_X <- log(X)
  sweep(log_X, 1, rowMeans(log_X), "-")
}

# STEP 4: 전처리 파이프라인
set.seed(SEED)
dev_rare  <- safe_rarefy(dev_counts, RARE_DEPTH)
dev_y2    <- dev_y[rownames(dev_rare)]
pf        <- prevalence_filter(dev_rare, PREV_CUT)
dev_filt  <- pf$filtered; keep_taxa <- pf$kept_taxa
dev_clr   <- clr_transform(dev_filt)
cat(sprintf("  DEV: %d samples, %d taxa (after rarefaction + prevalence filter)\n", nrow(dev_clr), ncol(dev_clr)))

ext_rare  <- safe_rarefy(ext_counts, RARE_DEPTH)
ext_y2    <- ext_y[rownames(ext_rare)]
miss_taxa <- setdiff(keep_taxa, colnames(ext_rare))
if (length(miss_taxa) > 0) {
  pad <- matrix(0L, nrow=nrow(ext_rare), ncol=length(miss_taxa), dimnames=list(rownames(ext_rare), miss_taxa))
  ext_rare <- cbind(ext_rare, pad)
}
ext_filt <- ext_rare[,keep_taxa,drop=FALSE]
ext_clr  <- clr_transform(ext_filt)
cat(sprintf("  EXT: %d samples, %d taxa (aligned to DEV taxa)\n\n", nrow(ext_clr), ncol(ext_clr)))

# STEP 5: Stratified 70/30 분할
cat("[5] Stratified train/test split (70/30)...\n")
stratified_split <- function(y, ratio=0.7, seed=42) {
  set.seed(seed)
  pos <- which(y==1); neg <- which(y==0)
  tr  <- c(sample(pos, floor(ratio*length(pos))), sample(neg, floor(ratio*length(neg))))
  list(train=sort(tr), test=setdiff(seq_along(y), sort(tr)))
}
sp      <- stratified_split(dev_y2, TRAIN_RATIO, SEED)
X_train <- dev_clr[sp$train,]; y_train <- dev_y2[sp$train]
X_test  <- dev_clr[sp$test, ]; y_test  <- dev_y2[sp$test]
cat(sprintf("  Train: %d (%d CRC, %d Control)\n", length(y_train), sum(y_train==1), sum(y_train==0)))
cat(sprintf("  Test:  %d (%d CRC, %d Control)\n\n", length(y_test), sum(y_test==1), sum(y_test==0)))

# STEP 6: 모델 학습
cat("[6] Training models...\n")
df_tr  <- data.frame(y=y_train, X_train)
lr_fit <- glm(y ~ ., data=df_tr, family=binomial())
lr_p_train <- as.numeric(predict(lr_fit, data.frame(X_train), type="response"))
lr_p_test  <- as.numeric(predict(lr_fit, data.frame(X_test),  type="response"))
lr_p_ext   <- as.numeric(predict(lr_fit, data.frame(ext_clr), type="response"))
cat("  [6.1] LR done.\n")

set.seed(SEED)
mtry0  <- max(1L, floor(sqrt(ncol(X_train))))
rf_fit <- randomForest(x=X_train, y=as.factor(y_train), ntree=500, mtry=mtry0, importance=TRUE)
rf_p_train <- as.numeric(predict(rf_fit, X_train, type="prob")[,"1"])
rf_p_test  <- as.numeric(predict(rf_fit, X_test,  type="prob")[,"1"])
rf_p_ext   <- as.numeric(predict(rf_fit, ext_clr, type="prob")[,"1"])
cat("  [6.2] RF done.\n")

dtrain <- xgb.DMatrix(as.matrix(X_train), label=y_train)
dtest  <- xgb.DMatrix(as.matrix(X_test),  label=y_test)
dext   <- xgb.DMatrix(as.matrix(ext_clr), label=ext_y2)
set.seed(SEED)
xgb_fit <- xgb.train(
  params=list(booster="gbtree", objective="binary:logistic", eval_metric="auc",
              max_depth=6, eta=0.1, subsample=0.8, colsample_bytree=0.8, min_child_weight=1),
  data=dtrain, nrounds=300, verbose=0
)
xgb_p_train <- as.numeric(predict(xgb_fit, dtrain))
xgb_p_test  <- as.numeric(predict(xgb_fit, dtest))
xgb_p_ext   <- as.numeric(predict(xgb_fit, dext))
cat("  [6.3] XGB done.\n\n")

# STEP 7: Youden J 임계값 선택
youden_threshold <- function(y_true, p_hat) {
  roc_obj <- pROC::roc(y_true, p_hat, quiet=TRUE, direction="<")
  as.numeric(pROC::coords(roc_obj,"best",best.method="youden",ret="threshold",transpose=FALSE))
}
lr_thr  <- youden_threshold(y_train, lr_p_train)
rf_thr  <- youden_threshold(y_train, rf_p_train)
xgb_thr <- youden_threshold(y_train, xgb_p_train)
cat(sprintf("  LR thr=%.4f | RF thr=%.4f | XGB thr=%.4f\n\n", lr_thr, rf_thr, xgb_thr))

# STEP 8: 평가 함수
calibration_ece <- function(y_true, p_hat, n_bins=10) {
  bins <- cut(p_hat, breaks=seq(0,1,length.out=n_bins+1), include.lowest=TRUE, right=TRUE)
  df   <- data.frame(y=y_true, p=p_hat, bin=bins)
  tab  <- df %>% group_by(bin) %>%
    summarise(n=n(), p_mean=mean(p), y_rate=mean(y), .groups="drop") %>%
    mutate(weight=n/sum(n), abs_gap=abs(y_rate-p_mean))
  list(table=tab, ece=sum(tab$weight*tab$abs_gap))
}
evaluate_model <- function(name, y_true, p_hat, threshold, boot_n=2000) {
  roc_obj <- pROC::roc(y_true, p_hat, quiet=TRUE, direction="<")
  auc_val <- as.numeric(pROC::auc(roc_obj))
  ci      <- tryCatch(as.numeric(pROC::ci.auc(roc_obj, method="bootstrap", boot.n=boot_n, conf.level=0.95)),
                      error=function(e) c(NA_real_, auc_val, NA_real_))
  fg     <- p_hat[y_true==1]; bg <- p_hat[y_true==0]
  pr_auc <- as.numeric(PRROC::pr.curve(scores.class0=fg, scores.class1=bg, curve=FALSE)$auc.integral)
  y_pred <- ifelse(p_hat >= threshold, 1L, 0L)
  tp <- sum(y_pred==1L & y_true==1L); tn <- sum(y_pred==0L & y_true==0L)
  fp <- sum(y_pred==1L & y_true==0L); fn <- sum(y_pred==0L & y_true==1L)
  acc    <- (tp+tn)/length(y_true)
  recall <- ifelse(tp+fn==0, NA_real_, tp/(tp+fn))
  prec   <- ifelse(tp+fp==0, NA_real_, tp/(tp+fp))
  f1     <- ifelse(is.na(recall)|is.na(prec)|(prec+recall)==0, NA_real_, 2*prec*recall/(prec+recall))
  brier  <- mean((p_hat-y_true)^2)
  cal    <- calibration_ece(y_true, p_hat, N_BINS_CAL)
  data.frame(Model=name, AUC=round(auc_val,4), AUC_CI_L=round(ci[1],4), AUC_CI_U=round(ci[3],4),
             PR_AUC=round(pr_auc,4), Accuracy=round(acc,4), Recall=round(recall,4),
             Precision=round(prec,4), F1=round(f1,4), Brier=round(brier,4), ECE=round(cal$ece,4),
             Threshold=round(threshold,4), stringsAsFactors=FALSE)
}

# STEP 9: 성능 평가 (Bootstrap CI)
cat("[9] Evaluating all models (bootstrap CI, ~2-5 min)...\n")
res_dev <- bind_rows(evaluate_model("LR",y_test,lr_p_test,lr_thr,BOOT_N),
                     evaluate_model("RF",y_test,rf_p_test,rf_thr,BOOT_N),
                     evaluate_model("XGB",y_test,xgb_p_test,xgb_thr,BOOT_N))
res_ext <- bind_rows(evaluate_model("LR",ext_y2,lr_p_ext,lr_thr,BOOT_N),
                     evaluate_model("RF",ext_y2,rf_p_ext,rf_thr,BOOT_N),
                     evaluate_model("XGB",ext_y2,xgb_p_ext,xgb_thr,BOOT_N))
print(res_dev[,c("Model","AUC","AUC_CI_L","AUC_CI_U","Recall","F1","Brier","ECE")])
print(res_ext[,c("Model","AUC","AUC_CI_L","AUC_CI_U","Recall","F1","Brier","ECE")])
fwrite(res_dev, file.path(OUTDIR,"metrics_dev_test.tsv"),  sep="\t")
fwrite(res_ext, file.path(OUTDIR,"metrics_external.tsv"), sep="\t")

# STEP 10: 통계 검정
cat("\n[10] Statistical testing...\n")
roc_lr_d  <- pROC::roc(y_test, lr_p_test,  quiet=TRUE, direction="<")
roc_rf_d  <- pROC::roc(y_test, rf_p_test,  quiet=TRUE, direction="<")
roc_xgb_d <- pROC::roc(y_test, xgb_p_test, quiet=TRUE, direction="<")
roc_lr_e  <- pROC::roc(ext_y2, lr_p_ext,   quiet=TRUE, direction="<")
roc_xgb_e <- pROC::roc(ext_y2, xgb_p_ext,  quiet=TRUE, direction="<")
dl_xlr_dev <- pROC::roc.test(roc_xgb_d, roc_lr_d,  method="delong")
dl_xrf_dev <- pROC::roc.test(roc_xgb_d, roc_rf_d,  method="delong")
dl_xlr_ext <- pROC::roc.test(roc_xgb_e, roc_lr_e,  method="delong")
pred_lr  <- ifelse(lr_p_test  >= lr_thr,  1L, 0L)
pred_rf  <- ifelse(rf_p_test  >= rf_thr,  1L, 0L)
pred_xgb <- ifelse(xgb_p_test >= xgb_thr, 1L, 0L)
tab_xlr  <- table(pred_xgb, pred_lr)
tab_xrf  <- table(pred_xgb, pred_rf)
mcn_xlr  <- mcnemar.test(tab_xlr)
mcn_xrf  <- mcnemar.test(tab_xrf)
sink(file.path(OUTDIR,"statistical_testing_summary.txt"))
cat("==== DeLong Tests ====\n"); print(dl_xlr_dev); print(dl_xrf_dev); print(dl_xlr_ext)
cat("==== McNemar Tests ====\n"); print(tab_xlr); print(mcn_xlr); print(tab_xrf); print(mcn_xrf)
sink()

# STEP 11: Calibration
ece_summary <- data.frame(Cohort=c(rep("DEV_Test",3),rep("External",3)),
                           Model=c("LR","RF","XGB","LR","RF","XGB"),
                           ECE=c(res_dev$ECE, res_ext$ECE))
fwrite(ece_summary, file.path(OUTDIR,"ece_summary.tsv"), sep="\t")

# STEP 12: 시각화
model_colors <- c(LR="#E07B54", RF="#5B8DB8", XGB="#3A7D44")
make_roc_df <- function(y_true, p_hat, model_name) {
  roc_obj <- pROC::roc(y_true, p_hat, quiet=TRUE, direction="<")
  auc_val <- as.numeric(pROC::auc(roc_obj))
  coords  <- pROC::coords(roc_obj,"all",ret=c("specificity","sensitivity"),transpose=FALSE)
  data.frame(FPR=1-coords$specificity, TPR=coords$sensitivity,
             Model=sprintf("%s (AUC=%.3f)", model_name, auc_val))
}
roc_dev_df <- bind_rows(make_roc_df(y_test,lr_p_test,"LR"), make_roc_df(y_test,rf_p_test,"RF"), make_roc_df(y_test,xgb_p_test,"XGB"))
roc_ext_df <- bind_rows(make_roc_df(ext_y2,lr_p_ext,"LR"), make_roc_df(ext_y2,rf_p_ext,"RF"), make_roc_df(ext_y2,xgb_p_ext,"XGB"))
plot_roc <- function(df, title) {
  ggplot(df, aes(FPR,TPR,color=Model)) +
    geom_abline(linetype="dashed",color="gray50",linewidth=0.7) + geom_line(linewidth=1.3) +
    coord_equal(xlim=c(0,1),ylim=c(0,1)) +
    labs(title=title, x="False Positive Rate", y="True Positive Rate", color=NULL) +
    theme_minimal(base_size=12) + theme(legend.position="inside", legend.position.inside=c(0.65,0.18),
                                        plot.title=element_text(face="bold"))
}
fig1 <- plot_roc(roc_dev_df,"Development Cohort (Test Set)") + plot_roc(roc_ext_df,"External Validation Cohort")
ggsave(file.path(OUTDIR,"Figure1_ROC_curves.png"), fig1, width=13, height=5.5, dpi=300)

make_cal_df <- function(y_true, p_hat, model_name, cohort) {
  tab <- calibration_ece(y_true, p_hat, N_BINS_CAL)$table
  tab$Model <- model_name; tab$Cohort <- cohort; tab
}
cal_df <- bind_rows(make_cal_df(y_test,lr_p_test,"LR","Development"),make_cal_df(y_test,rf_p_test,"RF","Development"),
                    make_cal_df(y_test,xgb_p_test,"XGB","Development"),make_cal_df(ext_y2,lr_p_ext,"LR","External"),
                    make_cal_df(ext_y2,rf_p_ext,"RF","External"),make_cal_df(ext_y2,xgb_p_ext,"XGB","External"))
fig2 <- ggplot(cal_df, aes(p_mean,y_rate,color=Model,group=Model)) +
  geom_abline(linetype="dashed",color="gray40",linewidth=0.7) + geom_line(linewidth=1.1) + geom_point(size=2.5) +
  scale_color_manual(values=model_colors) + coord_equal(xlim=c(0,1),ylim=c(0,1)) + facet_wrap(~Cohort) +
  labs(title="Calibration Reliability Curves", x="Mean Predicted Probability", y="Observed Event Rate", color=NULL) +
  theme_minimal(base_size=12)
ggsave(file.path(OUTDIR,"Figure2_Calibration.png"), fig2, width=11, height=5, dpi=300)

metrics_long <- bind_rows(res_dev%>%mutate(Cohort="Development"), res_ext%>%mutate(Cohort="External")) %>%
  select(Model,Cohort,AUC,PR_AUC,Recall,F1,Accuracy) %>%
  pivot_longer(cols=c(AUC,PR_AUC,Recall,F1,Accuracy), names_to="Metric", values_to="Value") %>%
  mutate(Model=factor(Model,levels=c("LR","RF","XGB")), Cohort=factor(Cohort,levels=c("Development","External")),
         Metric=factor(Metric,levels=c("AUC","PR_AUC","Recall","F1","Accuracy"),labels=c("AUC","PR-AUC","Recall","F1","Accuracy")))
fig3 <- ggplot(metrics_long, aes(Metric,Value,fill=Model)) +
  geom_bar(stat="identity",position=position_dodge(0.75),width=0.7,color="white",linewidth=0.3) +
  geom_text(aes(label=sprintf("%.3f",Value)),position=position_dodge(0.75),vjust=-0.4,size=2.5,fontface="bold") +
  scale_fill_manual(values=model_colors) + scale_y_continuous(limits=c(0,1.12),labels=scales::percent) +
  facet_wrap(~Cohort) + labs(title="Model Performance Comparison",x=NULL,y="Score",fill=NULL) + theme_minimal(base_size=12)
ggsave(file.path(OUTDIR,"Figure3_Performance.png"), fig3, width=13, height=5, dpi=300)

rf_imp <- importance(rf_fit,type=1); top_rf <- sort(rf_imp[,1],decreasing=TRUE)[1:20]
df_rf_imp <- data.frame(Taxa=factor(names(top_rf),levels=rev(names(top_rf))), Importance=as.numeric(top_rf), Model="Random Forest")
xgb_imp_data <- xgb.importance(model=xgb_fit)
df_xgb_imp <- xgb_imp_data[1:min(20,nrow(xgb_imp_data)),] %>% mutate(Model="XGBoost") %>%
  rename(Taxa=Feature, Importance=Gain) %>% mutate(Taxa=factor(Taxa,levels=rev(Taxa)))
fig4 <- (ggplot(df_rf_imp,aes(Importance,Taxa))+geom_bar(stat="identity",fill=model_colors["RF"],alpha=0.85)+
           labs(title="Random Forest\n(Mean Decrease Accuracy)",x="Importance",y=NULL)+theme_minimal(base_size=10)) +
         (ggplot(df_xgb_imp,aes(Importance,Taxa))+geom_bar(stat="identity",fill=model_colors["XGB"],alpha=0.85)+
           labs(title="XGBoost\n(Gain)",x="Importance",y=NULL)+theme_minimal(base_size=10)) +
  plot_annotation(title="Top 20 Feature Importance")
ggsave(file.path(OUTDIR,"Figure4_Feature_importance.png"), fig4, width=14, height=6, dpi=300)

auc_gen <- bind_rows(res_dev%>%mutate(Cohort="Development"), res_ext%>%mutate(Cohort="External")) %>%
  select(Model,Cohort,AUC,AUC_CI_L,AUC_CI_U) %>% mutate(Cohort=factor(Cohort,levels=c("Development","External")))
fig5 <- ggplot(auc_gen,aes(Cohort,AUC,color=Model,group=Model)) +
  geom_line(linewidth=0.8,linetype="dashed",alpha=0.6) + geom_point(size=4) +
  geom_errorbar(aes(ymin=AUC_CI_L,ymax=AUC_CI_U),width=0.12,linewidth=0.9) +
  scale_color_manual(values=model_colors) + scale_y_continuous(limits=c(0.60,1.02)) +
  labs(title="Model Generalizability",subtitle="Error bars: 95% bootstrap CI",x=NULL,y="AUC",color=NULL) +
  theme_minimal(base_size=12)
ggsave(file.path(OUTDIR,"Figure5_Generalizability.png"), fig5, width=7, height=5, dpi=300)

# STEP 13: JSON 저장
final_json <- list(cohort_info=cohort_info,
  split_info=list(train_n=length(y_train),train_crc=sum(y_train==1),train_ctrl=sum(y_train==0),
                  test_n=length(y_test),test_crc=sum(y_test==1),test_ctrl=sum(y_test==0),n_taxa=ncol(dev_clr)),
  thresholds=list(LR=lr_thr,RF=rf_thr,XGB=xgb_thr), dev_results=res_dev, ext_results=res_ext,
  statistical_tests=list(
    delong_dev_xgb_vs_lr=list(delta_auc=round(diff(dl_xlr_dev$estimate),4), p_value=round(dl_xlr_dev$p.value,4)),
    delong_dev_xgb_vs_rf=list(delta_auc=round(diff(dl_xrf_dev$estimate),4), p_value=round(dl_xrf_dev$p.value,4)),
    delong_ext_xgb_vs_lr=list(delta_auc=round(diff(dl_xlr_ext$estimate),4), p_value=round(dl_xlr_ext$p.value,4)),
    mcnemar_dev_xgb_vs_lr=list(p_value=round(mcn_xlr$p.value,4)),
    mcnemar_dev_xgb_vs_rf=list(p_value=round(mcn_xrf$p.value,4))
  ))
write(jsonlite::toJSON(final_json,pretty=TRUE,auto_unbox=TRUE), file.path(OUTDIR,"final_results.json"))

cat("\n========================================\n")
cat("  FINAL SUMMARY\n")
cat("========================================\n")
print(res_dev[,c("Model","AUC","AUC_CI_L","AUC_CI_U","Recall","F1","ECE")])
print(res_ext[,c("Model","AUC","AUC_CI_L","AUC_CI_U","Recall","F1","ECE")])
cat(sprintf("\n모든 결과가 저장되었습니다: %s\n", OUTDIR))

# EXTRA: add_figures.py 호환 파일 저장
FIG_DIR <- file.path(OUTDIR,"figures")
dir.create(FIG_DIR, showWarnings=FALSE, recursive=TRUE)
roc_rf_e <- pROC::roc(ext_y2, rf_p_ext, quiet=TRUE, direction="<")
saveRDS(list(dev=list(LR=roc_lr_d,RF=roc_rf_d,XGB=roc_xgb_d), ext=list(LR=roc_lr_e,RF=roc_rf_e,XGB=roc_xgb_e)),
        file.path(OUTDIR,"roc_curves.rds"))
saveRDS(xgb.importance(model=xgb_fit), file.path(OUTDIR,"xgb_importance.rds"))
fig_rename <- list("Figure1_ROC_curves.png"="fig1_roc_dev.png", "Figure2_Calibration.png"="fig_cal.png",
                   "Figure3_Performance.png"="fig3_perf_barplot.png",
                   "Figure4_Feature_importance.png"="fig4_feature_importance.png",
                   "Figure5_Generalizability.png"="fig5_generalizability.png")
for (src_name in names(fig_rename)) {
  src <- file.path(OUTDIR,src_name); dst <- file.path(FIG_DIR,fig_rename[[src_name]])
  if (file.exists(src)) file.copy(src, dst, overwrite=TRUE)
}
cat("\n[완료] 그래프 삽입 준비 완료!\n")
cat("다음 단계: source('R_scripts/generate_figures.R') → python add_figures.py\n")
