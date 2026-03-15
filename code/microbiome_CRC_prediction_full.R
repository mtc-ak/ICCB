# ============================================================
# Microbiome-Based CRC Prediction: Complete R Pipeline
# ============================================================
# 논문: 마이크로바이옴 기반 대장암 예측을 위한
#       머신러닝 모델의 통계적 비교 및 외부 검증 연구
# ------------------------------------------------------------
# 개발 코호트: ENA ERP005534 (ZellerG_2014, n=290)
# 외부 검증:   PRJDB4176    (YachidaS_2019, n=212)
# 모델:        LR, RF, XGBoost
# 지표:        AUC, PR-AUC, Accuracy, Recall, F1, Brier, ECE
# 검정:        DeLong (pROC), McNemar
# Bootstrap:   2,000 iterations, 95% CI
# Seed:        42
# ============================================================

# ─────────────────────────────────────────────
# [STEP 0]  패키지 설치 및 로드
# ─────────────────────────────────────────────
# 최초 실행 시 아래 주석을 해제하여 설치하세요.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("curatedMetagenomicData",
#                        "TreeSummarizedExperiment"))
# install.packages(c("data.table","dplyr","vegan","pROC",
#                    "PRROC","randomForest","xgboost",
#                    "ggplot2","patchwork","scales"))

suppressPackageStartupMessages({
  library(curatedMetagenomicData)   # 공개 마이크로바이옴 데이터
  library(TreeSummarizedExperiment) # SummarizedExperiment 확장
  library(data.table)
  library(dplyr)
  library(vegan)          # rrarefy (rarefaction)
  library(pROC)           # AUC, DeLong 검정
  library(PRROC)          # PR-AUC
  library(randomForest)   # Random Forest
  library(xgboost)        # XGBoost
  library(ggplot2)
  library(patchwork)      # 그래프 결합
  library(scales)
})

set.seed(42)

# ─────────────────────────────────────────────
# [STEP 1]  설정값 (User Config)
# ─────────────────────────────────────────────
RARE_DEPTH  <- 10000    # Rarefaction depth
PREV_CUT    <- 0.05     # Prevalence filter cutoff (5%)
N_BINS_CAL  <- 10       # Calibration bins
BOOT_N      <- 2000     # Bootstrap iterations
TRAIN_RATIO <- 0.70     # Train/Test split ratio
POS_LABEL   <- "CRC"    # 양성 클래스 레이블

OUTDIR <- "results_crc_microbiome"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

cat("========================================\n")
cat("  Microbiome CRC Prediction Pipeline\n")
cat("========================================\n\n")

# ─────────────────────────────────────────────
# [STEP 2]  데이터 로드 (curatedMetagenomicData)
# ─────────────────────────────────────────────
cat("[2] Loading public datasets...\n")

## 2-1. 개발 코호트: ZellerG_2014 (ERP005534)
dev_se <- curatedMetagenomicData(
  "ZellerG_2014.relative_abundance",
  dryrun = FALSE,
  counts  = TRUE   # 원시 count로 변환
)[[1]]

dev_counts <- t(assay(dev_se))   # rows=samples, cols=taxa
dev_meta   <- as.data.frame(colData(dev_se))
dev_y      <- ifelse(dev_meta$study_condition == POS_LABEL, 1L, 0L)
names(dev_y) <- rownames(dev_meta)

cat(sprintf("  DEV cohort: %d samples (%d CRC, %d Control), %d taxa\n",
            nrow(dev_counts),
            sum(dev_y == 1), sum(dev_y == 0),
            ncol(dev_counts)))

## 2-2. 외부 검증 코호트: YachidaS_2019 (PRJDB4176)
##      CRC + control 샘플만 선택 (adenoma 제외)
ext_se   <- curatedMetagenomicData(
  "YachidaS_2019.relative_abundance",
  dryrun = FALSE,
  counts  = TRUE
)[[1]]

ext_keep <- colData(ext_se)$study_condition %in% c("CRC", "control")
ext_se   <- ext_se[, ext_keep]

ext_counts <- t(assay(ext_se))
ext_meta   <- as.data.frame(colData(ext_se))
ext_y      <- ifelse(ext_meta$study_condition == POS_LABEL, 1L, 0L)
names(ext_y) <- rownames(ext_meta)

cat(sprintf("  EXT cohort: %d samples (%d CRC, %d Control), %d taxa\n\n",
            nrow(ext_counts),
            sum(ext_y == 1), sum(ext_y == 0),
            ncol(ext_counts)))

# ─────────────────────────────────────────────
# [STEP 3]  전처리 함수 정의
# ─────────────────────────────────────────────
cat("[3] Defining preprocessing functions...\n")

## Rarefaction
safe_rarefy <- function(counts_mat, depth) {
  rs   <- rowSums(counts_mat)
  keep <- rs >= depth
  if (!all(keep))
    message(sprintf("  [WARN] %d samples dropped (reads < %d)",
                    sum(!keep), depth))
  vegan::rrarefy(counts_mat[keep, , drop = FALSE], sample = depth)
}

## Prevalence filter  (DEV 기준으로 정의, EXT에 동일 적용)
prevalence_filter <- function(counts_mat, cut) {
  prev <- colSums(counts_mat > 0) / nrow(counts_mat)
  keep <- names(prev)[prev >= cut]
  list(filtered = counts_mat[, keep, drop = FALSE],
       kept_taxa = keep)
}

## CLR 변환 (pseudocount = 1)
clr_transform <- function(counts_mat, pseudocount = 1) {
  X     <- counts_mat + pseudocount
  log_X <- log(X)
  sweep(log_X, 1, rowMeans(log_X), "-")
}

# ─────────────────────────────────────────────
# [STEP 4]  전처리 파이프라인 실행
# ─────────────────────────────────────────────
cat("[4] Preprocessing...\n")

## 4-1. DEV: Rarefy → Prevalence filter → CLR
dev_rare   <- safe_rarefy(dev_counts, RARE_DEPTH)
dev_y2     <- dev_y[rownames(dev_rare)]          # 제거된 샘플 동기화

pf         <- prevalence_filter(dev_rare, PREV_CUT)
dev_filt   <- pf$filtered
keep_taxa  <- pf$kept_taxa
dev_clr    <- clr_transform(dev_filt)

cat(sprintf("  DEV: %d samples, %d taxa (after filtering)\n",
            nrow(dev_clr), ncol(dev_clr)))

## 4-2. EXT: 동일 파이프라인 + DEV taxa 집합 강제 정렬
ext_rare   <- safe_rarefy(ext_counts, RARE_DEPTH)
ext_y2     <- ext_y[rownames(ext_rare)]

# 누락 taxa → 0 패딩; 초과 taxa → 제거
miss_taxa  <- setdiff(keep_taxa, colnames(ext_rare))
if (length(miss_taxa) > 0) {
  message(sprintf("  [WARN] EXT missing %d taxa from DEV; padding with 0.",
                  length(miss_taxa)))
  pad <- matrix(0L, nrow = nrow(ext_rare), ncol = length(miss_taxa),
                dimnames = list(rownames(ext_rare), miss_taxa))
  ext_rare <- cbind(ext_rare, pad)
}
ext_filt   <- ext_rare[, keep_taxa, drop = FALSE]
ext_clr    <- clr_transform(ext_filt)

cat(sprintf("  EXT: %d samples, %d taxa (aligned to DEV)\n\n",
            nrow(ext_clr), ncol(ext_clr)))

# ─────────────────────────────────────────────
# [STEP 5]  Stratified 70 / 30 분할
# ─────────────────────────────────────────────
cat("[5] Stratified train/test split (70/30)...\n")

stratified_split <- function(y, ratio = 0.7, seed = 42) {
  set.seed(seed)
  pos <- which(y == 1); neg <- which(y == 0)
  tr  <- c(sample(pos, floor(ratio * length(pos))),
           sample(neg, floor(ratio * length(neg))))
  list(train = sort(tr),
       test  = setdiff(seq_along(y), sort(tr)))
}

sp      <- stratified_split(dev_y2, TRAIN_RATIO, seed = 42)
X_train <- dev_clr[sp$train, ]; y_train <- dev_y2[sp$train]
X_test  <- dev_clr[sp$test,  ]; y_test  <- dev_y2[sp$test]

cat(sprintf("  Train: %d (%d CRC) | Test: %d (%d CRC)\n\n",
            length(sp$train), sum(y_train),
            length(sp$test),  sum(y_test)))

# ─────────────────────────────────────────────
# [STEP 6]  모델 학습
# ─────────────────────────────────────────────
cat("[6] Training models...\n")

## 6-1. Logistic Regression (L2)
cat("  [6.1] Logistic Regression\n")
df_train  <- data.frame(y = y_train, X_train)
lr_fit    <- glm(y ~ ., data = df_train, family = binomial())

lr_p_train <- as.numeric(predict(lr_fit, data.frame(X_train), type = "response"))
lr_p_test  <- as.numeric(predict(lr_fit, data.frame(X_test),  type = "response"))
lr_p_ext   <- as.numeric(predict(lr_fit, data.frame(ext_clr), type = "response"))

## 6-2. Random Forest (500 trees)
cat("  [6.2] Random Forest (ntree=500)\n")
set.seed(42)
mtry0    <- max(1L, floor(sqrt(ncol(X_train))))
rf_fit   <- randomForest(
  x = X_train, y = as.factor(y_train),
  ntree = 500, mtry = mtry0,
  importance = TRUE
)

rf_p_train <- as.numeric(predict(rf_fit, X_train,  type = "prob")[, "1"])
rf_p_test  <- as.numeric(predict(rf_fit, X_test,   type = "prob")[, "1"])
rf_p_ext   <- as.numeric(predict(rf_fit, ext_clr,  type = "prob")[, "1"])

## 6-3. XGBoost (max_depth=6, lr=0.1, n=300)
cat("  [6.3] XGBoost (max_depth=6, eta=0.1, nrounds=300)\n")
dtrain <- xgb.DMatrix(as.matrix(X_train), label = y_train)
dtest  <- xgb.DMatrix(as.matrix(X_test),  label = y_test)
dext   <- xgb.DMatrix(as.matrix(ext_clr), label = ext_y2)

set.seed(42)
xgb_fit <- xgb.train(
  params = list(
    booster          = "gbtree",
    objective        = "binary:logistic",
    eval_metric      = "auc",
    max_depth        = 6,
    eta              = 0.1,
    subsample        = 0.8,
    colsample_bytree = 0.8,
    min_child_weight = 1
  ),
  data    = dtrain,
  nrounds = 300,
  verbose = 0
)

xgb_p_train <- as.numeric(predict(xgb_fit, dtrain))
xgb_p_test  <- as.numeric(predict(xgb_fit, dtest))
xgb_p_ext   <- as.numeric(predict(xgb_fit, dext))
cat("  Models trained.\n\n")

# ─────────────────────────────────────────────
# [STEP 7]  임계값 선택 (Youden J — 훈련셋 기준)
# ─────────────────────────────────────────────
cat("[7] Selecting thresholds (Youden J on train set)...\n")

youden_threshold <- function(y_true, p_hat) {
  roc_obj <- pROC::roc(y_true, p_hat, quiet = TRUE, direction = "<")
  as.numeric(pROC::coords(roc_obj, x = "best",
                           best.method = "youden",
                           ret = "threshold",
                           transpose = FALSE))
}

lr_thr  <- youden_threshold(y_train, lr_p_train)
rf_thr  <- youden_threshold(y_train, rf_p_train)
xgb_thr <- youden_threshold(y_train, xgb_p_train)

cat(sprintf("  LR thr=%.4f | RF thr=%.4f | XGB thr=%.4f\n\n",
            lr_thr, rf_thr, xgb_thr))

# ─────────────────────────────────────────────
# [STEP 8]  평가 함수 정의
# ─────────────────────────────────────────────

## ECE (Expected Calibration Error)
calibration_ece <- function(y_true, p_hat, n_bins = 10) {
  bins <- cut(p_hat,
              breaks = seq(0, 1, length.out = n_bins + 1),
              include.lowest = TRUE, right = TRUE)
  df  <- data.frame(y = y_true, p = p_hat, bin = bins)
  tab <- df %>%
    group_by(bin) %>%
    summarise(n     = n(),
              p_mean = mean(p),
              y_rate = mean(y),
              .groups = "drop") %>%
    mutate(weight  = n / sum(n),
           abs_gap = abs(y_rate - p_mean))
  list(table = tab, ece = sum(tab$weight * tab$abs_gap))
}

## 단일 모델 평가 (AUC CI 포함)
evaluate_model <- function(name, y_true, p_hat, threshold,
                            boot_n = 2000) {
  ## ROC & AUC
  roc_obj <- pROC::roc(y_true, p_hat, quiet = TRUE, direction = "<")
  auc_val <- as.numeric(pROC::auc(roc_obj))

  ## Bootstrap CI
  if (boot_n > 0) {
    ci <- tryCatch(
      as.numeric(pROC::ci.auc(roc_obj,
                               method  = "bootstrap",
                               boot.n  = boot_n,
                               conf.level = 0.95)),
      error = function(e) c(NA_real_, auc_val, NA_real_)
    )
  } else {
    ci <- c(NA_real_, auc_val, NA_real_)
  }

  ## PR-AUC
  fg     <- p_hat[y_true == 1]
  bg     <- p_hat[y_true == 0]
  pr_obj <- PRROC::pr.curve(scores.class0 = fg,
                              scores.class1 = bg,
                              curve = FALSE)
  pr_auc <- as.numeric(pr_obj$auc.integral)

  ## Classification metrics
  y_pred    <- ifelse(p_hat >= threshold, 1L, 0L)
  tp <- sum(y_pred == 1L & y_true == 1L)
  tn <- sum(y_pred == 0L & y_true == 0L)
  fp <- sum(y_pred == 1L & y_true == 0L)
  fn <- sum(y_pred == 0L & y_true == 1L)

  acc    <- (tp + tn) / length(y_true)
  recall <- ifelse(tp + fn == 0, NA_real_, tp / (tp + fn))
  prec   <- ifelse(tp + fp == 0, NA_real_, tp / (tp + fp))
  f1     <- ifelse(is.na(recall) | is.na(prec) | (prec + recall) == 0,
                   NA_real_,
                   2 * prec * recall / (prec + recall))
  brier  <- mean((p_hat - y_true)^2)

  ## ECE
  cal_res <- calibration_ece(y_true, p_hat, N_BINS_CAL)

  data.frame(
    Model      = name,
    AUC        = round(auc_val, 4),
    AUC_CI_L   = round(ci[1],   4),
    AUC_CI_U   = round(ci[3],   4),
    PR_AUC     = round(pr_auc,  4),
    Accuracy   = round(acc,     4),
    Recall     = round(recall,  4),
    Precision  = round(prec,    4),
    F1         = round(f1,      4),
    Brier      = round(brier,   4),
    ECE        = round(cal_res$ece, 4),
    Threshold  = round(threshold, 4),
    stringsAsFactors = FALSE
  )
}

# ─────────────────────────────────────────────
# [STEP 9]  성능 평가 실행
# ─────────────────────────────────────────────
cat("[9] Evaluating models (bootstrap CI takes a few minutes)...\n")

res_dev <- bind_rows(
  evaluate_model("LR",  y_test, lr_p_test,  lr_thr,  BOOT_N),
  evaluate_model("RF",  y_test, rf_p_test,  rf_thr,  BOOT_N),
  evaluate_model("XGB", y_test, xgb_p_test, xgb_thr, BOOT_N)
)
res_ext <- bind_rows(
  evaluate_model("LR",  ext_y2, lr_p_ext,  lr_thr,  BOOT_N),
  evaluate_model("RF",  ext_y2, rf_p_ext,  rf_thr,  BOOT_N),
  evaluate_model("XGB", ext_y2, xgb_p_ext, xgb_thr, BOOT_N)
)

cat("\n=== Development Cohort Test Results ===\n")
print(res_dev[, c("Model","AUC","AUC_CI_L","AUC_CI_U",
                  "PR_AUC","Recall","F1","Brier","ECE")])

cat("\n=== External Validation Results ===\n")
print(res_ext[, c("Model","AUC","AUC_CI_L","AUC_CI_U",
                  "PR_AUC","Recall","F1","Brier","ECE")])

## 결과 저장
fwrite(res_dev, file.path(OUTDIR, "metrics_dev_test.tsv"), sep = "\t")
fwrite(res_ext, file.path(OUTDIR, "metrics_external.tsv"), sep = "\t")
cat("\n  Metrics saved.\n")

# ─────────────────────────────────────────────
# [STEP 10]  통계 검정
# ─────────────────────────────────────────────
cat("\n[10] Statistical testing (DeLong + McNemar)...\n")

## ROC 객체
roc_lr_d  <- pROC::roc(y_test, lr_p_test,  quiet = TRUE, direction = "<")
roc_rf_d  <- pROC::roc(y_test, rf_p_test,  quiet = TRUE, direction = "<")
roc_xgb_d <- pROC::roc(y_test, xgb_p_test, quiet = TRUE, direction = "<")
roc_lr_e  <- pROC::roc(ext_y2, lr_p_ext,   quiet = TRUE, direction = "<")
roc_xgb_e <- pROC::roc(ext_y2, xgb_p_ext,  quiet = TRUE, direction = "<")

## DeLong 검정
dl_xlr_dev <- pROC::roc.test(roc_xgb_d, roc_lr_d, method = "delong")
dl_xrf_dev <- pROC::roc.test(roc_xgb_d, roc_rf_d, method = "delong")
dl_xlr_ext <- pROC::roc.test(roc_xgb_e, roc_lr_e, method = "delong")

## McNemar 검정 (DEV test set, 임계값 기반)
pred_lr  <- ifelse(lr_p_test  >= lr_thr,  1L, 0L)
pred_rf  <- ifelse(rf_p_test  >= rf_thr,  1L, 0L)
pred_xgb <- ifelse(xgb_p_test >= xgb_thr, 1L, 0L)

tab_xlr  <- table(pred_xgb, pred_lr)
tab_xrf  <- table(pred_xgb, pred_rf)
mcn_xlr  <- mcnemar.test(tab_xlr)
mcn_xrf  <- mcnemar.test(tab_xrf)

## 결과 출력
cat("\n--- DeLong Tests ---\n")
cat(sprintf("  DEV XGB vs LR : ΔAUC=%+.4f, p=%s\n",
            diff(dl_xlr_dev$estimate),
            format.pval(dl_xlr_dev$p.value, digits = 4)))
cat(sprintf("  DEV XGB vs RF : ΔAUC=%+.4f, p=%s\n",
            diff(dl_xrf_dev$estimate),
            format.pval(dl_xrf_dev$p.value, digits = 4)))
cat(sprintf("  EXT XGB vs LR : ΔAUC=%+.4f, p=%s\n",
            diff(dl_xlr_ext$estimate),
            format.pval(dl_xlr_ext$p.value, digits = 4)))

cat("\n--- McNemar Tests (DEV test set) ---\n")
cat(sprintf("  XGB vs LR : p=%s  (b=%d, c=%d)\n",
            format.pval(mcn_xlr$p.value, digits = 4),
            tab_xlr[2, 1], tab_xlr[1, 2]))
cat(sprintf("  XGB vs RF : p=%s  (b=%d, c=%d)\n",
            format.pval(mcn_xrf$p.value, digits = 4),
            tab_xrf[2, 1], tab_xrf[1, 2]))

## 통계 검정 결과 저장
sink(file.path(OUTDIR, "statistical_testing_summary.txt"))
cat("==== DeLong Tests ====\n")
print(dl_xlr_dev); cat("\n")
print(dl_xrf_dev); cat("\n")
print(dl_xlr_ext); cat("\n")
cat("==== McNemar Tests ====\n")
print(tab_xlr); print(mcn_xlr); cat("\n")
print(tab_xrf); print(mcn_xrf)
sink()
cat("\n  Statistical tests saved.\n")

# ─────────────────────────────────────────────
# [STEP 11]  Calibration 분석
# ─────────────────────────────────────────────
cat("\n[11] Calibration analysis...\n")

## ECE 요약표
ece_summary <- data.frame(
  Cohort = c(rep("DEV_Test", 3), rep("External", 3)),
  Model  = c("LR","RF","XGB","LR","RF","XGB"),
  ECE    = c(res_dev$ECE, res_ext$ECE)
)
fwrite(ece_summary, file.path(OUTDIR, "ece_summary.tsv"), sep = "\t")
cat("  ECE Summary:\n")
print(ece_summary)

# ─────────────────────────────────────────────
# [STEP 12]  시각화
# ─────────────────────────────────────────────
cat("\n[12] Generating plots...\n")

model_colors <- c(LR = "#E07B54", RF = "#5B8DB8", XGB = "#3A7D44")

## ── Figure 1: ROC Curves ──────────────────────
make_roc_df <- function(y_true, p_hat, model_name) {
  roc_obj <- pROC::roc(y_true, p_hat, quiet = TRUE, direction = "<")
  auc_val <- as.numeric(pROC::auc(roc_obj))
  coords  <- pROC::coords(roc_obj, "all",
                           ret = c("specificity","sensitivity"),
                           transpose = FALSE)
  data.frame(
    FPR   = 1 - coords$specificity,
    TPR   = coords$sensitivity,
    Model = sprintf("%s (AUC=%.3f)", model_name, auc_val)
  )
}

roc_dev <- bind_rows(
  make_roc_df(y_test, lr_p_test,  "LR"),
  make_roc_df(y_test, rf_p_test,  "RF"),
  make_roc_df(y_test, xgb_p_test, "XGB")
)
roc_ext <- bind_rows(
  make_roc_df(ext_y2, lr_p_ext,  "LR"),
  make_roc_df(ext_y2, rf_p_ext,  "RF"),
  make_roc_df(ext_y2, xgb_p_ext, "XGB")
)

plot_roc <- function(df, title) {
  ggplot(df, aes(FPR, TPR, color = Model)) +
    geom_abline(linetype = "dashed", color = "gray50", linewidth = 0.7) +
    geom_line(linewidth = 1.3) +
    scale_color_manual(
      values = setNames(
        model_colors[c("LR","RF","XGB")],
        unique(grep("LR|RF|XGB", df$Model, value = TRUE))
      )
    ) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title,
         x = "False Positive Rate (1 – Specificity)",
         y = "True Positive Rate (Sensitivity)",
         color = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.65, 0.18),
          plot.title = element_text(face = "bold"))
}

fig1 <- plot_roc(roc_dev, "Development Cohort (Test Set)") +
        plot_roc(roc_ext, "External Validation Cohort")

ggsave(file.path(OUTDIR, "Figure1_ROC_curves.png"),
       fig1, width = 13, height = 5.5, dpi = 300)
cat("  Figure 1 (ROC) saved.\n")

## ── Figure 2: Calibration Reliability Curves ──
make_cal_df <- function(y_true, p_hat, model_name, cohort) {
  cal <- calibration_ece(y_true, p_hat, N_BINS_CAL)$table
  cal$Model  <- model_name
  cal$Cohort <- cohort
  cal
}

cal_df <- bind_rows(
  make_cal_df(y_test, lr_p_test,  "LR",  "Development"),
  make_cal_df(y_test, rf_p_test,  "RF",  "Development"),
  make_cal_df(y_test, xgb_p_test, "XGB", "Development"),
  make_cal_df(ext_y2, lr_p_ext,  "LR",  "External"),
  make_cal_df(ext_y2, rf_p_ext,  "RF",  "External"),
  make_cal_df(ext_y2, xgb_p_ext, "XGB", "External")
)

fig2 <- ggplot(cal_df, aes(p_mean, y_rate, color = Model, group = Model)) +
  geom_abline(linetype = "dashed", color = "gray40", linewidth = 0.7) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = model_colors) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  facet_wrap(~ Cohort) +
  labs(title = "Calibration Reliability Curves",
       x = "Mean Predicted Probability",
       y = "Observed Event Rate",
       color = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

ggsave(file.path(OUTDIR, "Figure2_Calibration.png"),
       fig2, width = 11, height = 5, dpi = 300)
cat("  Figure 2 (Calibration) saved.\n")

## ── Figure 3: Multi-metric Bar Chart ──────────
metrics_long <- bind_rows(
  res_dev %>% mutate(Cohort = "Development"),
  res_ext %>% mutate(Cohort = "External")
) %>%
  select(Model, Cohort, AUC, PR_AUC, Recall, F1, Accuracy) %>%
  tidyr::pivot_longer(cols = c(AUC, PR_AUC, Recall, F1, Accuracy),
                      names_to = "Metric", values_to = "Value") %>%
  mutate(Model  = factor(Model,  levels = c("LR","RF","XGB")),
         Cohort = factor(Cohort, levels = c("Development","External")),
         Metric = factor(Metric, levels = c("AUC","PR_AUC","Recall","F1","Accuracy"),
                         labels = c("AUC","PR-AUC","Recall","F1","Accuracy")))

fig3 <- ggplot(metrics_long, aes(Metric, Value, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(0.75),
           width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", Value)),
            position = position_dodge(0.75),
            vjust = -0.4, size = 2.5, fontface = "bold") +
  scale_fill_manual(values = model_colors) +
  scale_y_continuous(limits = c(0, 1.1), labels = scales::percent) +
  facet_wrap(~ Cohort) +
  labs(title = "Model Performance Comparison",
       x = NULL, y = "Score", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        strip.text    = element_text(face = "bold"),
        axis.text.x   = element_text(angle = 15, hjust = 1))

ggsave(file.path(OUTDIR, "Figure3_Performance.png"),
       fig3, width = 13, height = 5, dpi = 300)
cat("  Figure 3 (Performance bar chart) saved.\n")

## ── Figure 4: Feature Importance (RF + XGB) ───
## RF importance (mean decrease Gini)
rf_imp    <- importance(rf_fit, type = 1)   # MeanDecreaseAccuracy
top_rf    <- sort(rf_imp[, 1], decreasing = TRUE)[1:20]
df_rf_imp <- data.frame(
  Taxa       = factor(names(top_rf), levels = rev(names(top_rf))),
  Importance = as.numeric(top_rf),
  Model      = "Random Forest"
)

## XGB importance (gain)
xgb_imp    <- xgb.importance(model = xgb_fit)
df_xgb_imp <- xgb_imp[1:min(20, nrow(xgb_imp)), ] %>%
  mutate(Model = "XGBoost") %>%
  rename(Taxa = Feature, Importance = Gain) %>%
  mutate(Taxa = factor(Taxa, levels = rev(Taxa)))

fig4_rf  <- ggplot(df_rf_imp, aes(Importance, Taxa)) +
  geom_bar(stat = "identity", fill = model_colors["RF"], alpha = 0.85) +
  labs(title = "Random Forest\n(Mean Decrease Accuracy)",
       x = "Importance", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10))

fig4_xgb <- ggplot(df_xgb_imp, aes(Importance, Taxa)) +
  geom_bar(stat = "identity", fill = model_colors["XGB"], alpha = 0.85) +
  labs(title = "XGBoost\n(Gain)",
       x = "Importance", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 10))

fig4 <- fig4_rf + fig4_xgb +
  plot_annotation(title = "Top 20 Feature Importance",
                  theme = theme(plot.title = element_text(face = "bold")))

ggsave(file.path(OUTDIR, "Figure4_Feature_importance.png"),
       fig4, width = 14, height = 6, dpi = 300)
cat("  Figure 4 (Feature importance) saved.\n")

## ── Figure 5: Generalizability (AUC drop DEV→EXT) ──
auc_gen <- bind_rows(
  res_dev %>% mutate(Cohort = "Development"),
  res_ext %>% mutate(Cohort = "External")
) %>%
  select(Model, Cohort, AUC, AUC_CI_L, AUC_CI_U) %>%
  mutate(Cohort = factor(Cohort, levels = c("Development","External")))

fig5 <- ggplot(auc_gen,
               aes(Cohort, AUC, color = Model, group = Model)) +
  geom_line(linewidth = 0.8, linetype = "dashed", alpha = 0.6) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = AUC_CI_L, ymax = AUC_CI_U),
                width = 0.12, linewidth = 0.9) +
  scale_color_manual(values = model_colors) +
  scale_y_continuous(limits = c(0.70, 1.02),
                     labels = scales::number_format(accuracy = 0.01)) +
  labs(title = "Model Generalizability: Development → External Validation",
       subtitle = "Error bars: 95% bootstrap CI",
       x = NULL, y = "AUC", color = NULL) +
  theme_minimal(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(color = "gray50"))

ggsave(file.path(OUTDIR, "Figure5_Generalizability.png"),
       fig5, width = 7, height = 5, dpi = 300)
cat("  Figure 5 (Generalizability) saved.\n")

# ─────────────────────────────────────────────
# [STEP 13]  최종 요약 출력
# ─────────────────────────────────────────────
cat("\n========================================\n")
cat("  FINAL SUMMARY\n")
cat("========================================\n")
cat("\n[Development Cohort Test Set]\n")
print(res_dev[, c("Model","AUC","AUC_CI_L","AUC_CI_U",
                  "Recall","F1","Brier","ECE")])
cat("\n[External Validation Cohort]\n")
print(res_ext[, c("Model","AUC","AUC_CI_L","AUC_CI_U",
                  "Recall","F1","Brier","ECE")])
cat("\n[Statistical Tests]\n")
cat(sprintf("  DeLong  DEV XGB vs LR: p=%.4f\n", dl_xlr_dev$p.value))
cat(sprintf("  DeLong  DEV XGB vs RF: p=%.4f\n", dl_xrf_dev$p.value))
cat(sprintf("  DeLong  EXT XGB vs LR: p=%.4f\n", dl_xlr_ext$p.value))
cat(sprintf("  McNemar DEV XGB vs LR: p=%.4f\n", mcn_xlr$p.value))
cat(sprintf("  McNemar DEV XGB vs RF: p=%.4f\n", mcn_xrf$p.value))
cat(sprintf("\nAll outputs saved in: %s\n", OUTDIR))
cat("========================================\n")

# ─────────────────────────────────────────────
# END OF SCRIPT
# ─────────────────────────────────────────────
