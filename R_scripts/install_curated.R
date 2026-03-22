# =============================================================================
# install_curated.R
# curatedMetagenomicData 설치 및 데이터 로드 확인
# bluster 설치 완료 후 이 스크립트를 실행하세요.
# =============================================================================

cat("=== curatedMetagenomicData 설치 시작 ===\n")
cat("R 버전:", R.version$version.string, "\n\n")

is_installed <- function(pkg) requireNamespace(pkg, quietly = TRUE)

# ─── STEP 1: curatedMetagenomicData 설치 ─────────────────────────────────────
if (!is_installed("curatedMetagenomicData")) {
  cat("[STEP 1] curatedMetagenomicData 설치 중...\n")
  BiocManager::install("curatedMetagenomicData", update = FALSE, ask = FALSE)
} else {
  cat("[STEP 1] curatedMetagenomicData 이미 설치됨\n")
}

# ─── STEP 2: 기타 필요 패키지 확인 ──────────────────────────────────────────
cat("\n[STEP 2] 분석 필요 패키지 확인 중...\n")
required <- c("pROC", "PRROC", "randomForest", "xgboost",
              "ggplot2", "vegan", "dplyr", "data.table",
              "tidyr", "jsonlite")
missing <- required[!sapply(required, is_installed)]
if (length(missing) > 0) {
  cat("  설치 필요:", paste(missing, collapse=", "), "\n")
  install.packages(missing,
                   repos = "https://cloud.r-project.org",
                   type  = "binary")
} else {
  cat("  모든 패키지 설치됨 ✓\n")
}

# ─── STEP 3: 데이터 로드 테스트 ──────────────────────────────────────────────
cat("\n[STEP 3] 데이터 로드 테스트...\n")
library(curatedMetagenomicData)

tryCatch({
  # 개발 코호트 - ZellerG_2014
  cat("  ZellerG_2014 로드 중...\n")
  zeller <- curatedMetagenomicData("ZellerG_2014.relative_abundance",
                                   dryrun = FALSE)
  cat("  → ZellerG_2014 샘플 수:", ncol(zeller[[1]]), "\n")

  # 외부 검증 코호트 - YachidaS_2019
  cat("  YachidaS_2019 로드 중...\n")
  yachida <- curatedMetagenomicData("YachidaS_2019.relative_abundance",
                                    dryrun = FALSE)
  cat("  → YachidaS_2019 샘플 수:", ncol(yachida[[1]]), "\n")

  cat("\n[SUCCESS] 두 코호트 로드 완료! 분석 실행 준비가 되었습니다.\n")
  cat("다음 명령을 실행하세요:\n")
  cat("  source('C:/Users/minta/ICBBS/R_scripts/run_analysis.R')\n")

}, error = function(e) {
  cat("[ERROR] 데이터 로드 실패:", conditionMessage(e), "\n")
  cat("인터넷 연결을 확인하고 다시 시도하세요.\n")
})
