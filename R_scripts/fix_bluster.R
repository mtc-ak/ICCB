# =============================================================================
# fix_bluster.R
# bluster 패키지 설치 문제 해결 스크립트
#
# 문제: curatedMetagenomicData 의존성인 bluster 가
#       Bioconductor 3.22 에서 바이너리 제공 안 됨
#       → 소스 컴파일 필요 (Rtools 4.5 required on Windows)
#
# 실행 순서: 이 스크립트를 RStudio에서 source() 로 실행하세요.
# =============================================================================

cat("=== bluster 설치 수리 스크립트 시작 ===\n")
cat("R 버전:", R.version$version.string, "\n\n")

is_installed <- function(pkg) requireNamespace(pkg, quietly = TRUE)

cat("[1] 현재 설치 상태 확인:\n")
pkgs_check <- c("bluster", "BiocNeighbors", "curatedMetagenomicData")
for (p in pkgs_check) {
  status <- if (is_installed(p)) "✓ 설치됨" else "✗ 없음"
  cat(sprintf("    %-30s %s\n", p, status))
}
cat("\n")

if (is_installed("bluster")) {
  cat("[OK] bluster 이미 설치되어 있습니다. 종료합니다.\n")
  quit(save = "no")
}

if (!is_installed("BiocManager")) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
library(BiocManager)
cat("[2] Bioconductor 버전:", as.character(BiocManager::version()), "\n\n")

# ─── 방법 1: POSIT RSPM 바이너리 ───────────────────────────────────────────────
cat("[방법 1] POSIT RSPM 바이너리 시도...\n")
tryCatch({
  rspm_url <- "https://packagemanager.posit.co/bioconductor/latest"
  BiocManager::install("bluster", site_repository=rspm_url,
                       update=FALSE, ask=FALSE, type="binary")
  if (is_installed("bluster")) {
    cat("[OK] 방법 1 성공\n"); stop("done")
  }
}, error=function(e) cat("[FAIL] 방법 1 실패:", conditionMessage(e), "\n"))

# ─── 방법 2: Bioconductor 기본 바이너리 ──────────────────────────────────────────
if (!is_installed("bluster")) {
  cat("[방법 2] Bioconductor 기본 바이너리 시도...\n")
  tryCatch({
    BiocManager::install("bluster", update=FALSE, ask=FALSE, type="binary")
    if (is_installed("bluster")) cat("[OK] 방법 2 성공\n")
  }, error=function(e) cat("[FAIL] 방법 2 실패:\n"))
}

# ─── 방법 3: 소스 컴파일 (Rtools 필요) ─────────────────────────────────────
if (!is_installed("bluster")) {
  cat("[방법 3] 소스 컴파일 시도 (Rtools 4.5 필요)...\n")
  tryCatch({
    BiocManager::install("bluster", update=FALSE, ask=FALSE, type="source")
    if (is_installed("bluster")) cat("[OK] 방법 3 성공\n")
  }, error=function(e) cat("[FAIL] 방법 3 실패\n"))
}

# ─── 방법 4: GitHub (최신 소스) ──────────────────────────────────────────────────
if (!is_installed("bluster")) {
  cat("[방법 4] GitHub 소스 시도 (remotes 필요)...\n")
  tryCatch({
    if (!is_installed("remotes")) install.packages("remotes", repos="https://cloud.r-project.org")
    remotes::install_github("bioc/bluster", upgrade="never", build_vignettes=FALSE)
    if (is_installed("bluster")) cat("[OK] 방법 4 성공: GitHub 소스 설치 완료\n")
  }, error=function(e) cat("[FAIL] 방법 4 실패:", conditionMessage(e), "\n"))
}

# ─── 방법 5: 의존성을 먼저 설치 후 재시도 ──────────────────────────────────────
if (!is_installed("bluster")) {
  cat("[방법 5] 의존성 선행 설치 후 재시도...\n")
  deps <- c("BiocNeighbors","BiocParallel","igraph","S4Vectors","BiocGenerics","Matrix","Rcpp")
  for (d in deps) {
    if (!is_installed(d)) {
      tryCatch(BiocManager::install(d, update=FALSE, ask=FALSE, type="binary"),
               error=function(e) install.packages(d, repos="https://cloud.r-project.org", type="binary"))
    }
  }
  tryCatch({
    BiocManager::install("bluster", update=FALSE, ask=FALSE)
    if (is_installed("bluster")) cat("[OK] 방법 5 성공\n")
  }, error=function(e) cat("[FAIL] 방법 5 실패\n"))
}

cat("\n=== 최종 결과 확인 ===\n")
for (p in pkgs_check) {
  status <- if (is_installed(p)) "✓ 설치됨" else "✗ 설치 실패"
  cat(sprintf("    %-30s %s\n", p, status))
}

if (is_installed("bluster")) {
  cat("\n[SUCCESS] bluster 설치 완료!\n")
  cat("이제 curatedMetagenomicData 를 설치하세요:\n")
  cat("  BiocManager::install('curatedMetagenomicData')\n")
} else {
  cat("\n[FAIL] bluster 설치 실패.\n")
  cat("Rtools 4.5 설치: https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html\n")
}
