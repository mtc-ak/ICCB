"""
============================================================
Microbiome-Based CRC Prediction: Complete Python Pipeline
============================================================
논문: 마이크로바이옴 기반 대장암 예측을 위한
      머신러닝 모델의 통계적 비교 및 외부 검증 연구
------------------------------------------------------------
개발 코호트 : ENA ERP005534 (ZellerG_2014, n=290)
외부 검증   : PRJDB4176    (YachidaS_2019, n=212)
모델        : Logistic Regression, Random Forest, XGBoost
지표        : AUC, PR-AUC, Accuracy, Recall, F1, Brier, ECE
검정        : DeLong (AUC), McNemar (분류 패턴)
Bootstrap   : 2,000 iterations, 95% CI
------------------------------------------------------------
필수 패키지:
  pip install numpy pandas scikit-learn xgboost scipy
              matplotlib seaborn
선택(실제 데이터):
  pip install rpy2   # R curatedMetagenomicData 연동 시
============================================================
"""

# ─────────────────────────────────────────────────────────
# [0]  Import
# ─────────────────────────────────────────────────────────
import os
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path

# Scikit-learn
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    roc_auc_score, average_precision_score,
    accuracy_score, recall_score, precision_score,
    f1_score, brier_score_loss,
    roc_curve, precision_recall_curve,
)

# XGBoost
import xgboost as xgb

# SciPy (통계 검정)
from scipy import stats
from scipy.special import erfc

warnings.filterwarnings("ignore")
np.random.seed(42)

# ─────────────────────────────────────────────────────────
# [1]  설정값 (User Config)
# ─────────────────────────────────────────────────────────
RARE_DEPTH  = 10_000   # Rarefaction depth
PREV_CUT    = 0.05     # Prevalence filter (5%)
N_BINS_CAL  = 10       # Calibration bins
BOOT_N      = 2_000    # Bootstrap iterations
TRAIN_RATIO = 0.70     # Train/Test split
SEED        = 42

OUTDIR = Path("results_crc_python")
OUTDIR.mkdir(parents=True, exist_ok=True)

MODEL_COLORS = {"LR": "#E07B54", "RF": "#5B8DB8", "XGB": "#3A7D44"}

print("=" * 50)
print("  Microbiome CRC Prediction Pipeline (Python)")
print("=" * 50)


# ─────────────────────────────────────────────────────────
# [2]  데이터 로드
#      옵션 A: 시뮬레이션 (기본값)
#      옵션 B: curatedMetagenomicData via rpy2 (주석 해제)
# ─────────────────────────────────────────────────────────

# ── 옵션 B: rpy2로 실제 공개 데이터 로드 ─────────────────
# (rpy2 설치 후 아래 주석 해제)
#
# import rpy2.robjects as ro
# from rpy2.robjects import pandas2ri
# from rpy2.robjects.packages import importr
# pandas2ri.activate()
# cmd_pkg = importr("curatedMetagenomicData")
# tse_pkg = importr("TreeSummarizedExperiment")
#
# def load_cmd_dataset(dataset_name):
#     """curatedMetagenomicData에서 count matrix + 메타데이터 반환"""
#     se = ro.r(f"""
#         curatedMetagenomicData("{dataset_name}",
#                                dryrun=FALSE, counts=TRUE)[[1]]
#     """)
#     counts_r = ro.r("function(x) t(assay(x))")(se)
#     meta_r   = ro.r("function(x) as.data.frame(colData(x))")(se)
#     counts   = pd.DataFrame(
#         np.array(counts_r),
#         index   = list(ro.r("rownames")(counts_r)),
#         columns = list(ro.r("colnames")(counts_r))
#     )
#     meta = pandas2ri.rpy2py(meta_r)
#     return counts, meta
#
# dev_counts, dev_meta = load_cmd_dataset("ZellerG_2014.relative_abundance")
# ext_counts, ext_meta = load_cmd_dataset("YachidaS_2019.relative_abundance")
#
# dev_y = (dev_meta["study_condition"] == "CRC").astype(int).values
# ext_mask = ext_meta["study_condition"].isin(["CRC", "control"])
# ext_counts = ext_counts.loc[ext_mask]
# ext_meta   = ext_meta.loc[ext_mask]
# ext_y = (ext_meta["study_condition"] == "CRC").astype(int).values
#
# ─────────────────────────────────────────────────────────

# ── 옵션 A: 시뮬레이션 데이터 (기본값) ───────────────────
print("\n[2] Simulating realistic microbiome data...")
print("    (실제 데이터: 상단 옵션 B 주석 해제)")

N_TAXA = 4812   # 전처리 전 taxa 수 (논문과 동일)

# CRC / Control 관련 미생물 (문헌 기반)
CRC_TAXA = [
    "Fusobacterium_nucleatum", "Peptostreptococcus_stomatis",
    "Parvimonas_micra", "Gemella_morbillorum", "Clostridium_symbiosum",
    "Bacteroides_fragilis", "Porphyromonas_asaccharolytica",
    "Dialister_pneumosintes", "Peptoniphilus_harei",
    "Lachnoclostridium_sp",
]
CTRL_TAXA = [
    "Faecalibacterium_prausnitzii", "Bifidobacterium_longum",
    "Roseburia_intestinalis", "Blautia_obeum",
    "Lachnospiraceae_bacterium", "Ruminococcus_champanellensis",
    "Coprococcus_eutactus", "Eubacterium_hallii",
    "Butyrivibrio_fibrisolvens",
]


def simulate_cohort(n_crc: int, n_ctrl: int,
                    n_taxa: int = N_TAXA, seed: int = 42) -> tuple:
    """
    Dirichlet-Multinomial 모델 기반 마이크로바이옴 count 시뮬레이션.
    - CRC: CRC-관련 taxa 농도 증가 / Control-관련 taxa 감소
    - Control: 반대 패턴
    반환: (count_matrix [samples × taxa], labels [0/1], taxa_names)
    """
    rng = np.random.default_rng(seed)
    n   = n_crc + n_ctrl
    y   = np.array([1] * n_crc + [0] * n_ctrl)

    alpha_base = rng.exponential(0.01, size=n_taxa).clip(1e-5, 10)
    counts = np.zeros((n, n_taxa), dtype=np.float64)

    crc_idx  = np.arange(10)
    ctrl_idx = np.arange(10, 19)

    for i in range(n):
        alpha = alpha_base.copy()
        if y[i] == 1:
            alpha[crc_idx]  *= rng.uniform(5, 20, len(crc_idx))
            alpha[ctrl_idx] *= rng.uniform(0.1, 0.4, len(ctrl_idx))
        else:
            alpha[ctrl_idx] *= rng.uniform(5, 15, len(ctrl_idx))
            alpha[crc_idx]  *= rng.uniform(0.1, 0.5, len(crc_idx))

        depth = int(rng.integers(15_000, 35_000))
        probs = rng.dirichlet(alpha)
        counts[i] = rng.multinomial(depth, probs)

    # taxa 이름 생성
    taxa_names = [f"taxa_{j}" for j in range(n_taxa)]
    for k, name in enumerate(CRC_TAXA):
        taxa_names[k] = name
    for k, name in enumerate(CTRL_TAXA):
        taxa_names[10 + k] = name

    return counts, y, np.array(taxa_names)


dev_counts_raw, dev_y, taxa_all = simulate_cohort(141, 149, N_TAXA, seed=42)
ext_counts_raw, ext_y, _        = simulate_cohort(102, 110, N_TAXA, seed=123)

print(f"    DEV raw: {dev_counts_raw.shape}  ({dev_y.sum()} CRC)")
print(f"    EXT raw: {ext_counts_raw.shape}  ({ext_y.sum()} CRC)")


# ─────────────────────────────────────────────────────────
# [3]  전처리 함수 정의
# ─────────────────────────────────────────────────────────
print("\n[3] Preprocessing functions defined.")


def safe_rarefy(counts: np.ndarray, depth: int,
                seed: int = 42) -> tuple[np.ndarray, np.ndarray]:
    """
    Rarefaction: depth 미만 샘플 제거 후 다항분포 재샘플링.
    반환: (rarefied_counts, kept_mask)
    """
    rng  = np.random.default_rng(seed)
    rs   = counts.sum(axis=1)
    keep = rs >= depth
    if not keep.all():
        print(f"  [WARN] {(~keep).sum()} samples dropped (reads < {depth})")
    mat = counts[keep].astype(float)
    out = np.zeros_like(mat, dtype=float)
    for i, row in enumerate(mat):
        total = row.sum()
        probs = row / total
        out[i] = rng.multinomial(depth, probs).astype(float)
    return out, keep


def prevalence_filter(counts: np.ndarray,
                      cutoff: float = 0.05) -> np.ndarray:
    """
    Prevalence filter: cutoff 비율 이상의 샘플에서 검출된 taxa만 유지.
    반환: bool mask (taxa 차원)
    """
    prev = (counts > 0).mean(axis=0)
    return prev >= cutoff


def clr_transform(counts: np.ndarray,
                  pseudocount: float = 1.0) -> np.ndarray:
    """
    Centered Log-Ratio 변환.
    CLR(xᵢ) = log(xᵢ + ps) − mean(log(x + ps))
    """
    X    = counts + pseudocount
    logX = np.log(X)
    return logX - logX.mean(axis=1, keepdims=True)


# ─────────────────────────────────────────────────────────
# [4]  전처리 파이프라인 실행
# ─────────────────────────────────────────────────────────
print("\n[4] Running preprocessing pipeline...")

# DEV
dev_rare, dev_keep_samp = safe_rarefy(dev_counts_raw, RARE_DEPTH, SEED)
dev_y2 = dev_y[dev_keep_samp]

taxa_mask  = prevalence_filter(dev_rare, PREV_CUT)
dev_filt   = dev_rare[:, taxa_mask]
dev_clr    = clr_transform(dev_filt)
keep_taxa  = taxa_all[taxa_mask]

print(f"  DEV: {dev_clr.shape[0]} samples × {dev_clr.shape[1]} taxa")

# EXT (DEV 기준 taxa 집합 정렬)
ext_rare, ext_keep_samp = safe_rarefy(ext_counts_raw, RARE_DEPTH, SEED + 1)
ext_y2 = ext_y[ext_keep_samp]

# taxa 정렬: DEV와 같은 taxa 집합 사용
n_keep = taxa_mask.sum()
ext_filt = ext_rare[:, :n_keep] if ext_rare.shape[1] >= n_keep else \
    np.hstack([ext_rare,
               np.zeros((ext_rare.shape[0], n_keep - ext_rare.shape[1]))])
ext_filt = ext_filt[:, :n_keep]
ext_clr  = clr_transform(ext_filt)

print(f"  EXT: {ext_clr.shape[0]} samples × {ext_clr.shape[1]} taxa")


# ─────────────────────────────────────────────────────────
# [5]  Stratified Train / Test 분할 (70/30)
# ─────────────────────────────────────────────────────────
print("\n[5] Stratified train/test split (70/30)...")


def stratified_split(X: np.ndarray, y: np.ndarray,
                     ratio: float = 0.7,
                     seed: int = 42) -> tuple:
    rng      = np.random.default_rng(seed)
    pos_idx  = np.where(y == 1)[0]
    neg_idx  = np.where(y == 0)[0]
    rng.shuffle(pos_idx); rng.shuffle(neg_idx)

    n_pos_tr = int(len(pos_idx) * ratio)
    n_neg_tr = int(len(neg_idx) * ratio)

    tr = np.sort(np.concatenate([pos_idx[:n_pos_tr], neg_idx[:n_neg_tr]]))
    te = np.sort(np.setdiff1d(np.arange(len(y)), tr))
    return X[tr], y[tr], X[te], y[te]


X_train, y_train, X_test, y_test = stratified_split(
    dev_clr, dev_y2, TRAIN_RATIO, SEED
)
print(f"  Train: {len(y_train)} ({y_train.sum()} CRC) | "
      f"Test: {len(y_test)} ({y_test.sum()} CRC)")


# ─────────────────────────────────────────────────────────
# [6]  모델 학습
# ─────────────────────────────────────────────────────────
print("\n[6] Training models...")

# 6-1. Logistic Regression (L2, solver=lbfgs)
print("  [6.1] Logistic Regression")
lr_model = LogisticRegression(
    C=1.0, penalty="l2", solver="lbfgs",
    max_iter=1000, random_state=SEED
)
lr_model.fit(X_train, y_train)
lr_p_train = lr_model.predict_proba(X_train)[:, 1]
lr_p_test  = lr_model.predict_proba(X_test)[:, 1]
lr_p_ext   = lr_model.predict_proba(ext_clr)[:, 1]

# 6-2. Random Forest (500 trees)
print("  [6.2] Random Forest (n_estimators=500)")
rf_model = RandomForestClassifier(
    n_estimators=500,
    max_features="sqrt",
    min_samples_leaf=1,
    n_jobs=-1,
    random_state=SEED
)
rf_model.fit(X_train, y_train)
rf_p_train = rf_model.predict_proba(X_train)[:, 1]
rf_p_test  = rf_model.predict_proba(X_test)[:, 1]
rf_p_ext   = rf_model.predict_proba(ext_clr)[:, 1]

# 6-3. XGBoost
print("  [6.3] XGBoost (max_depth=6, eta=0.1, n=300)")
xgb_model = xgb.XGBClassifier(
    n_estimators=300,
    max_depth=6,
    learning_rate=0.1,
    subsample=0.8,
    colsample_bytree=0.8,
    min_child_weight=1,
    objective="binary:logistic",
    eval_metric="auc",
    use_label_encoder=False,
    random_state=SEED,
    verbosity=0,
)
xgb_model.fit(X_train, y_train,
              eval_set=[(X_test, y_test)],
              verbose=False)
xgb_p_train = xgb_model.predict_proba(X_train)[:, 1]
xgb_p_test  = xgb_model.predict_proba(X_test)[:, 1]
xgb_p_ext   = xgb_model.predict_proba(ext_clr)[:, 1]
print("  Models trained.")


# ─────────────────────────────────────────────────────────
# [7]  임계값 선택 (Youden J — 훈련셋 기준)
# ─────────────────────────────────────────────────────────
print("\n[7] Selecting thresholds (Youden J on train set)...")


def youden_threshold(y_true: np.ndarray, p_hat: np.ndarray) -> float:
    """ROC 기반 Youden J 최대화 임계값."""
    fpr, tpr, thresholds = roc_curve(y_true, p_hat)
    j_scores = tpr - fpr
    best_idx = np.argmax(j_scores)
    return float(thresholds[best_idx])


lr_thr  = youden_threshold(y_train, lr_p_train)
rf_thr  = youden_threshold(y_train, rf_p_train)
xgb_thr = youden_threshold(y_train, xgb_p_train)

print(f"  LR thr={lr_thr:.4f} | RF thr={rf_thr:.4f} | XGB thr={xgb_thr:.4f}")


# ─────────────────────────────────────────────────────────
# [8]  평가 함수 정의
# ─────────────────────────────────────────────────────────

def calibration_ece(y_true: np.ndarray, p_hat: np.ndarray,
                    n_bins: int = 10) -> tuple[pd.DataFrame, float]:
    """
    ECE (Expected Calibration Error) + Reliability 테이블 반환.
    """
    bins    = np.linspace(0, 1, n_bins + 1)
    p_means, y_rates, ns = [], [], []

    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (p_hat >= lo) & (p_hat <= hi) if hi == 1 else \
               (p_hat >= lo) & (p_hat <  hi)
        if mask.sum() == 0:
            continue
        p_means.append(p_hat[mask].mean())
        y_rates.append(y_true[mask].mean())
        ns.append(mask.sum())

    p_means = np.array(p_means)
    y_rates = np.array(y_rates)
    weights = np.array(ns) / len(y_true)
    ece     = float(np.sum(weights * np.abs(y_rates - p_means)))

    table = pd.DataFrame({
        "p_mean": p_means,
        "y_rate": y_rates,
        "n":      np.array(ns),
        "weight": weights,
    })
    return table, ece


def bootstrap_auc_ci(y_true: np.ndarray, p_hat: np.ndarray,
                     n_boot: int = 2000,
                     alpha: float = 0.05,
                     seed: int = 42) -> tuple[float, float]:
    """
    Bootstrap으로 AUC 95% 신뢰구간 산출.
    """
    rng  = np.random.default_rng(seed)
    n    = len(y_true)
    aucs = []
    for _ in range(n_boot):
        idx  = rng.integers(0, n, size=n)
        ys   = y_true[idx]
        ps   = p_hat[idx]
        if len(np.unique(ys)) < 2:
            continue
        aucs.append(roc_auc_score(ys, ps))
    aucs = np.array(aucs)
    return (float(np.percentile(aucs, 100 * alpha / 2)),
            float(np.percentile(aucs, 100 * (1 - alpha / 2))))


def evaluate_model(name: str,
                   y_true: np.ndarray, p_hat: np.ndarray,
                   threshold: float,
                   boot_n: int = 2000) -> dict:
    """단일 모델 전체 지표 산출."""
    auc_val = roc_auc_score(y_true, p_hat)
    pr_auc  = average_precision_score(y_true, p_hat)

    ci_lo, ci_hi = bootstrap_auc_ci(y_true, p_hat, boot_n) \
                   if boot_n > 0 else (np.nan, np.nan)

    y_pred = (p_hat >= threshold).astype(int)
    acc    = accuracy_score(y_true, y_pred)
    rec    = recall_score(y_true, y_pred, zero_division=0)
    prec   = precision_score(y_true, y_pred, zero_division=0)
    f1     = f1_score(y_true, y_pred, zero_division=0)
    brier  = brier_score_loss(y_true, p_hat)

    _, ece = calibration_ece(y_true, p_hat, N_BINS_CAL)

    return {
        "Model":     name,
        "AUC":       round(auc_val, 4),
        "AUC_CI_L":  round(ci_lo,   4),
        "AUC_CI_U":  round(ci_hi,   4),
        "PR_AUC":    round(pr_auc,  4),
        "Accuracy":  round(acc,     4),
        "Recall":    round(rec,     4),
        "Precision": round(prec,    4),
        "F1":        round(f1,      4),
        "Brier":     round(brier,   4),
        "ECE":       round(ece,     4),
        "Threshold": round(threshold, 4),
    }


# ─────────────────────────────────────────────────────────
# [9]  성능 평가 실행
# ─────────────────────────────────────────────────────────
print("\n[9] Evaluating models (bootstrap CI = 2,000 iters)...")

rows_dev, rows_ext = [], []
for name, p_te, p_ex, thr in [
    ("LR",  lr_p_test,  lr_p_ext,  lr_thr),
    ("RF",  rf_p_test,  rf_p_ext,  rf_thr),
    ("XGB", xgb_p_test, xgb_p_ext, xgb_thr),
]:
    print(f"  {name}...")
    rows_dev.append(evaluate_model(name, y_test, p_te, thr, BOOT_N))
    rows_ext.append(evaluate_model(name, ext_y2, p_ex, thr, BOOT_N))

res_dev = pd.DataFrame(rows_dev)
res_ext = pd.DataFrame(rows_ext)

SHOW_COLS = ["Model","AUC","AUC_CI_L","AUC_CI_U",
             "PR_AUC","Recall","F1","Brier","ECE"]
print("\n=== Development Cohort Test Results ===")
print(res_dev[SHOW_COLS].to_string(index=False))
print("\n=== External Validation Results ===")
print(res_ext[SHOW_COLS].to_string(index=False))

res_dev.to_csv(OUTDIR / "metrics_dev_test.tsv", sep="\t", index=False)
res_ext.to_csv(OUTDIR / "metrics_external.tsv", sep="\t", index=False)
print("\n  Metrics saved.")


# ─────────────────────────────────────────────────────────
# [10]  통계 검정
# ─────────────────────────────────────────────────────────
print("\n[10] Statistical testing...")


# ── DeLong 검정 (FastDeLong 근사) ─────────────────────────
def delong_test(y_true: np.ndarray,
                score1: np.ndarray,
                score2: np.ndarray) -> tuple[float, float]:
    """
    DeLong et al. (1988) 기반 두 상관 AUC 비교.
    반환: (AUC 차이, p-value)
    """
    pos_idx = np.where(y_true == 1)[0]
    neg_idx = np.where(y_true == 0)[0]
    n_pos, n_neg = len(pos_idx), len(neg_idx)

    def placement_values(scores):
        s_pos = scores[pos_idx]
        s_neg = scores[neg_idx]
        # V10: 각 양성 샘플이 음성 샘플보다 높을 확률
        v10 = np.array([
            np.mean(s_neg < sp) + 0.5 * np.mean(s_neg == sp)
            for sp in s_pos
        ])
        # V01: 각 음성 샘플이 양성 샘플보다 낮을 확률
        v01 = np.array([
            np.mean(s_pos > sn) + 0.5 * np.mean(s_pos == sn)
            for sn in s_neg
        ])
        return v10, v01

    v10_1, v01_1 = placement_values(score1)
    v10_2, v01_2 = placement_values(score2)

    auc1, auc2 = v10_1.mean(), v10_2.mean()

    # 공분산 행렬
    S10 = np.cov(np.stack([v10_1, v10_2])) / n_pos
    S01 = np.cov(np.stack([v01_1, v01_2])) / n_neg
    S   = S10 + S01                             # 2×2

    var_diff = S[0, 0] + S[1, 1] - 2 * S[0, 1]
    if var_diff <= 0:
        return auc1 - auc2, 1.0

    z    = (auc1 - auc2) / np.sqrt(var_diff)
    pval = float(erfc(abs(z) / np.sqrt(2)))     # two-tailed
    return float(auc1 - auc2), pval


# ── McNemar 검정 ──────────────────────────────────────────
def mcnemar_test(y_true: np.ndarray,
                 pred1: np.ndarray,
                 pred2: np.ndarray) -> tuple[float, int, int]:
    """
    McNemar 검정: 두 분류기의 예측 불일치 비교.
    반환: (p-value, b, c)
      b = model1 맞고 model2 틀린 수
      c = model1 틀리고 model2 맞는 수
    """
    b = int(((pred1 == y_true) & (pred2 != y_true)).sum())
    c = int(((pred1 != y_true) & (pred2 == y_true)).sum())
    if b + c == 0:
        return 1.0, b, c
    # 연속성 보정 카이제곱
    chi2 = (abs(b - c) - 1) ** 2 / (b + c)
    pval = float(stats.chi2.sf(chi2, df=1))
    return pval, b, c


# 예측 레이블
pred_lr  = (lr_p_test  >= lr_thr).astype(int)
pred_rf  = (rf_p_test  >= rf_thr).astype(int)
pred_xgb = (xgb_p_test >= xgb_thr).astype(int)

# DeLong
dl_xlr_dev_diff, dl_xlr_dev_p = delong_test(y_test, xgb_p_test, lr_p_test)
dl_xrf_dev_diff, dl_xrf_dev_p = delong_test(y_test, xgb_p_test, rf_p_test)
dl_xlr_ext_diff, dl_xlr_ext_p = delong_test(ext_y2, xgb_p_ext,  lr_p_ext)

# McNemar
mcn_xlr_p, mcn_b_xlr, mcn_c_xlr = mcnemar_test(y_test, pred_xgb, pred_lr)
mcn_xrf_p, mcn_b_xrf, mcn_c_xrf = mcnemar_test(y_test, pred_xgb, pred_rf)

print("\n--- DeLong Tests ---")
print(f"  DEV XGB vs LR: ΔAUC={dl_xlr_dev_diff:+.4f}, p={dl_xlr_dev_p:.4f}")
print(f"  DEV XGB vs RF: ΔAUC={dl_xrf_dev_diff:+.4f}, p={dl_xrf_dev_p:.4f}")
print(f"  EXT XGB vs LR: ΔAUC={dl_xlr_ext_diff:+.4f}, p={dl_xlr_ext_p:.4f}")

print("\n--- McNemar Tests (DEV test set) ---")
print(f"  XGB vs LR: p={mcn_xlr_p:.4f}  (b={mcn_b_xlr}, c={mcn_c_xlr})")
print(f"  XGB vs RF: p={mcn_xrf_p:.4f}  (b={mcn_b_xrf}, c={mcn_c_xrf})")

# 저장
stats_out = {
    "DeLong_DEV_XGB_vs_LR": {"delta_auc": dl_xlr_dev_diff, "p_value": dl_xlr_dev_p},
    "DeLong_DEV_XGB_vs_RF": {"delta_auc": dl_xrf_dev_diff, "p_value": dl_xrf_dev_p},
    "DeLong_EXT_XGB_vs_LR": {"delta_auc": dl_xlr_ext_diff, "p_value": dl_xlr_ext_p},
    "McNemar_DEV_XGB_vs_LR": {"p_value": mcn_xlr_p, "b": mcn_b_xlr, "c": mcn_c_xlr},
    "McNemar_DEV_XGB_vs_RF": {"p_value": mcn_xrf_p, "b": mcn_b_xrf, "c": mcn_c_xrf},
}
with open(OUTDIR / "statistical_testing_summary.json", "w") as f:
    json.dump(stats_out, f, indent=2)
print("\n  Statistical tests saved.")


# ─────────────────────────────────────────────────────────
# [11]  Calibration 분석
# ─────────────────────────────────────────────────────────
print("\n[11] Calibration analysis...")

ece_rows = []
for name, p_te, p_ex in [
    ("LR",  lr_p_test,  lr_p_ext),
    ("RF",  rf_p_test,  rf_p_ext),
    ("XGB", xgb_p_test, xgb_p_ext),
]:
    _, ece_d = calibration_ece(y_test, p_te, N_BINS_CAL)
    _, ece_e = calibration_ece(ext_y2, p_ex, N_BINS_CAL)
    ece_rows.append({"Model": name, "ECE_DEV": round(ece_d, 4),
                     "ECE_EXT": round(ece_e, 4)})

ece_df = pd.DataFrame(ece_rows)
ece_df.to_csv(OUTDIR / "ece_summary.tsv", sep="\t", index=False)
print(ece_df.to_string(index=False))


# ─────────────────────────────────────────────────────────
# [12]  시각화
# ─────────────────────────────────────────────────────────
print("\n[12] Generating plots...")

plt.rcParams.update({
    "font.family":       "DejaVu Sans",
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "axes.grid":         True,
    "grid.alpha":        0.3,
    "figure.dpi":        150,
    "font.size":         11,
})


# ── Figure 1: ROC Curves (DEV + EXT) ──────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

for ax, (y_true, probs_dict), title in [
    (axes[0],
     (y_test, {"LR": lr_p_test, "RF": rf_p_test, "XGB": xgb_p_test}),
     "Development Cohort (Test Set, n=88)"),
    (axes[1],
     (ext_y2, {"LR": lr_p_ext,  "RF": rf_p_ext,  "XGB": xgb_p_ext}),
     "External Validation (n=212)"),
]:
    ax.plot([0, 1], [0, 1], "k--", alpha=0.4, lw=1)
    for m, p in probs_dict.items():
        fpr, tpr, _ = roc_curve(y_true, p)
        auc_val = roc_auc_score(y_true, p)
        ax.plot(fpr, tpr, color=MODEL_COLORS[m], lw=2.2,
                label=f"{m}  AUC={auc_val:.3f}")
    ax.set(xlim=(-0.01, 1.01), ylim=(-0.01, 1.01),
           xlabel="False Positive Rate (1 – Specificity)",
           ylabel="True Positive Rate (Sensitivity)",
           title=title)
    ax.set_aspect("equal")
    ax.legend(loc="lower right", fontsize=9.5)
    ax.set_title(title, fontweight="bold")

fig.suptitle("ROC Curve Comparison: LR vs RF vs XGBoost",
             fontsize=14, fontweight="bold", y=1.01)
plt.tight_layout()
fig.savefig(OUTDIR / "Figure1_ROC_curves.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("  Figure 1 (ROC) saved.")


# ── Figure 2: Calibration Reliability Curves ──────────────
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
model_pairs = [
    (lr_p_test,  "LR",  y_test, 0, 0),
    (rf_p_test,  "RF",  y_test, 0, 1),
    (xgb_p_test, "XGB", y_test, 0, 2),
    (lr_p_ext,   "LR",  ext_y2, 1, 0),
    (rf_p_ext,   "RF",  ext_y2, 1, 1),
    (xgb_p_ext,  "XGB", ext_y2, 1, 2),
]
row_labels = ["Development Cohort (Test)", "External Validation"]

for p_hat, model_name, y_true, row, col in model_pairs:
    ax = axes[row, col]
    cal_tbl, ece_val = calibration_ece(y_true, p_hat, N_BINS_CAL)
    ax.plot([0, 1], [0, 1], "k--", alpha=0.5, lw=1.2)
    ax.plot(cal_tbl["p_mean"], cal_tbl["y_rate"],
            "o-", color=MODEL_COLORS[model_name],
            lw=2, ms=7, mfc="white", mew=2,
            label=f"ECE={ece_val:.3f}")
    ax.fill_between(cal_tbl["p_mean"], cal_tbl["y_rate"],
                    cal_tbl["p_mean"],
                    where=cal_tbl["y_rate"] < cal_tbl["p_mean"],
                    alpha=0.12, color="red")
    ax.fill_between(cal_tbl["p_mean"], cal_tbl["y_rate"],
                    cal_tbl["p_mean"],
                    where=cal_tbl["y_rate"] >= cal_tbl["p_mean"],
                    alpha=0.12, color="steelblue")
    ax.set(xlim=(0, 1), ylim=(0, 1),
           xlabel="Mean Predicted Probability",
           ylabel="Observed Event Rate",
           title=f"{model_name} – {row_labels[row]}\nECE={ece_val:.3f}")
    ax.set_aspect("equal")
    ax.legend(loc="upper left", fontsize=9)
    ax.set_title(f"{model_name} – {row_labels[row]}\nECE={ece_val:.3f}",
                 fontweight="bold", fontsize=10)

fig.suptitle("Calibration Analysis: Reliability Curves",
             fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(OUTDIR / "Figure2_Calibration.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("  Figure 2 (Calibration) saved.")


# ── Figure 3: Multi-metric Bar Chart ──────────────────────
metrics_keys   = ["AUC", "PR_AUC", "Recall", "F1", "Accuracy"]
metrics_labels = ["AUC", "PR-AUC", "Recall", "F1", "Accuracy"]

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
x  = np.arange(len(metrics_keys))
bw = 0.25

for ax, res_df, title in [
    (axes[0], res_dev, "Development Cohort (Test Set)"),
    (axes[1], res_ext, "External Validation Cohort"),
]:
    for i, (model, offset) in enumerate(
            zip(["LR", "RF", "XGB"], [-bw, 0, bw])):
        row  = res_df[res_df["Model"] == model].iloc[0]
        vals = [row[k] for k in metrics_keys]
        bars = ax.bar(x + offset, vals, bw * 0.9,
                      label=model, color=MODEL_COLORS[model],
                      alpha=0.85, edgecolor="white")
        for bar, v in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.005,
                    f"{v:.3f}", ha="center", va="bottom",
                    fontsize=7.5, fontweight="bold")
    # AUC CI
    for model, offset in zip(["LR", "RF", "XGB"], [-bw, 0, bw]):
        row   = res_df[res_df["Model"] == model].iloc[0]
        ci_lo = row["AUC"] - row["AUC_CI_L"]
        ci_hi = row["AUC_CI_U"] - row["AUC"]
        ax.errorbar(0 + offset, row["AUC"],
                    yerr=[[ci_lo], [ci_hi]],
                    fmt="none", color="black",
                    capsize=4, capthick=1.5, lw=1.5)
    ax.set_xticks(x)
    ax.set_xticklabels(metrics_labels, fontsize=11)
    ax.set_ylim([0.0, 1.1])
    ax.set_ylabel("Score", fontsize=12)
    ax.set_title(title, fontweight="bold", fontsize=13)
    ax.legend(loc="upper right", fontsize=10)

fig.suptitle("Model Performance Comparison",
             fontsize=14, fontweight="bold")
plt.tight_layout()
fig.savefig(OUTDIR / "Figure3_Performance.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("  Figure 3 (Performance bar chart) saved.")


# ── Figure 4: Feature Importance ──────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(14, 7))

## RF: Mean Decrease Impurity
rf_imp    = pd.Series(rf_model.feature_importances_,
                      index=keep_taxa)
top_rf    = rf_imp.nlargest(20).sort_values()
axes[0].barh(top_rf.index, top_rf.values,
             color=MODEL_COLORS["RF"], alpha=0.85)
axes[0].set_xlabel("Mean Decrease Impurity", fontsize=11)
axes[0].set_title("Random Forest\nFeature Importance (Top 20)",
                  fontweight="bold", fontsize=11)

## XGB: Gain
xgb_fi  = xgb_model.get_booster().get_score(importance_type="gain")
xgb_imp = (pd.Series(xgb_fi)
           .rename(index=lambda i: keep_taxa[int(i[1:])]
                   if i.startswith("f") else i)
           .nlargest(20).sort_values())
axes[1].barh(xgb_imp.index, xgb_imp.values,
             color=MODEL_COLORS["XGB"], alpha=0.85)
axes[1].set_xlabel("Gain", fontsize=11)
axes[1].set_title("XGBoost\nFeature Importance (Top 20, Gain)",
                  fontweight="bold", fontsize=11)

for ax in axes:
    ax.tick_params(axis="y", labelsize=9)

fig.suptitle("Top 20 Microbiome Feature Importance",
             fontsize=13, fontweight="bold")
plt.tight_layout()
fig.savefig(OUTDIR / "Figure4_Feature_importance.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("  Figure 4 (Feature importance) saved.")


# ── Figure 5: AUC Generalizability (DEV → EXT) ────────────
models  = ["LR", "RF", "XGB"]
auc_d   = [res_dev[res_dev["Model"] == m]["AUC"].values[0]     for m in models]
auc_e   = [res_ext[res_ext["Model"] == m]["AUC"].values[0]     for m in models]
ci_lo_d = [res_dev[res_dev["Model"] == m]["AUC_CI_L"].values[0] for m in models]
ci_hi_d = [res_dev[res_dev["Model"] == m]["AUC_CI_U"].values[0] for m in models]
ci_lo_e = [res_ext[res_ext["Model"] == m]["AUC_CI_L"].values[0] for m in models]
ci_hi_e = [res_ext[res_ext["Model"] == m]["AUC_CI_U"].values[0] for m in models]

fig, ax = plt.subplots(figsize=(8, 5))
x  = np.array([0, 1.3, 2.6])
bw = 0.45
for i, m in enumerate(models):
    ax.bar(x[i] - bw/2, auc_d[i], bw, color=MODEL_COLORS[m],
           alpha=0.9, label=f"{m} (Dev)")
    ax.bar(x[i] + bw/2, auc_e[i], bw, color=MODEL_COLORS[m],
           alpha=0.45, edgecolor=MODEL_COLORS[m], lw=1.5)
    ax.errorbar(x[i] - bw/2, auc_d[i],
                yerr=[[auc_d[i]-ci_lo_d[i]], [ci_hi_d[i]-auc_d[i]]],
                fmt="none", color="black", capsize=4, capthick=1.5)
    ax.errorbar(x[i] + bw/2, auc_e[i],
                yerr=[[auc_e[i]-ci_lo_e[i]], [ci_hi_e[i]-auc_e[i]]],
                fmt="none", color="black", capsize=4, capthick=1.5)
    drop = auc_d[i] - auc_e[i]
    ax.text(x[i], max(auc_d[i], auc_e[i]) + 0.012,
            f"Δ={drop:.3f}", ha="center", fontsize=9,
            color="gray", style="italic")

ax.set_xticks(x)
ax.set_xticklabels(models, fontsize=12)
ax.set_ylim([0.65, 1.05])
ax.set_ylabel("AUC (95% CI)", fontsize=12)
ax.set_title("Model Generalizability: Development → External Validation",
             fontweight="bold", fontsize=12)
dev_patch  = mpatches.Patch(facecolor="#888", alpha=0.9,
                             label="Development cohort")
ext_patch  = mpatches.Patch(facecolor="#888", alpha=0.4,
                             label="External validation")
ax.legend(handles=[dev_patch, ext_patch], fontsize=10)
plt.tight_layout()
fig.savefig(OUTDIR / "Figure5_Generalizability.png",
            dpi=300, bbox_inches="tight")
plt.close()
print("  Figure 5 (Generalizability) saved.")


# ─────────────────────────────────────────────────────────
# [13]  최종 요약 출력
# ─────────────────────────────────────────────────────────
print("\n" + "=" * 50)
print("  FINAL SUMMARY")
print("=" * 50)

print("\n[Development Cohort Test Set]")
print(res_dev[SHOW_COLS].to_string(index=False))
print("\n[External Validation Cohort]")
print(res_ext[SHOW_COLS].to_string(index=False))

print("\n[Statistical Tests]")
print(f"  DeLong  DEV XGB vs LR: p={dl_xlr_dev_p:.4f}")
print(f"  DeLong  DEV XGB vs RF: p={dl_xrf_dev_p:.4f}")
print(f"  DeLong  EXT XGB vs LR: p={dl_xlr_ext_p:.4f}")
print(f"  McNemar DEV XGB vs LR: p={mcn_xlr_p:.4f}  "
      f"(b={mcn_b_xlr}, c={mcn_c_xlr})")
print(f"  McNemar DEV XGB vs RF: p={mcn_xrf_p:.4f}  "
      f"(b={mcn_b_xrf}, c={mcn_c_xrf})")

print(f"\nAll outputs saved in: {OUTDIR.resolve()}")
print("=" * 50)

# ─────────────────────────────────────────────────────────
# END OF SCRIPT
# ─────────────────────────────────────────────────────────
