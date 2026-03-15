# 마이크로바이옴 기반 대장암 예측 머신러닝 연구

**Statistical Comparison and External Validation of Machine Learning Models for Microbiome-Based Colorectal Cancer Prediction**

🌐 **홈페이지**: [https://mtc-ak.github.io/ICCB](https://mtc-ak.github.io/ICCB)

---

## 개요

장내 마이크로바이옴 샷건 메타게노믹스 데이터를 활용하여 대장암(CRC)을 예측하는 로지스틱 회귀, 랜덤 포레스트, XGBoost 모델을 비교 분석하고 외부 코호트에서 검증한 연구입니다.

| 항목 | 내용 |
|------|------|
| 개발 코호트 | ZellerG_2014 (n=290) |
| 외부 검증 코호트 | YachidaS_2019 (n=212) |
| 데이터 출처 | curatedMetagenomicData (Bioconductor) |
| 최고 성능 모델 | XGBoost (AUC 0.998 개발 / 0.870 외부) |

## 저장소 구조

```
ICCB/
├── index.html
├── README.md
├── code/
│   ├── microbiome_CRC_prediction_full.R
│   └── microbiome_CRC_prediction_full.py
├── figures/
│   ├── Figure1_ROC_curves.png
│   ├── Figure2_Calibration_plots.png
│   ├── Figure3_Performance_barplot.png
│   ├── Figure4_Feature_importance.png
│   └── Figure5_Generalizability.png
├── results/
│   ├── final_results.json
│   ├── metrics_dev_test.tsv
│   ├── metrics_external.tsv
│   ├── ece_summary.tsv
│   └── statistical_testing_summary.txt
└── dashboard/
    └── results_dashboard.html
```

## 핵심 결과

### 개발 코호트 (ZellerG_2014, n=290)

| 모델 | AUC | Accuracy | Recall | F1 |
|------|-----|----------|--------|----||
| LR | 0.997 | 0.739 | 0.465 | 0.635 |
| RF | 0.907 | 0.886 | 0.814 | 0.875 |
| **XGB** | **0.998** | **0.977** | **0.977** | **0.977** |

### 외부 검증 코호트 (YachidaS_2019, n=212)

| 모델 | AUC | Accuracy | Recall | F1 |
|------|-----|----------|--------|----||
| LR | 0.864 | 0.528 | **0.020** ⚠️ | 0.038 |
| RF | 0.802 | 0.764 | 0.539 | 0.688 |
| **XGB** | **0.870** | 0.684 | 0.353 | 0.518 |

## 주요 발견

1. **임계값 일반화 실패**: LR 최적 임계값(0.986)이 외부 코호트에서 Recall=0.020 극단적 저하
2. **DeLong-McNemar 불일치**: AUC 기준 비유의(p=0.811) ↔ 분류 패턴 기준 고도유의(p<0.001)
3. **다중 지표 필요성**: AUC 단독 비교는 실제 임상 성능을 과대평가할 수 있음

## 실행 방법

### R
```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("curatedMetagenomicData")
source("code/microbiome_CRC_prediction_full.R")
```

### Python
```bash
pip install scikit-learn xgboost scipy numpy pandas matplotlib
python code/microbiome_CRC_prediction_full.py
```

## 라이센스
MIT License
- ZellerG_2014: [PubMed 25432777](https://pubmed.ncbi.nlm.nih.gov/25432777/)
- YachidaS_2019: [PubMed 31273294](https://pubmed.ncbi.nlm.nih.gov/31273294/)