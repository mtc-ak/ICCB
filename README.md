# 마이크로바이옴 기반 대장암 예측 머신러닝 연구

**Statistical Comparison and External Validation of Machine Learning Models for Microbiome-Based Colorectal Cancer Prediction**

🌐 **연구 홈페이지**: [https://mtc-ak.github.io/ICCB](https://mtc-ak.github.io/ICCB)

---

## 개요

장내 마이크로바이옴 샷건 메타게노믹스 데이터를 활용하여 대장암(CRC)을 예측하는 로지스틱 회귀, 랜덤 포레스트, XGBoost 모델을 비교 분석하고 외부 코호트에서 검증한 연구입니다.

| 항목 | 내용 |
|------|------|
| 개발 코호트 | ZellerG_2014 (n=290) |
| 외부 검증 코호트 | YachidaS_2019 (n=212) |
| 데이터 출처 | curatedMetagenomicData (Bioconductor) |
| 최고 성능 모델 | XGBoost (AUC 0.998 개발 / 0.870 외부) |

---

## 저장소 구조

이 저장소(`ICCB`)는 **공개 홈페이지 전용**입니다.

```
ICCB/ (public - 홈페이지)
└── index.html    # GitHub Pages 연구 소개 홈페이지
```

분석 코드, 결과 데이터, 그래프 파일은 **비공개 저장소(`ICCB-data`)**에서 관리됩니다.

---

## 핵심 결과

### 개발 코호트 (ZellerG_2014, n=290)

| 모델 | AUC | Accuracy | Recall | F1 |
|------|-----|----------|--------|----|n| LR   | 0.997 | 0.739 | 0.465 | 0.635 |
| RF   | 0.907 | 0.886 | 0.814 | 0.875 |
| **XGB** | **0.998** | **0.977** | **0.977** | **0.977** |

### 외부 검증 코호트 (YachidaS_2019, n=212)

| 모델 | AUC | Accuracy | Recall | F1 |
|------|-----|----------|--------|----|n| LR   | 0.864 | 0.528 | **0.020** ⚠️ | 0.038 |
| RF   | 0.802 | 0.764 | 0.539 | 0.688 |
| **XGB** | **0.870** | 0.684 | 0.353 | 0.518 |

---

## 코드 & 데이터 접근

분석 코드 및 결과 데이터 접근이 필요하신 경우, [GitHub 이슈](https://github.com/mtc-ak/ICCB/issues)로 요청해주세요.

---

## 데이터 출처

- ZellerG_2014: [PubMed 25432777](https://pubmed.ncbi.nlm.nih.gov/25432777/)
- YachidaS_2019: [PubMed 31273294](https://pubmed.ncbi.nlm.nih.gov/31273294/)
- curatedMetagenomicData: [Bioconductor](https://bioconductor.org/packages/curatedMetagenomicData/)