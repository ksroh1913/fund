# Team 2 — NPS 2027 SAA / MVO / Robust / Policy Black–Litterman

BAF.60080 「연기금 운용전략과 성과평가」 4주차 후속 프로젝트의 재현 가능한 계산 저장소입니다.

## Reproducibility

핵심 분석은 네트워크 없이 재현됩니다.

1. `src/step2_build_corr.py`
   - 저장소에 고정된 2016-09~2026-08 월별 KRW 기준 프록시 수익률을 사용합니다.
   - `data/monthly_returns_krw_part1~4.csv` → 상관행렬 → `Sigma = D_sigma rho D_sigma`.
   - `src/team2_step2_corr_package.py`가 분석용 CSV와 xlsx package를 다시 생성합니다.
2. `src/analyze_steps5_8.py`
3. `src/step9_michaud.py`
4. `src/step10_method_synthesis.py`
5. `src/step11_team_target.py`
6. `src/step12_policy_bl.py`
7. `src/audit_sensitivities.py`
8. `src/step13_stress_test.py`
9. `src/step14_candidate_comparison.py`
10. `src/step14_final_decision.py`

통합 실행은 `submission/team2_NPS_MVO_BL.ipynb`에서 가능합니다. 저장소 root 또는 `submission/`에서 실행해도 root를 자동 탐지합니다.

## Data provenance

고정 월수익률의 원가격 출처는 Yahoo Finance입니다. 사용 프록시는 069500.KS, URTH, EEM, 152380.KS, LQD, PSP, IGF, BIZD, 환율 KRW=X이며, 해외자산은 원화 비헤지 기준으로 환산한 표본입니다. 라이브 Yahoo 다운로드는 결과 재현에 필요하지 않습니다.

2026H1 실제비중과 총자산은 NPS 공개 공시를 이행가능성 예시에 사용합니다. 공식 2027 목표비중은 과제 원문이 직접 제공한 값을 사용합니다.

## CMA

| Asset | mu | sigma |
|---|---:|---:|
| 국내주식 | 6.2% | 19.5% |
| 글로벌 선진국 DM | 5.6% | 16.2% |
| 글로벌 신흥국 EM | 8.4% | 21.0% |
| 국내 국채(초장기) | 3.6% | 6.5% |
| 글로벌 IG 크레딧 | 5.3% | 7.8% |
| 사모주식 PE/VC | 6.8% | 22.0% |
| 실물 인프라 | 6.8% | 11.0% |
| 사모대출 PD | 7.8% | 9.5% |

## Declared Team assumptions

- official mapped 해외주식: DM 92% / EM 8%.
- official mapped 대체투자: PE / 인프라 / PD 동일배분.
- 위 두 규칙은 CMA와 독립적으로 사전에 정한 neutral mapping rule이며 NPS 공식 세부목표가 아닙니다.
- EM, PE, 인프라, PD 상한 10%는 Team detailed concentration guardrail입니다. NPS 공식 세부한도가 아닙니다.
- 현금 0.1%는 baseline에서 유동성·지급 목적의 정책슬리브로 고정합니다. 0~2% 최적화는 sensitivity로 별도 보고합니다.
- gamma=4는 baseline입니다. gamma=2와 4는 약 10.1%의 유사한 정책-MVO 위험수준을 보이며 공식 mapped 약 10.7%와 가깝고, gamma=6은 약 9.8%의 낮은 위험수준 sensitivity입니다.
- BL의 delta=2.5와 MVO의 gamma=4는 역할이 다릅니다. delta는 prior-implied return 역산, gamma는 포트폴리오 효용 위험회피계수입니다.

## Method changes after independent audit

### Michaud
Table B의 Michaud는 **pure Michaud**입니다.
- 월수익률 bootstrap B=300, seed=60080
- CMA 중심 재조정
- 각 draw의 sample correlation × CMA sigma
- 동일 policy constraints

기존의 per-draw Ledoit–Wolf + Michaud는 sensitivity로 유지합니다.

### Ledoit–Wolf comparison
최적화에는 LW correlation × CMA sigma covariance를 사용합니다. 방법론 간 최종 위험비교는 모두 공통 CMA Sigma로 측정하며, LW 자체 covariance 기준 volatility/Sharpe는 보조열로 남깁니다.

### Box / Ellipsoid
- Box baseline: 과제 기본폭(주식 1.0%p, 채권 0.5%p, 대체 1.5%p).
- Box bootstrap-SE 폭: sensitivity only.
- Ellipsoid baseline: proxy bootstrap S_mu.
- Ellipsoid Sigma/T: sensitivity only.

두 uncertainty set의 크기와 구조가 다르므로 Box-vs-Ellipsoid 결과는 set geometry뿐 아니라 calibration에도 영향을 받습니다.

### T assumption
W3 CMA는 forward-looking building-block 전망이라 통계적 sample mean이 아닙니다. T=10은 Omega=Sigma/T를 적용하기 위해 120개월 risk/correlation sample을 10년 equivalent로 준용한 분석가정입니다. T=5/10/20에서 mu_BL, BL-MVO, Robust BL, Stress 1을 모두 다시 계산합니다.

## Stress definitions

### Stress 1
CMA가 최종 BL active 방향에 불리하게 1SE 빗나간다고 보고,
`SE_i = sqrt(Omega_ii)`를 사용합니다. 포트폴리오 비중은 고정하며 CMA에 직접 충격을 적용합니다. posterior를 다시 통과시켜 충격을 20%로 축소하지 않습니다.

### Stress 2
Baseline equity = 국내주식 / DM / EM. PE는 정책 mapping상 대체투자이므로 baseline에서 제외하고, PE 포함 버전을 sensitivity로 보고합니다.

### Stress 3
대체 σ×1.5, 대체-글로벌주식 corr +0.20. raw shocked correlation이 PSD가 아니므로 eigenvalue-clipping과 Higham nearest-correlation 결과를 모두 보고합니다.

## Final IC decision status

최종 IC 안은 **아직 확정하지 않았습니다**. 최종 후보는 다음 세 가지입니다.

- Team: direct CMA view miss에 강함.
- baseline Robust BL (TE 1.0%): equity/correlation 및 alternative-risk shock에 강함.
- 50:50 Team / Robust BL committee overlay: 두 오류원 사이의 절충 후보.

TE 0.95% Robust BL과 25/75 blend는 민감도로만 유지합니다. `src/step14_final_decision.py`의 `FINAL_CANDIDATE=None` 상태에서도 전체 파이프라인이 정상 종료되도록 설계했으며, 최종안을 정할 때 설정값만 변경하고 전체 결과를 다시 생성합니다.

`results/step14_candidate_comparison.csv`에는 CMA 기준 기대수익률과 공식목표 대비 기대수익률 차이, 공통 Sigma 기준 위험, 효용, TE, turnover, Stress 1~3, 국내주식 이행부담을 함께 둡니다. `results/step14_decision_aid.csv`의 손익분기 확률은 시나리오 발생확률을 추정하거나 주장하는 값이 아니라 IC 판단 보조자료입니다.

재심의 조건:
1. 공식 목표 대비 ex-ante TE > 1.0%.
2. 최근 8분기 view hit rate < 55%. 이는 통계적 유의성 기준이 아니라 Team governance monitoring rule입니다.
3. 공식 이행경로 대비 추가 국내주식 거래참여율 > ADV 1%.

이행가능성 계산에서는 **ADV 20조원, 연간 거래일수 240일을 가정**합니다. 전체 2026→2027 이행부담과 공식 목표 대비 추가 부담을 별도로 보고합니다.

## Submission files

- `submission/team2_NPS_MVO_BL.ipynb`
- `submission/team2_weights.csv`
- `submission/team2_views.csv`
- 최종 PPT는 별도 생성·검증 후 제출합니다.

## Checkpoint

감사 수정 전 상태는 commit `3fd5f450f9236461bfde85ded33bc4049d73c5b4`로 보존했습니다.
