# fund

국민연금 CMA 기반 2027 목표비중 및 MVO/Robust/Policy Black-Litterman 종합과제 작업 저장소입니다.

현재 구현 범위: **STEP 2, STEP 5~14**

## 구조

```
src/
  step2_build_corr.py
  analyze_steps5_8.py
  step9_michaud.py
  step10_method_synthesis.py
  step11_team_target.py
  step12_policy_bl.py\n  step13_stress_test.py

results/
  step5_mvo.csv
  step6_lw_diagnostics.csv
  step7_box.csv
  step8_S_mu.csv
  step8_ellipsoid.csv
  step9_michaud_summary.csv
  step9_michaud_meta.json
  step10_tableB_weights.csv
  step10_tableB_metrics.csv
  step10_reaggregated.csv
  step10_mu_50bp_sensitivity_summary.csv
  step10_small_eigenvectors.csv
  step10_meta.json
  step11_tableC_detailed.csv
  step11_tableC_common.csv
  step11_transition_2026H1_to_team.csv
  step11_team_target_meta.json
  step12_tableD.csv
  step12_allocations.csv
  step12_metrics.csv
  step12_common_buckets.csv
  step12_tauSigma.csv
  step12_Omega.csv
  step12_VBL.csv
  step12_meta.json\n  step13_stress_results.csv\n  step13_stress1_view_miss.csv\n  step13_stress1_active_loss.csv\n  step13_covariance_diagnostics.csv\n  step13_relative_to_official.csv\n  step13_meta.json

figures/
  step9_michaud_intervals.svg
```

## 실행 순서

```bash
pip install -r requirements.txt
python src/step2_build_corr.py
python src/analyze_steps5_8.py
python src/step9_michaud.py
python src/step10_method_synthesis.py
python src/step11_team_target.py
python src/step12_policy_bl.py\npython src/step13_stress_test.py
```

## 공통 분석 가정

- 자산: 국내주식, DM, EM, 국내국채, 글로벌 IG, PE/VC, 인프라, PD
- CMA 기대수익률·변동성: 2팀 W3 CMA 유지
- 상관구조: 2016-09~2026-08, 120개월 공개시장 프록시 원화수익률
- 기준 위험회피계수: gamma = 4
- 세부 상한: EM 10%, PE 10%, 인프라 10%, PD 10% (팀 분석 가정)
- 현금 0.1%: 8자산 모델 밖에서 고정

## STEP 11 팀 목표비중 시나리오

| 자산군 | 공식 2027 | 2026.6 실제 | 팀 시나리오 |
|---|---:|---:|---:|
| 국내주식 | 20.8% | 29.1% | 22.0% |
| 해외주식 | 35.6% | 35.4% | 34.0% |
| 국내채권 | 21.8% | 15.4% | 20.5% |
| 해외채권 | 7.4% | 5.9% | 8.5% |
| 대체투자 | 14.3% | 14.0% | 14.9% |
| 단기자금 | 0.1% | 0.2% | 0.1% |

## STEP 12 Policy Black-Litterman

- prior: **2027 official mapped**, Step 11 team target은 prior로 사용하지 않음
- delta = 2.5, tau = 0.025
- CMA는 총수익률이므로 rf=3.0% 차감 후 excess return으로 Q 구성
- P = I
- T = 120개월 / 12 = 10년
- Omega = Sigma / T
- posterior weight: prior 80%, CMA 20%
- Sigma_BL: 기준분석에서는 Sigma 그대로 사용
- BL 배분 제약: 자산별 공식목표 대비 |active| <= 3%p, ex-ante TE <= 1.0%

### BL posterior total expected return

| 자산 | CMA | BL posterior |
|---|---:|---:|
| 국내주식 | 6.20% | 6.80% |
| 글로벌 DM | 5.60% | 6.72% |
| 글로벌 EM | 8.40% | 7.58% |
| 국내국채 | 3.60% | 3.59% |
| 글로벌 IG | 5.30% | 4.31% |
| PE/VC | 6.80% | 7.55% |
| 인프라 | 6.80% | 5.45% |
| PD | 7.80% | 5.25% |

### BL 자산배분

| 자산 | 공식 mapped | BL-MVO | Robust BL |
|---|---:|---:|---:|
| 국내주식 | 20.80% | 18.57% | 19.17% |
| 글로벌 DM | 32.75% | 29.75% | 30.09% |
| 글로벌 EM | 2.85% | 3.68% | 2.75% |
| 국내국채 | 21.80% | 23.84% | 24.80% |
| 글로벌 IG | 7.40% | 10.40% | 9.25% |
| PE/VC | 4.77% | 1.77% | 1.77% |
| 인프라 | 4.77% | 4.13% | 4.32% |
| PD | 4.77% | 7.77% | 7.77% |
| 단기자금 | 0.10% | 0.10% | 0.10% |

BL-MVO가 여러 active cap과 TE 1.0% 한도에 붙어, 수작업으로 수정하지 않고 V_BL을 불확실성 행렬로 사용한 Ellipsoidal Robust BL을 추가 적용했다. 최종 스트레스 분석에서는 **Robust BL**을 BL 최종비중으로 사용한다.

## 현재 상태

- STEP 2: 상관/공분산 구축
- STEP 5: 기준 MVO
- STEP 6: Ledoit-Wolf
- STEP 7: Box Robust
- STEP 8: Ellipsoidal Robust
- STEP 9: Michaud Resampling
- STEP 10: 방법론 종합비교
- STEP 11: 2027 팀 목표비중 시나리오
- STEP 12: Policy Black-Litterman 완료
- STEP 13: Stress Test 완료
- **STEP 14: 최종 의결안·재심의 조건 완료**
- 다음: **최종 노트북 통합 → PPT 작성**

## STEP 14 최종 의결안

최종 분석 시나리오는 Robust BL 비중을 사용한다.

| 자산군 | 공식 2027 | 최종 |
|---|---:|---:|
| 국내주식 | 20.8% | 19.17% |
| 해외주식 | 35.6% | 32.83% |
| 국내채권 | 21.8% | 24.80% |
| 해외채권 | 7.4% | 9.25% |
| 대체투자 | 14.3% | 13.85% |
| 단기자금 | 0.1% | 0.10% |

재심의 조건:
1. 공식목표 대비 사전 TE > 1.0%
2. 최근 8분기 뷰 적중률 < 55%
3. 국내주식 목표 이행에 필요한 시장참여율 > 일평균 거래대금의 1%

제출용 중간파일:
- submission/team2_weights.csv
- submission/team2_views.csv
- docs/step14_final_resolution.md

## STEP 13 Stress Test

### Stress 1 — CMA 뷰 1SE 실패
최종 Robust BL의 active 부호에 불리하도록 각 CMA 기대수익률을 Step 8 bootstrap 표준오차 1개만큼 충격하고, Policy BL의 CMA 가중치 20%를 통해 posterior에 반영했다. 이 충격은 active bet을 일부러 불리하게 만드는 상대성과 스트레스이므로 절대 기대수익률보다 공식목표 대비 active alpha 변화를 본다.

- Team active loss vs official: 약 **-0.65bp**
- Robust BL active loss vs official: 약 **-15.78bp**

### Stress 2 — 주식 기대수익률 -2%p, 주식 상관 +0.15

| 포트폴리오 | 기대수익률 변화 | 변동성 변화 | 효용 변화 |
|---|---:|---:|---:|
| Official | -1.128%p | +0.373%p | -1.291%p |
| Team | -1.120%p | +0.430%p | -1.306%p |
| Robust BL | **-1.040%p** | **+0.348%p** | **-1.179%p** |

### Stress 3 — 대체 변동성 +50%, 대체-글로벌주식 상관 +0.20

| 포트폴리오 | 변동성 | 변동성 변화 | 효용 변화 |
|---|---:|---:|---:|
| Official | 11.49% | +0.778%p | -0.345%p |
| Team | 11.32% | +0.716%p | -0.314%p |
| Robust BL | **10.40%** | **+0.585%p** | **-0.237%p** |

Stress 3의 직접 상관충격 행렬은 최소 고유값이 -0.201로 비정상(비-PSD)이 되어 nearest-PSD 보정을 적용했다. 보정 적용 여부와 최종 최소 고유값은 results/step13_covariance_diagnostics.csv에 기록했다.

세 스트레스 모두에서 Robust BL은 위험 증가 또는 효용 감소 측면에서 공식목표와 팀 시나리오보다 완충력이 높게 나타났다. 단, Stress 1에서는 CMA 뷰 실패가 BL active bet을 겨냥하기 때문에 공식목표 대비 active alpha가 추가로 약 15.8bp 악화되며, 이는 최종 재심의 조건의 뷰 적중률 기준과 연결한다.


## STEP 15 통합 재현 노트북

제출용 노트북:
- `submission/team2_NPS_MVO_BL.ipynb`

노트북은 STEP 2 → 5~14를 순서대로 실행하고 핵심 중간표를 표시한다. 마지막 셀에서 `team2_weights.csv`와 `team2_views.csv`를 다시 생성한다.

다음 단계: 과제 지정 12장 구조의 최종 PPT 작성.
