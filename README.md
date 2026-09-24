# fund

국민연금 CMA 기반 2027 목표비중 및 MVO/Robust/Black-Litterman 종합과제 작업 저장소입니다.

현재 커밋에는 **STEP 2 및 STEP 5~11 재현 코드와 중간 산출물**을 정리했습니다.

## 구조

```
src/
  step2_build_corr.py
  analyze_steps5_8.py
  step9_michaud.py
  step10_method_synthesis.py\n  step11_team_target.py

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
  step10_meta.json\n  step11_tableC_detailed.csv\n  step11_tableC_common.csv\n  step11_transition_2026H1_to_team.csv\n  step11_team_target_meta.json

figures/
  step9_michaud_intervals.svg
```

## 분석 자산 순서

1. 국내주식
2. 글로벌 선진국 DM
3. 글로벌 신흥국 EM
4. 국내 국채(초장기)
5. 글로벌 IG 크레딧
6. 사모주식 PE/VC
7. 실물 인프라
8. 사모대출 PD

## 실행

```bash
pip install -r requirements.txt
python src/step2_build_corr.py
python src/analyze_steps5_8.py
python src/step9_michaud.py
python src/step10_method_synthesis.py\npython src/step11_team_target.py
```

## 주요 분석 가정

- 2팀 3주차 CMA의 기대수익률과 변동성은 유지
- 상관구조는 2016-09~2026-08, 120개월 공개시장 프록시의 원화 기준 월간수익률에서 추정
- 해외주식 공식비중의 세부 매핑: DM 92%, EM 8%
- 대체투자 공식비중의 세부 매핑: PE/VC, 인프라, PD 동일비중
- 기준 위험회피계수: gamma = 4
- 세부 상한: EM 10%, PE 10%, 인프라 10%, PD 10%
  - 공식 국민연금 세부 한도가 아닌 팀 분석 가정
- Box 폭: 주식 1.0%p, 채권 0.5%p, 대체 1.5%p
- Ellipsoid: bootstrap 5,000회, seed 60080
- Michaud: B=300, seed=60080, 매 반복 Ledoit-Wolf 재추정 후 평균

## STEP 10 방법론 종합

### Table B

| 자산 | Official mapped | Policy MVO | LW | Box | Ellipsoid | Michaud |
|---|---:|---:|---:|---:|---:|---:|
| 국내주식 | 20.82% | 21.38% | 19.99% | 20.65% | 15.00% | 20.05% |
| 글로벌 DM | 32.78% | 20.62% | 22.01% | 21.35% | 27.00% | 26.01% |
| 글로벌 EM | 2.85% | 10.00% | 10.00% | 10.00% | 10.00% | 7.75% |
| 국내국채 | 21.82% | 18.00% | 18.00% | 21.65% | 25.00% | 20.08% |
| 글로벌 IG | 7.41% | 12.00% | 12.00% | 12.00% | 10.00% | 10.26% |
| PE/VC | 4.77% | 0.00% | 0.00% | 0.00% | 0.00% | 2.00% |
| 인프라 | 4.77% | 8.00% | 8.00% | 4.35% | 10.00% | 6.34% |
| PD | 4.77% | 10.00% | 10.00% | 10.00% | 3.00% | 7.51% |

### 50bp 기대수익률 민감도

- 국내주식 CMA를 -0.5%p 낮추면 국내주식 비중은 약 **-3.83%p**
- 글로벌 DM CMA를 +0.5%p 높이면 DM 비중은 약 **+3.83%p**
- 나머지 자산은 현재 정책제약/세부상한에 붙어 있어 ±0.5%p 범위에서 자체 비중 변화가 거의 없음
- 따라서 이 결과는 “민감하지 않다”라기보다 **제약이 국소 민감도를 가리고 있음**으로 해석

### 작은 고유값

최소 고유값은 약 0.001493. 최소 고유벡터는 부호를 뒤집어도 동일한 방향이므로, 경제적으로는 대략

- 한쪽: 국내국채·PD·DM
- 반대쪽: 글로벌 IG·PE·인프라

의 상대가치 조합을 뜻한다. 작은 분산으로 추정된 이 상대조합은 역공분산을 사용하는 MVO에서 작은 입력오차를 큰 비중 이동으로 증폭시킬 수 있으며, 공매도 금지·정책제약 아래에서는 상·하한 코너해 형태로 나타난다.

### 주된 참고모형

**Michaud(LW+재표본) 평균을 주된 참고모형으로 사용**한다.

- 단일 표본 MVO의 코너해를 평균화
- μ와 Σ의 표본불확실성을 함께 반영
- 각 반복에서 LW를 적용해 공분산 잡음도 완화
- 공식 mapped 대비 TE와 회전율이 기준 MVO보다 낮음
- Ellipsoid는 기대수익률 오차의 방향성과 크기를 점검하는 강건성 검증모형으로 병행

단, 이후 w2027 team은 Michaud 비중을 그대로 복사하지 않고 규모·시장지분·대체투자 집행속도·환위험·2026→2027 이행가능성을 추가 판단한다.

## 현재 상태

- STEP 2: 상관/공분산 구축
- STEP 5: 기준 MVO
- STEP 6: Ledoit-Wolf
- STEP 7: Box Robust
- STEP 8: Ellipsoidal Robust
- STEP 9: Michaud Resampling
- **STEP 10: 방법론 종합비교 완료**
- 다음: **2027 팀 목표비중 제안 → Policy Black-Litterman → Stress Test**
