# fund

국민연금 CMA 기반 2027 목표비중 및 MVO/Robust/Black-Litterman 종합과제 작업 저장소입니다.

현재 커밋에는 **STEP 2 및 STEP 5~9 재현 코드와 중간 산출물**을 정리했습니다.

## 구조

```
src/
  step2_build_corr.py       # 8자산 실제 프록시 월간수익률 → 상관/공분산 구축
  analyze_steps5_8.py       # 기준 MVO, Ledoit-Wolf, Box, Ellipsoid 재현
  step9_michaud.py          # Michaud 300회 재표본 + LW

results/
  step5_mvo.csv
  step6_lw_diagnostics.csv
  step7_box.csv
  step8_S_mu.csv
  step8_ellipsoid.csv
  step9_michaud_summary.csv
  step9_michaud_meta.json

figures/
  step9_michaud_intervals.svg

requirements.txt
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
```

## 주요 분석 가정

- 2팀 3주차 CMA의 기대수익률과 변동성은 유지
- 상관구조는 2016-09~2026-08, 120개월 공개시장 프록시의 원화 기준 월간수익률에서 추정
- 해외주식 공식비중의 세부 매핑: DM 92%, EM 8%
- 대체투자 공식비중의 세부 매핑: PE/VC, 인프라, PD 동일비중
- 기준 위험회피계수: gamma = 4
- 세부 상한: EM 10%, PE 10%, 인프라 10%, PD 10%
  - 이 한도는 공식 국민연금 한도가 아니라 **팀 분석 가정**이며 추후 민감도 검토 대상
- Box 폭: 주식 1.0%p, 채권 0.5%p, 대체 1.5%p
- Ellipsoid: bootstrap 5,000회, seed 60080
- Michaud:
  - B=300, seed=60080
  - T=120개월을 복원추출
  - 각 표본의 역사적 평균 오차를 2팀 CMA 기대수익률 중심으로 재중심화
  - 각 반복에서 Ledoit-Wolf 상관구조를 재추정하고 CMA 변동성으로 재스케일
  - gamma=4 및 동일 정책제약으로 300회 최적화 후 평균

## STEP 9 Michaud 결과

| 자산 | 기준 MVO | Michaud 평균 | 5% | 50% | 95% |
|---|---:|---:|---:|---:|---:|
| 국내주식 | 21.38% | 20.05% | 15.0% | 22.0% | 25.0% |
| 글로벌 DM | 20.62% | 26.01% | 20.0% | 27.0% | 37.05% |
| 글로벌 EM | 10.00% | 7.75% | 0.0% | 10.0% | 10.0% |
| 국내국채 | 18.00% | 20.08% | 15.0% | 20.16% | 25.0% |
| 글로벌 IG | 12.00% | 10.26% | 5.0% | 12.0% | 12.0% |
| PE/VC | 0.00% | 2.00% | 0.0% | 0.0% | 9.63% |
| 인프라 | 8.00% | 6.34% | 0.0% | 8.0% | 10.0% |
| PD | 10.00% | 7.51% | 0.0% | 10.0% | 10.0% |

Michaud 평균 포트폴리오는 기대수익률 5.77%, 변동성 10.31%, 공식 mapped 대비 TE 1.06%, 회전율 12.07%이다. 기준 정책 MVO의 TE 1.80%, 회전율 20.75%보다 공식 목표에 가까워졌고, PE가 0%에서 평균 2.0%로 복원되는 등 단일 코너해가 일부 완화되었다. 다만 EM·PD 등은 반복별로 0% 또는 상한에 자주 붙어 분포 폭이 크므로 단일 수치보다 5/50/95% 구간을 함께 해석한다.

## 현재 상태

- STEP 2: 상관/공분산 구축
- STEP 5: 기준 MVO
- STEP 6: Ledoit-Wolf
- STEP 7: Box Robust
- STEP 8: Ellipsoidal Robust
- **STEP 9: Michaud Resampling 완료**
- 다음: 방법론 종합비교 → 2027 팀 목표비중 → Policy Black-Litterman → Stress Test
