# fund

국민연금 CMA 기반 2027 목표비중 및 MVO/Robust/Black-Litterman 종합과제 작업 저장소입니다.

현재 커밋에는 **STEP 2 및 STEP 5~8 재현 코드와 중간 산출물**을 정리했습니다.

## 구조

```
src/
  step2_build_corr.py       # 8자산 실제 프록시 월간수익률 → 상관/공분산 구축
  analyze_steps5_8.py       # 기준 MVO, Ledoit-Wolf, Box, Ellipsoid 재현

results/
  step5_mvo.csv
  step6_lw_diagnostics.csv
  step7_box.csv
  step8_S_mu.csv
  step8_ellipsoid.csv

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

## 현재 상태

- STEP 2: 상관/공분산 구축
- STEP 5: 기준 MVO
- STEP 6: Ledoit-Wolf
- STEP 7: Box Robust
- STEP 8: Ellipsoidal Robust

이후 STEP 9 Michaud Resampling, 2027 팀 목표비중, Policy Black-Litterman, Stress Test를 같은 구조에 추가할 예정입니다.
