# STEP 14 Final IC Resolution

## Final common target

| 자산군 | Official 2027 | Final | Delta |
|---|---:|---:|---:|
| 국내주식 | 20.80% | 20.58% | -0.22%p |
| 해외주식 | 35.60% | 33.42% | -2.18%p |
| 국내채권 | 21.80% | 22.65% | +0.85%p |
| 해외채권 | 7.40% | 8.87% | +1.47%p |
| 대체투자 | 14.30% | 14.38% | +0.08%p |
| 단기자금 | 0.10% | 0.10% | 0.00%p |

## Decision

위원회는 2027년 공식 목표비중을 수정하여 국내주식 20.6%, 해외주식 33.4%, 국내채권 22.7%, 해외채권 8.9%, 대체투자 14.4%, 단기자금 0.1%로 설정한다. 공식 대비 수정폭은 각각 -0.2%p, -2.2%p, +0.9%p, +1.5%p, +0.1%p, 0.0%p다. Team안은 CMA 뷰 실패에 강하고 Robust BL은 주식·대체 충격에 강해 두 안을 50:50 결합했으며, 단순 25·50·75% 혼합안 중 50:50이 TE와 회전율도 가장 낮았다. 사전 TE가 1.0%를 초과하거나 최근 8분기 뷰 적중률이 55% 미만이거나 공식 이행경로 대비 추가 국내주식 거래참여율이 ADV 1%를 초과하면 공식 목표로 복귀하거나 재심의한다.

## Why 50:50

- Team direct-CMA Stress 1 active loss: 약 -3.5bp.
- Robust BL direct-CMA Stress 1 active loss: 약 -66.6bp.
- Robust BL is better under equity/correlation and alternative-risk stresses.
- 50:50 blend: CMA ER 5.61%, common-Sigma volatility 10.20%, TE 0.73%, turnover 6.56%.
- 50:50 blend Stress 1 active loss: 약 -35.0bp; Stress 2 utility change: -1.24%p; Stress 3 Higham utility change: -0.32%p.
- No scenario probabilities are supplied, so the equal blend is a transparent committee overlay rather than a claim of statistical optimality.

## Implementation

- 2026H1 domestic equity 29.1% → final 20.58%: static-balance equivalent sell amount 약 158.9조원.
- Official 20.8% path 대비 additional sell amount: 약 4.0조원.
- W04 educational assumptions (ADV 20조원, 240 trading days) 기준 additional participation: 약 0.08%.
- Full transition burden is still disclosed separately because 2026 actual-to-2027 feasibility is required in Part 3.

## Governance

- TE trigger: ex-ante TE > 1.0%.
- View trigger: recent 8-quarter hit rate < 55%; this is a Team governance monitoring rule, not a statistical-significance threshold.
- Execution trigger: incremental domestic-equity participation caused by modifying the official path > 1% of ADV.
- Minority view: maintain the official target until enough view-hit history is accumulated and use BL only as a validation tool.
