# -*- coding: utf-8 -*-
"""STEP 14 - final IC decision and submission CSV generation.

The final decision is deliberately made *after* the common-basis candidate
comparison.  Stress evidence is split:
- Team is much more resilient to direct CMA-view failure.
- Robust BL is more resilient to equity/correlation and alternative-risk shocks.

With no scenario probabilities supplied by the assignment, the committee uses
a transparent 50:50 blend of the independently constructed Team target and the
baseline Robust-BL portfolio.  Among the simple 25/50/75% RBL blends tested in
Step 14 candidate comparison, the 50:50 blend also has the lowest TE and
turnover.  It is therefore a committee overlay, not another optimizer output.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data";RESULTS=ROOT/"results";SUB=ROOT/"submission";DOCS=ROOT/"docs"
for p in (RESULTS,SUB,DOCS):p.mkdir(exist_ok=True)
RF=.03;CASH=.001;AUM=1865.6;ACTUAL_KR=.291
ADV_TRN=20.0      # W04 course-case educational assumption
TRADING_DAYS=240  # W04 course-case educational assumption
PARTICIPATION_TRIGGER=.01
HIT_RATE_TRIGGER=.55
TE_TRIGGER=.01

proxy=pd.read_csv(DATA/"proxy_mapping.csv");assets=proxy["asset"].tolist()
alloc=pd.read_csv(RESULTS/"step12_allocations.csv").set_index("asset")
tc=pd.read_csv(RESULTS/"step11_tableC_detailed.csv",index_col=0)
tableD=pd.read_csv(RESULTS/"step12_tableD.csv")
comparison=pd.read_csv(RESULTS/"step14_candidate_comparison.csv")

w_off=alloc.loc[assets,"official_mapped_total"].astype(float).values
w_blm=alloc.loc[assets,"BL_MVO"].astype(float).values
w_rbl=alloc.loc[assets,"Robust_BL"].astype(float).values
w_team=tc.loc[assets,"Team scenario"].astype(float).values
w_final=.5*w_team+.5*w_rbl

def common(w,cash=CASH):
    s=pd.Series(w,index=assets)
    return pd.Series({
      "국내주식":s["국내주식"],
      "해외주식":s["글로벌 선진국 DM"]+s["글로벌 신흥국 EM"],
      "국내채권":s["국내 국채(초장기)"],
      "해외채권":s["글로벌 IG 크레딧"],
      "대체투자":s[["사모주식 PE/VC","실물 인프라","사모대출 PD"]].sum(),
      "단기자금":cash,
    })
official_common=pd.Series({"국내주식":.208,"해외주식":.356,"국내채권":.218,"해외채권":.074,"대체투자":.143,"단기자금":.001})
team_common=common(w_team);rbl_common=common(w_rbl);final_common=common(w_final)
delta=final_common-official_common

# The exact model diagnostics come from the 50% RBL row generated before this decision.
row50=comparison.loc[comparison["candidate"]=="Blend 50% RBL"].iloc[0]
full_sell=(ACTUAL_KR-final_common["국내주식"])*AUM
incremental_sell=(official_common["국내주식"]-final_common["국내주식"])*AUM
full_participation=full_sell/(ADV_TRN*TRADING_DAYS)
incremental_participation=incremental_sell/(ADV_TRN*TRADING_DAYS)

# Detailed + common submission weights in one auditable file.
detail=pd.DataFrame({
 "level":"detailed","asset":assets,
 "official_or_mapped":w_off,"team_2027":w_team,"BL_MVO":w_blm,
 "Robust_BL":w_rbl,"final_decision":w_final,
})
detail["delta_final_vs_official"]=detail["final_decision"]-detail["official_or_mapped"]
cc=pd.DataFrame({
 "level":"common","asset":official_common.index,
 "official_or_mapped":official_common.values,"team_2027":team_common.values,
 "BL_MVO":common(w_blm).values,"Robust_BL":rbl_common.values,
 "final_decision":final_common.values,
})
cc["delta_final_vs_official"]=cc["final_decision"]-cc["official_or_mapped"]
weights=pd.concat([detail,cc],ignore_index=True)
weights.to_csv(SUB/"team2_weights.csv",index=False,encoding="utf-8-sig")
weights.to_csv(RESULTS/"step14_final_weights.csv",index=False,encoding="utf-8-sig")

views=tableD[["asset","CMA_total","Q_excess","pi_excess","view_error_Q_minus_pi",
              "mu_BL_total","Omega_diag","V_BL_diag"]].copy()
views.to_csv(SUB/"team2_views.csv",index=False,encoding="utf-8-sig")

# Korean final resolution: target, delta, rationale, and three numeric review triggers.
resolution=(
"위원회는 2027년 공식 목표비중을 수정하여 국내주식 20.6%, 해외주식 33.4%, 국내채권 22.7%, "
"해외채권 8.9%, 대체투자 14.4%, 단기자금 0.1%로 설정한다. 공식 대비 수정폭은 각각 -0.2%p, "
"-2.2%p, +0.9%p, +1.5%p, +0.1%p, 0.0%p다. Team안은 CMA 뷰 실패에 강하고 Robust BL은 "
"주식·대체 충격에 강해 두 안을 50:50 결합했으며, 단순 25·50·75% 혼합안 중 50:50이 TE와 "
"회전율도 가장 낮았다. 사전 TE가 1.0%를 초과하거나 최근 8분기 뷰 적중률이 55% 미만이거나 "
"공식 이행경로 대비 추가 국내주식 거래참여율이 ADV 1%를 초과하면 공식 목표로 복귀하거나 재심의한다."
)

summary={
 "decision":"50:50 Team / baseline Robust BL committee blend",
 "final_common_weights":{k:float(v) for k,v in final_common.items()},
 "final_common_delta_vs_official":{k:float(v) for k,v in delta.items()},
 "diagnostics":{
   "expected_return_CMA":float(row50["expected_return"]),
   "volatility_common_Sigma":float(row50["volatility_common_Sigma"]),
   "utility":float(row50["utility"]),"TE":float(row50["TE"]),"turnover":float(row50["turnover"]),
   "stress1_active_loss":float(row50["stress1_active_loss"]),
   "stress2_utility_change":float(row50["stress2_utility_change"]),
   "stress3_Higham_utility_change":float(row50["stress3_Higham_utility_change"]),
 },
 "implementation":{
   "2026H1_to_final_domestic_equity_sell_trn":float(full_sell),
   "incremental_vs_official_domestic_equity_sell_trn":float(incremental_sell),
   "W04_case_ADV_trn":ADV_TRN,"W04_case_trading_days":TRADING_DAYS,
   "full_path_participation_rate":float(full_participation),
   "incremental_vs_official_participation_rate":float(incremental_participation),
   "interpretation":"The full transition burden is reported for feasibility; the review trigger focuses on incremental burden caused by modifying the official target."
 },
 "review_conditions":{
   "ex_ante_TE_max":TE_TRIGGER,
   "view_hit_rate_min":HIT_RATE_TRIGGER,
   "view_hit_rate_definition":"Team governance monitoring rule, not a statistical-significance threshold; monitor whether realized relative-return direction agrees with the sign of Q-pi over the recent 8 quarters.",
   "incremental_domestic_equity_participation_max":PARTICIPATION_TRIGGER,
 },
 "minority_view":"Maintain the official target until sufficient view-hit history is accumulated; use BL only as a validation tool.",
 "resolution":resolution,
}
(RESULTS/"step14_final_resolution.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
md=f"""# STEP 14 Final IC Resolution

## Final common target
{pd.DataFrame({"Official":official_common,"Final":final_common,"Delta":delta}).to_markdown()}

## Decision
{resolution}

## Why 50:50
- Team materially limits the direct-CMA view-miss active loss.
- Robust BL materially reduces equity/correlation and alternative-risk shock utility losses.
- With no assignment-provided scenario probabilities, an equal-weight committee overlay is transparent.
- In the pre-specified 25/50/75 blend sensitivity, 50:50 has the lowest TE and turnover.

## Implementation
- Full 2026H1 domestic-equity reduction to final target: {full_sell:.1f} trillion KRW static-balance equivalent.
- Incremental reduction versus implementing the official 20.8% target: {incremental_sell:.1f} trillion KRW.
- Using the W04 educational assumptions ADV={ADV_TRN:.0f} trillion KRW and {TRADING_DAYS} trading days, incremental participation is {incremental_participation*100:.2f}% of ADV.
- The full transition burden is shown separately because 2026 actual-to-2027 implementation feasibility remains a Part-3 consideration.

## Governance note
The 55% view-hit threshold is a Team monitoring rule, not a claim of statistical significance.
"""
(DOCS/"step14_final_resolution.md").write_text(md,encoding="utf-8")
print(resolution)
print(pd.DataFrame({"Official":official_common,"Final":final_common,"Delta":delta}))
