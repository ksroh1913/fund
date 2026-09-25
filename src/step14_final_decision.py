# -*- coding: utf-8 -*-
"""STEP 14 - final IC decision and submission CSV generation.

The final IC choice is Robust BL with the assignment TE cap of 1.0%.

The decision follows the assignment stress table and a supplemental fairness
diagnostic that shocks each candidate against its own active direction.
Team and 50:50 remain comparison alternatives; TE 0.95% remains a sensitivity.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data";RESULTS=ROOT/"results";SUB=ROOT/"submission";DOCS=ROOT/"docs"
for p in (RESULTS,SUB,DOCS):p.mkdir(exist_ok=True)
RF=.03;CASH=.001;AUM=1865.6;ACTUAL_KR=.291
ADV_TRN=20.0
TRADING_DAYS=240
PARTICIPATION_TRIGGER=.01
HIT_RATE_TRIGGER=.55
TE_TRIGGER=.01
FINAL_CANDIDATE="Robust BL TE1.00"

proxy=pd.read_csv(DATA/"proxy_mapping.csv");assets=proxy["asset"].tolist()
alloc=pd.read_csv(RESULTS/"step12_allocations.csv").set_index("asset")
tc=pd.read_csv(RESULTS/"step11_tableC_detailed.csv",index_col=0)
tableD=pd.read_csv(RESULTS/"step12_tableD.csv")
comparison=pd.read_csv(RESULTS/"step14_candidate_comparison.csv").set_index("candidate")

w_off=alloc.loc[assets,"official_mapped_total"].astype(float).values
w_blm=alloc.loc[assets,"BL_MVO"].astype(float).values
w_rbl=alloc.loc[assets,"Robust_BL"].astype(float).values
w_team=tc.loc[assets,"Team scenario"].astype(float).values
candidate_weights={
 "Team":w_team,
 "Robust BL TE1.00":w_rbl,
 "Blend 50% RBL":.5*w_team+.5*w_rbl,
}
if FINAL_CANDIDATE not in ({None}|set(candidate_weights)):
    raise ValueError(f"Unknown FINAL_CANDIDATE: {FINAL_CANDIDATE}")
w_final=None if FINAL_CANDIDATE is None else candidate_weights[FINAL_CANDIDATE]

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
team_common=common(w_team);rbl_common=common(w_rbl)

# Candidate columns are always populated; final_decision remains blank while undecided.
detail=pd.DataFrame({
 "level":"detailed","asset":assets,
 "official_or_mapped":w_off,"team_2027":w_team,"BL_MVO":w_blm,
 "Robust_BL":w_rbl,"Blend_50_RBL":candidate_weights["Blend 50% RBL"],
})
cc=pd.DataFrame({
 "level":"common","asset":official_common.index,
 "official_or_mapped":official_common.values,"team_2027":team_common.values,
 "BL_MVO":common(w_blm).values,"Robust_BL":rbl_common.values,
 "Blend_50_RBL":common(candidate_weights["Blend 50% RBL"]).values,
})
if w_final is None:
    detail["final_decision"]=np.nan;cc["final_decision"]=np.nan
    detail["delta_final_vs_official"]=np.nan;cc["delta_final_vs_official"]=np.nan
else:
    final_common=common(w_final)
    detail["final_decision"]=w_final
    detail["delta_final_vs_official"]=detail["final_decision"]-detail["official_or_mapped"]
    cc["final_decision"]=final_common.values
    cc["delta_final_vs_official"]=cc["final_decision"]-cc["official_or_mapped"]
weights=pd.concat([detail,cc],ignore_index=True)
weights.to_csv(SUB/"team2_weights.csv",index=False,encoding="utf-8-sig")
weights.to_csv(RESULTS/"step14_final_weights.csv",index=False,encoding="utf-8-sig")

views=tableD[["asset","CMA_total","Q_excess","pi_excess","view_error_Q_minus_pi",
              "mu_BL_total","Omega_diag","V_BL_diag"]].copy()
views.to_csv(SUB/"team2_views.csv",index=False,encoding="utf-8-sig")

if FINAL_CANDIDATE is None:
    resolution=""
    diagnostics={}
    implementation={}
else:
    final_common=common(w_final);delta=final_common-official_common
    row=comparison.loc[FINAL_CANDIDATE]
    full_sell=(ACTUAL_KR-final_common["국내주식"])*AUM
    incremental_sell=(official_common["국내주식"]-final_common["국내주식"])*AUM
    diagnostics={k:float(row[k]) for k in [
      "expected_return","expected_return_vs_official","volatility_common_Sigma",
      "utility","TE","turnover","stress1_active_loss",
      "stress2_utility_change","stress3_Higham_utility_change"]}
    implementation={
      "2026H1_to_final_domestic_equity_sell_trn":float(full_sell),
      "incremental_vs_official_domestic_equity_sell_trn":float(incremental_sell),
      "assumed_ADV_trn":ADV_TRN,"assumed_trading_days":TRADING_DAYS,
      "full_path_participation_rate":float(full_sell/(ADV_TRN*TRADING_DAYS)),
      "incremental_vs_official_participation_rate":float(incremental_sell/(ADV_TRN*TRADING_DAYS)),
    }
    resolution=(
      "위원회는 2027년 목표비중을 수정한다. 수정폭은 자산별로 국내주식 -1.6%p, 해외주식 -2.8%p, "
      "국내채권 +3.0%p, 해외채권 +1.8%p, 대체투자 -0.4%p, 단기자금 0.0%p다. "
      "이 결정은 기준 CMA, 공분산 추정, CMA 절대수익률 뷰와 데이터 기반 뷰 불확실성에 근거한다. "
      "공식 비중을 사전비중으로 둔 Robust BL은 CMA 기준 기대수익률이 공식보다 약 5bp 낮지만 변동성을 0.9%p 낮춘다. "
      "다만 운용 중 실제 비중 드리프트로 사전 TE가 1.0%를 넘으면 리밸런싱하고, 리밸런싱 후에도 초과하거나 최근 "
      "8분기 뷰 적중률이 55% 미만이거나 공식 경로 대비 추가 국내주식 거래가 ADV의 1%를 넘으면 공식 목표로 복귀하거나 재심의한다."
    )

summary={
 "status":"UNDECIDED" if FINAL_CANDIDATE is None else "DECIDED",
 "decision":FINAL_CANDIDATE,
 "final_common_weights":None if w_final is None else {k:float(v) for k,v in common(w_final).items()},
 "diagnostics":diagnostics,
 "binding_constraints":{
   "국내채권_active_pp":3.0,
   "사모주식_PE_active_pp":-3.0,
   "사모대출_PD_active_pp":3.0,
   "target_weight_ex_ante_TE_pct":1.0
 },
 "implementation":implementation,
 "assumptions":{"ADV_trn":ADV_TRN,"trading_days":TRADING_DAYS},
 "review_conditions":{
   "target_weight_ex_ante_TE":1.0e-2,
   "ex_ante_TE_limit":TE_TRIGGER,
   "TE_monitoring_rule":"Target-weight ex-ante TE is 1.00% at the limit. During implementation, if drifted actual-weight ex-ante TE exceeds 1.0%, rebalance; if it remains above 1.0% after rebalancing, reconsider the decision.",
   "view_hit_rate_min":HIT_RATE_TRIGGER,
   "view_hit_rate_definition":"Team governance monitoring rule, not a statistical-significance threshold; realized relative-return direction versus sign(Q-pi), recent 8 quarters.",
   "incremental_domestic_equity_participation_max":PARTICIPATION_TRIGGER
 },
 "minority_view":"Team안 채택. 팀 CMA 기준 기대수익률이 공식보다 약 15bp 높고, 공식 경로 대비 국내주식 매도 부담이 약 22조원 적으며, 6개 공통분류 기준 수정폭이 모두 ±3%p 이내라는 점을 근거로 한다.",
 "resolution_template":resolution,
}
(RESULTS/"step14_final_resolution.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
md=f"""# STEP 14 IC Decision Status

**Status:** {"UNDECIDED" if FINAL_CANDIDATE is None else FINAL_CANDIDATE}

## Final decision
- Robust BL (assignment TE cap 1.0%)
- Final common weights are generated directly from the Step 12 Robust BL result.
- Team and 50:50 remain comparison alternatives; TE 0.95% Robust BL is a sensitivity.

## Resolution
{resolution}

## Implementation assumptions
- ADV: {ADV_TRN:.0f} trillion KRW
- Trading days: {TRADING_DAYS}
- These are analysis assumptions.
- Full 2026H1-to-target transition burden and incremental burden versus the official target are reported separately once a final candidate is selected.

## Governance note
- Target-weight ex-ante TE is 1.00%, inside the assignment limit but exactly on the boundary.
- During implementation, if drifted actual-weight ex-ante TE exceeds 1.0%, rebalance; if it still exceeds 1.0%, reconsider.
- The 55% view-hit threshold is a Team monitoring rule, not a claim of statistical significance.
- With ADV 20 trillion KRW and 240 trading days assumed, the additional domestic-equity participation versus the official path is about 0.64%, below the 1% review threshold.

## Minority view
Team안을 채택해야 한다는 소수의견을 유지한다. 팀 CMA 기준 기대수익률이 공식보다 약 15bp 높고, 공식 경로 대비 국내주식 매도 부담이 약 22조원 적으며, 6개 공통분류 기준 수정폭이 모두 ±3%p 이내라는 점이 근거다.
"""
(DOCS/"step14_final_resolution.md").write_text(md,encoding="utf-8")
print(resolution)
print(weights)
