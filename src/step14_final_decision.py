# -*- coding: utf-8 -*-
"""STEP 14 - final IC decision and submission CSV generation.

The calculation pipeline is complete, but the final IC choice is intentionally
left undecided until the common-basis candidate table is reviewed.

Set FINAL_CANDIDATE to one of:
- "Team"
- "Robust BL TE1.00"
- "Blend 50% RBL"
or leave it as None.  None must remain a valid, reproducible pipeline state.
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
FINAL_CANDIDATE=None

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
    resolution=(
      "위원회는 2027년 공식 목표비중의 유지 또는 수정 여부와 자산별 수정폭을 후보 비교 후 확정한다. "
      "이 결정은 기준 CMA, 공분산 추정, CMA 절대수익률 뷰와 데이터 기반 뷰 불확실성에 근거한다. "
      "최종 후보는 Team, Robust BL(TE 1.0%), 50:50 Team/Robust BL 위원회 overlay이며, "
      "사전 TE 1.0%, 최근 8분기 뷰 적중률 55%, 공식 이행경로 대비 추가 국내주식 거래참여율 ADV 1%를 "
      "재심의 기준으로 사용한다. 뷰 적중률 55%는 통계적 유의성 기준이 아니라 팀의 모니터링 기준이다."
    )
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
      "위원회는 2027년 공식 목표비중을 수정한다. 수정폭은 자산별로 최종 비중표와 같다. "
      "이 결정은 기준 CMA, 공분산 추정, CMA 절대수익률 뷰와 데이터 기반 뷰 불확실성에 근거한다. "
      "다만 사전 TE, 뷰 적중률 및 공식 이행경로 대비 추가 거래참여율이 재심의 기준을 벗어나면 "
      "공식 목표 유지 또는 목표비중 재조정을 검토한다."
    )

summary={
 "status":"UNDECIDED" if FINAL_CANDIDATE is None else "DECIDED",
 "decision":FINAL_CANDIDATE,
 "final_common_weights":None if w_final is None else {k:float(v) for k,v in common(w_final).items()},
 "diagnostics":diagnostics,
 "implementation":implementation,
 "assumptions":{"ADV_trn":ADV_TRN,"trading_days":TRADING_DAYS},
 "review_conditions":{
   "ex_ante_TE_max":TE_TRIGGER,
   "view_hit_rate_min":HIT_RATE_TRIGGER,
   "view_hit_rate_definition":"Team governance monitoring rule, not a statistical-significance threshold; realized relative-return direction versus sign(Q-pi), recent 8 quarters.",
   "incremental_domestic_equity_participation_max":PARTICIPATION_TRIGGER,
 },
 "minority_view":"Maintain the official target until sufficient view-hit history is accumulated; use BL only as a validation tool.",
 "resolution_template":resolution,
}
(RESULTS/"step14_final_resolution.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding="utf-8")
md=f"""# STEP 14 IC Decision Status

**Status:** {"UNDECIDED" if FINAL_CANDIDATE is None else FINAL_CANDIDATE}

## Final candidates
- Team
- Robust BL (assignment TE cap 1.0%)
- 50:50 Team / Robust BL committee overlay

TE 0.95% Robust BL and 25/75 blends are sensitivities, not final candidates.

## Resolution template
{resolution}

## Implementation assumptions
- ADV: {ADV_TRN:.0f} trillion KRW
- Trading days: {TRADING_DAYS}
- These are analysis assumptions.
- Full 2026H1-to-target transition burden and incremental burden versus the official target are reported separately once a final candidate is selected.

## Governance note
The 55% view-hit threshold is a Team monitoring rule, not a claim of statistical significance.
"""
(DOCS/"step14_final_resolution.md").write_text(md,encoding="utf-8")
print(resolution)
print(weights)
