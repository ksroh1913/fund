# -*- coding: utf-8 -*-
"""STEP 11 - 2027 Team Target Scenario for Team 2 NPS assignment.

This script does NOT claim to identify a uniquely optimal public-policy target.
It constructs the team's academic scenario by combining:
- Official 2027 common-bucket targets
- Step 10 method synthesis (MVO/LW/Box/Ellipsoid/Michaud)
- 2026.6 actual portfolio implementation gap
- FX-risk and alternatives implementation considerations

Cash is restored explicitly at 0.1%. Steps 5-10 modeled the 99.9% risky sleeve.
For Table C, each model portfolio is scaled by 0.999 and cash 0.001 is appended.
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

RF = 0.03
TOTAL_ASSETS_2026H1_TRN = 1865.6

# Official common-bucket target (2027).
OFFICIAL_COMMON = {
    "국내주식": 0.208,
    "해외주식": 0.356,
    "국내채권": 0.218,
    "해외채권": 0.074,
    "대체투자": 0.143,
    "단기자금": 0.001,
}

# 2026.6 actual NPS portfolio.
ACTUAL_2026H1 = {
    "국내주식": 0.291,
    "해외주식": 0.354,
    "국내채권": 0.154,
    "해외채권": 0.059,
    "대체투자": 0.140,
    "단기자금": 0.002,
}

# Team 2027 scenario.
# This is an IC scenario, not a mechanical copy of any single method.
TEAM_DETAILED = {
    "국내주식": 0.220,
    "글로벌 선진국 DM": 0.275,
    "글로벌 신흥국 EM": 0.065,
    "국내 국채(초장기)": 0.205,
    "글로벌 IG 크레딧": 0.085,
    "사모주식 PE/VC": 0.030,
    "실물 인프라": 0.065,
    "사모대출 PD": 0.054,
    "단기자금": 0.001,
}

# Exact official mapped vector with cash kept explicitly.
OFFICIAL_MAPPED = {
    "국내주식": 0.208,
    "글로벌 선진국 DM": 0.356 * 0.92,
    "글로벌 신흥국 EM": 0.356 * 0.08,
    "국내 국채(초장기)": 0.218,
    "글로벌 IG 크레딧": 0.074,
    "사모주식 PE/VC": 0.143 / 3,
    "실물 인프라": 0.143 / 3,
    "사모대출 PD": 0.143 / 3,
    "단기자금": 0.001,
}

assert abs(sum(TEAM_DETAILED.values()) - 1.0) < 1e-12
assert abs(sum(OFFICIAL_MAPPED.values()) - 1.0) < 1e-12
assert abs(sum(OFFICIAL_COMMON.values()) - 1.0) < 1e-12
assert abs(sum(ACTUAL_2026H1.values()) - 1.0) < 1e-12

proxy = pd.read_csv(DATA / "proxy_mapping.csv")
cov_cma = pd.read_csv(DATA / "covariance_cma.csv", index_col=0)
assets8 = proxy["asset"].tolist()
mu = proxy.set_index("asset").loc[assets8, "mu_cma"].astype(float).values
Sigma = cov_cma.loc[assets8, assets8].astype(float).values

# Step 10 methods are risky-sleeve portfolios summing to 1.0.
table_b = pd.read_csv(RESULTS / "step10_tableB_weights.csv", index_col=0)
table_c = pd.DataFrame(index=assets8 + ["단기자금"])

for method in table_b.columns:
    table_c[method] = 0.0
    table_c.loc[assets8, method] = table_b.loc[assets8, method].values * 0.999
    table_c.loc["단기자금", method] = 0.001

table_c["Team scenario"] = pd.Series(TEAM_DETAILED)
table_c.to_csv(RESULTS / "step11_tableC_detailed.csv", encoding="utf-8-sig")

def reaggregate(s):
    return pd.Series({
        "국내주식": s["국내주식"],
        "해외주식": s["글로벌 선진국 DM"] + s["글로벌 신흥국 EM"],
        "국내채권": s["국내 국채(초장기)"],
        "해외채권": s["글로벌 IG 크레딧"],
        "대체투자": s["사모주식 PE/VC"] + s["실물 인프라"] + s["사모대출 PD"],
        "단기자금": s["단기자금"],
    })

common = pd.DataFrame({c: reaggregate(table_c[c]) for c in table_c.columns})
common.insert(0, "Official 2027", pd.Series(OFFICIAL_COMMON))
common.to_csv(RESULTS / "step11_tableC_common.csv", encoding="utf-8-sig")

# Transition from 2026.6 actual to team scenario.
team_common = reaggregate(pd.Series(TEAM_DETAILED))
transition = pd.DataFrame({
    "2026H1_actual": pd.Series(ACTUAL_2026H1),
    "2027_team": team_common,
})
transition["change_pp"] = (transition["2027_team"] - transition["2026H1_actual"]) * 100
transition["static_balance_equivalent_trn"] = (
    transition["2027_team"] - transition["2026H1_actual"]
) * TOTAL_ASSETS_2026H1_TRN
transition["official_2027"] = pd.Series(OFFICIAL_COMMON)
transition["official_change_pp"] = (transition["official_2027"] - transition["2026H1_actual"]) * 100
transition["official_static_balance_equivalent_trn"] = (
    transition["official_2027"] - transition["2026H1_actual"]
) * TOTAL_ASSETS_2026H1_TRN
transition["incremental_team_vs_official_trn"] = (
    transition["2027_team"] - transition["official_2027"]
) * TOTAL_ASSETS_2026H1_TRN
transition.to_csv(
    RESULTS / "step11_transition_2026H1_to_team.csv",
    encoding="utf-8-sig",
)

# Team scenario model metrics.
w_team = np.array([TEAM_DETAILED[a] for a in assets8])
w_off = np.array([OFFICIAL_MAPPED[a] for a in assets8])
cash = TEAM_DETAILED["단기자금"]

expected_return = float(mu @ w_team + RF * cash)
volatility = float(np.sqrt(w_team @ Sigma @ w_team))
sharpe = float((expected_return - RF) / volatility)
active = w_team - w_off
te = float(np.sqrt(active @ Sigma @ active))
turnover_vs_official = float(
    0.5 * (
        np.abs(active).sum()
        + abs(TEAM_DETAILED["단기자금"] - OFFICIAL_MAPPED["단기자금"])
    )
)

detail_diff = pd.DataFrame({
    "official_mapped": pd.Series(OFFICIAL_MAPPED),
    "team": pd.Series(TEAM_DETAILED),
})
detail_diff["diff_pp"] = (detail_diff["team"] - detail_diff["official_mapped"]) * 100
detail_diff.to_csv(
    RESULTS / "step11_detail_diff_vs_official.csv",
    encoding="utf-8-sig",
)

metadata = {
    "team_common": {k: float(v) for k, v in team_common.items()},
    "team_detailed": TEAM_DETAILED,
    "metrics": {
        "expected_return": expected_return,
        "volatility": volatility,
        "sharpe": sharpe,
        "TE_vs_official_mapped": te,
        "turnover_vs_official_mapped": turnover_vs_official,
    },
    "notes": [
        "Cash 0.1% is fixed outside the 8-asset risky-sleeve optimization.",
        "Static-balance equivalents are not required transaction amounts; market moves, cash flows, and natural drift can change the path.",
        "DM (-5.25%p) and EM (+3.65%p) differ from official mapped weights by more than 2%p. The team deliberately reduces the neutral benchmark's DM concentration while preserving the total foreign-equity allocation close to policy.",
        "The 2026H1-to-2027 implementation table reports both the full transition burden and the incremental burden versus simply implementing the official 2027 target.",
    ],
    "mapping_rule": {
        "foreign_equity": "DM 92% / EM 8%, pre-set neutral mapping independent of Team CMA",
        "alternatives": "PE / Infrastructure / Private Debt equal thirds, pre-set neutral mapping independent of Team CMA",
    },
    "team_judgment": [
        "The Team target is an IC scenario, not a mechanical copy of any Part-2 optimizer.",
        "Domestic equity is kept above the official 20.8% to reduce the 2026H1-to-2027 implementation gap.",
        "Foreign equity remains near the official total while the internal DM/EM split is less concentrated in DM than the neutral mapped benchmark.",
        "Alternatives stay close to the official aggregate because implementation speed and illiquidity constrain rapid changes.",
        "Foreign-asset risk is evaluated on a KRW-unhedged basis in the proxy covariance.",
    ],
    "sources": [
        "Yahoo Finance: original proxy-price source for the frozen 2016-09 to 2026-08 return sample; live download not required by the analysis pipeline.",
        "NPS public allocation disclosure: source for 2026H1 actual weights and total assets used in the implementation-gap illustration.",
        "Official 2027 target weights: provided directly by the assignment, not counted as an additional public-data input."
    ],
}

(RESULTS / "step11_team_target_meta.json").write_text(
    json.dumps(metadata, ensure_ascii=False, indent=2),
    encoding="utf-8",
)

print("Team common target")
print((team_common * 100).round(2))
print("\nMetrics")
print(json.dumps(metadata["metrics"], ensure_ascii=False, indent=2))
