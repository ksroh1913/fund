# -*- coding: utf-8 -*-
"""STEP 10 - Method synthesis for Team 2 NPS assignment.

Outputs:
- Table B: Official mapped, policy MVO, Ledoit-Wolf, Box, Ellipsoid, Michaud
- 50bp expected-return sensitivity
- Small-eigenvalue/eigenvector analysis
- Reaggregated common buckets
"""

from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.covariance import LedoitWolf

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

GAMMA = 4.0
RF = 0.03

OFFICIAL_TOP = {
    "국내주식": 0.208,
    "해외주식": 0.356,
    "국내채권": 0.218,
    "해외채권": 0.074,
    "대체투자": 0.143,
}
DM_SHARE = 0.92
EM_SHARE = 0.08
ALT_SPLIT = {
    "사모주식 PE/VC": 1 / 3,
    "실물 인프라": 1 / 3,
    "사모대출 PD": 1 / 3,
}
DETAIL_BOUNDS = {
    "글로벌 신흥국 EM": (0.0, 0.10),
    "사모주식 PE/VC": (0.0, 0.10),
    "실물 인프라": (0.0, 0.10),
    "사모대출 PD": (0.0, 0.10),
}

proxy = pd.read_csv(DATA / "proxy_mapping.csv")
rets = pd.read_csv(DATA / "monthly_returns_krw.csv", index_col=0)
cov_cma = pd.read_csv(DATA / "covariance_cma.csv", index_col=0)

assets = proxy["asset"].tolist()
idx = {a: i for i, a in enumerate(assets)}
mu = proxy.set_index("asset").loc[assets, "mu_cma"].astype(float).values
sigma = proxy.set_index("asset").loc[assets, "sigma_cma"].astype(float).values
Sigma = cov_cma.loc[assets, assets].astype(float).values
R = rets[assets].astype(float).values
n = len(assets)

scale = 1.0 / sum(OFFICIAL_TOP.values())
w_off = pd.Series(0.0, index=assets)
w_off["국내주식"] = OFFICIAL_TOP["국내주식"] * scale
w_off["글로벌 선진국 DM"] = OFFICIAL_TOP["해외주식"] * DM_SHARE * scale
w_off["글로벌 신흥국 EM"] = OFFICIAL_TOP["해외주식"] * EM_SHARE * scale
w_off["국내 국채(초장기)"] = OFFICIAL_TOP["국내채권"] * scale
w_off["글로벌 IG 크레딧"] = OFFICIAL_TOP["해외채권"] * scale
for asset, share in ALT_SPLIT.items():
    w_off[asset] = OFFICIAL_TOP["대체투자"] * share * scale

def group_sum(names):
    ids = [idx[a] for a in names]
    return lambda w: w[ids].sum()

policy_cons = [
    {"type":"eq", "fun":lambda w:w.sum()-1.0},
    {"type":"ineq", "fun":lambda w:w[idx["국내주식"]]-0.15},
    {"type":"ineq", "fun":lambda w:0.25-w[idx["국내주식"]]},
    {"type":"ineq", "fun":lambda w:group_sum(["글로벌 선진국 DM","글로벌 신흥국 EM"])(w)-0.30},
    {"type":"ineq", "fun":lambda w:0.42-group_sum(["글로벌 선진국 DM","글로벌 신흥국 EM"])(w)},
    {"type":"ineq", "fun":lambda w:w[idx["국내 국채(초장기)"]]-0.15},
    {"type":"ineq", "fun":lambda w:0.25-w[idx["국내 국채(초장기)"]]},
    {"type":"ineq", "fun":lambda w:w[idx["글로벌 IG 크레딧"]]-0.05},
    {"type":"ineq", "fun":lambda w:0.12-w[idx["글로벌 IG 크레딧"]]},
    {"type":"ineq", "fun":lambda w:group_sum(["사모주식 PE/VC","실물 인프라","사모대출 PD"])(w)-0.12},
    {"type":"ineq", "fun":lambda w:0.18-group_sum(["사모주식 PE/VC","실물 인프라","사모대출 PD"])(w)},
    {"type":"ineq", "fun":lambda w:group_sum(["국내주식","글로벌 선진국 DM","글로벌 신흥국 EM"])(w)-0.52},
    {"type":"ineq", "fun":lambda w:0.60-group_sum(["국내주식","글로벌 선진국 DM","글로벌 신흥국 EM"])(w)},
    {"type":"ineq", "fun":lambda w:group_sum(["국내 국채(초장기)","글로벌 IG 크레딧"])(w)-0.25},
    {"type":"ineq", "fun":lambda w:0.35-group_sum(["국내 국채(초장기)","글로벌 IG 크레딧"])(w)},
]

bounds = [(0.0, 1.0)] * n
for asset, bound in DETAIL_BOUNDS.items():
    bounds[idx[asset]] = bound

def solve(mu_vec, cov, x0=None):
    if x0 is None:
        x0 = w_off.values.copy()
    res = minimize(
        lambda w: -(mu_vec @ w - GAMMA / 2 * (w @ cov @ w)),
        x0,
        method="SLSQP",
        bounds=bounds,
        constraints=policy_cons,
        options={"ftol":1e-12, "maxiter":5000},
    )
    if not res.success:
        raise RuntimeError(res.message)
    return res.x

def metrics(w, risk_cov=Sigma):
    er = float(mu @ w)
    vol = float(np.sqrt(w @ risk_cov @ w))
    active = w - w_off.values
    return {
        "expected_return": er,
        "volatility": vol,
        "sharpe": float((er - RF) / vol),
        "TE": float(np.sqrt(active @ Sigma @ active)),
        "turnover": float(0.5 * np.abs(active).sum()),
    }

# Recompute point-estimate methods.
w_mvo = solve(mu, Sigma)

lw = LedoitWolf().fit(R)
S_lw = lw.covariance_ * 12
sd_lw = np.sqrt(np.diag(S_lw))
rho_lw = S_lw / np.outer(sd_lw, sd_lw)
Sigma_lw = np.outer(sigma, sigma) * rho_lw
w_lw = solve(mu, Sigma_lw, w_mvo)

box = pd.read_csv(RESULTS / "step7_box.csv")
w_box = box.set_index("자산").loc[assets, "Box비중"].values

ell = pd.read_csv(RESULTS / "step8_ellipsoid.csv")
ell_row = ell.iloc[(ell["p"] - 0.5).abs().argsort()[:1]]
w_ell = ell_row[assets].iloc[0].values

mich = pd.read_csv(RESULTS / "step9_michaud_summary.csv")
w_mich = mich.set_index("asset").loc[assets, "michaud_mean"].values

methods = {
    "Official mapped": w_off.values,
    "Policy MVO": w_mvo,
    "Ledoit-Wolf": w_lw,
    "Box Robust": w_box,
    "Ellipsoid Robust": w_ell,
    "Michaud mean": w_mich,
}

table_b = pd.DataFrame(methods, index=assets)
table_b.to_csv(RESULTS / "step10_tableB_weights.csv", encoding="utf-8-sig")

metric_rows = []
for name, w in methods.items():
    metric_rows.append({
        "method": name,
        **metrics(w, Sigma_lw if name == "Ledoit-Wolf" else Sigma),
    })
pd.DataFrame(metric_rows).to_csv(
    RESULTS / "step10_tableB_metrics.csv", index=False, encoding="utf-8-sig"
)

reagg = []
for name, w in methods.items():
    s = pd.Series(w, index=assets)
    reagg.append({
        "method": name,
        "국내주식": s["국내주식"],
        "해외주식": s["글로벌 선진국 DM"] + s["글로벌 신흥국 EM"],
        "국내채권": s["국내 국채(초장기)"],
        "해외채권": s["글로벌 IG 크레딧"],
        "대체투자": s[["사모주식 PE/VC","실물 인프라","사모대출 PD"]].sum(),
    })
pd.DataFrame(reagg).to_csv(
    RESULTS / "step10_reaggregated.csv", index=False, encoding="utf-8-sig"
)

# +/-50bp expected-return sensitivity, one asset at a time.
sens = []
for j, asset in enumerate(assets):
    for shock in (-0.005, 0.005):
        mu2 = mu.copy()
        mu2[j] += shock
        w = solve(mu2, Sigma, w_mvo)
        sens.append({
            "asset_shocked": asset,
            "shock": shock,
            "own_weight_base": w_mvo[j],
            "own_weight_new": w[j],
            "own_weight_change": w[j] - w_mvo[j],
            "portfolio_L1_change": float(np.abs(w - w_mvo).sum()),
            "portfolio_turnover_change": float(0.5 * np.abs(w - w_mvo).sum()),
        })

sens_df = pd.DataFrame(sens)
sens_df.to_csv(
    RESULTS / "step10_mu_50bp_sensitivity.csv", index=False, encoding="utf-8-sig"
)

summary_rows = []
for asset, g in sens_df.groupby("asset_shocked", sort=False):
    own_idx = g["own_weight_change"].abs().idxmax()
    l1_idx = g["portfolio_L1_change"].idxmax()
    own = g.loc[own_idx]
    l1 = g.loc[l1_idx]
    summary_rows.append({
        "asset": asset,
        "max_abs_own_weight_change": abs(own["own_weight_change"]),
        "shock_for_max_own": own["shock"],
        "max_portfolio_L1_change": l1["portfolio_L1_change"],
        "shock_for_max_L1": l1["shock"],
    })
pd.DataFrame(summary_rows).sort_values(
    "max_abs_own_weight_change", ascending=False
).to_csv(
    RESULTS / "step10_mu_50bp_sensitivity_summary.csv",
    index=False,
    encoding="utf-8-sig",
)

# Small-eigenvalue analysis.
evals, evecs = np.linalg.eigh(Sigma)
order = np.argsort(evals)
evals = evals[order]
evecs = evecs[:, order]

rows = []
for k in range(min(4, n)):
    v = evecs[:, k]
    v = v / np.max(np.abs(v))
    rows.append({
        "rank_smallest": k + 1,
        "eigenvalue": evals[k],
        **{assets[i]: v[i] for i in range(n)},
    })
pd.DataFrame(rows).to_csv(
    RESULTS / "step10_small_eigenvectors.csv", index=False, encoding="utf-8-sig"
)

v = evecs[:, 0] / np.max(np.abs(evecs[:, 0]))
meta = {
    "baseline_gamma": GAMMA,
    "lw_shrinkage": float(lw.shrinkage_),
    "smallest_eigenvalue": float(evals[0]),
    "smallest_eigenvector": {assets[i]: float(v[i]) for i in range(n)},
    "primary_reference_method":
        "Michaud mean (per-draw Ledoit-Wolf + resampling)",
    "reason":
        "Smooths single-sample corner solutions, incorporates joint mu/Sigma "
        "sampling uncertainty, and produced lower TE/turnover versus point-estimate MVO; "
        "Ellipsoid is retained as a robustness check.",
}
(RESULTS / "step10_meta.json").write_text(
    json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8"
)

print((table_b * 100).round(2))
print(pd.DataFrame(metric_rows))
