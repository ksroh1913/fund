# -*- coding: utf-8 -*-
"""STEP 9 - Michaud resampling for Team 2 NPS 8-asset assignment.

Method
------
- B = 300, seed = 60080.
- Resample all T=120 monthly KRW proxy returns with replacement.
- Recenter each bootstrap mean shock around the Team 2 forward-looking CMA mean:
      mu_b = mu_CMA + 12 * (mean(R_b) - mean(R_full))
  Historical data therefore supplies estimation-error shape, while CMA remains
  the central expected-return forecast.
- In every draw, estimate a Ledoit-Wolf covariance from the resampled returns,
  convert it to a correlation matrix, and re-scale it by Team 2 CMA volatilities.
- Optimize the same gamma=4 policy-constrained MVO used in Steps 5-8.
- Report mean and 5/50/95 percentiles of weights.
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

SEED = 60080
B = 300
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

# Team analysis assumptions; not official NPS detailed limits.
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
mu = proxy.set_index("asset").loc[assets, "mu_cma"].astype(float).values
sigma = proxy.set_index("asset").loc[assets, "sigma_cma"].astype(float).values
Sigma = cov_cma.loc[assets, assets].astype(float).values
R = rets[assets].astype(float).values
T, n = R.shape
idx = {a: i for i, a in enumerate(assets)}

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
    {"type": "eq", "fun": lambda w: w.sum() - 1.0},
    {"type": "ineq", "fun": lambda w: w[idx["국내주식"]] - 0.15},
    {"type": "ineq", "fun": lambda w: 0.25 - w[idx["국내주식"]]},
    {"type": "ineq", "fun": lambda w: group_sum(["글로벌 선진국 DM", "글로벌 신흥국 EM"])(w) - 0.30},
    {"type": "ineq", "fun": lambda w: 0.42 - group_sum(["글로벌 선진국 DM", "글로벌 신흥국 EM"])(w)},
    {"type": "ineq", "fun": lambda w: w[idx["국내 국채(초장기)"]] - 0.15},
    {"type": "ineq", "fun": lambda w: 0.25 - w[idx["국내 국채(초장기)"]]},
    {"type": "ineq", "fun": lambda w: w[idx["글로벌 IG 크레딧"]] - 0.05},
    {"type": "ineq", "fun": lambda w: 0.12 - w[idx["글로벌 IG 크레딧"]]},
    {"type": "ineq", "fun": lambda w: group_sum(["사모주식 PE/VC", "실물 인프라", "사모대출 PD"])(w) - 0.12},
    {"type": "ineq", "fun": lambda w: 0.18 - group_sum(["사모주식 PE/VC", "실물 인프라", "사모대출 PD"])(w)},
    {"type": "ineq", "fun": lambda w: group_sum(["국내주식", "글로벌 선진국 DM", "글로벌 신흥국 EM"])(w) - 0.52},
    {"type": "ineq", "fun": lambda w: 0.60 - group_sum(["국내주식", "글로벌 선진국 DM", "글로벌 신흥국 EM"])(w)},
    {"type": "ineq", "fun": lambda w: group_sum(["국내 국채(초장기)", "글로벌 IG 크레딧"])(w) - 0.25},
    {"type": "ineq", "fun": lambda w: 0.35 - group_sum(["국내 국채(초장기)", "글로벌 IG 크레딧"])(w)},
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
        options={"ftol": 1e-12, "maxiter": 5000},
    )
    if not res.success:
        raise RuntimeError(res.message)
    return res.x

def metrics(w):
    expected_return = float(mu @ w)
    volatility = float(np.sqrt(w @ Sigma @ w))
    active = w - w_off.values
    return {
        "expected_return": expected_return,
        "volatility": volatility,
        "sharpe": float((expected_return - RF) / volatility),
        "TE": float(np.sqrt(active @ Sigma @ active)),
        "turnover": float(0.5 * np.abs(active).sum()),
    }

w_base = solve(mu, Sigma)

rng = np.random.default_rng(SEED)
full_mean = R.mean(axis=0)
W = np.empty((B, n))
mu_draws = np.empty((B, n))
shrinkage = np.empty(B)
previous = w_base.copy()

for b in range(B):
    bootstrap_index = rng.integers(0, T, T)
    R_b = R[bootstrap_index]

    mu_b = mu + 12.0 * (R_b.mean(axis=0) - full_mean)

    lw = LedoitWolf().fit(R_b)
    S_lw = lw.covariance_ * 12.0
    sd_lw = np.sqrt(np.diag(S_lw))
    rho_lw = S_lw / np.outer(sd_lw, sd_lw)
    Sigma_b = np.outer(sigma, sigma) * rho_lw

    try:
        w_b = solve(mu_b, Sigma_b, previous)
    except RuntimeError:
        w_b = solve(mu_b, Sigma_b, w_off.values.copy())

    W[b] = w_b
    mu_draws[b] = mu_b
    shrinkage[b] = lw.shrinkage_
    previous = w_b

Wdf = pd.DataFrame(W, columns=assets)
summary = pd.DataFrame(
    {
        "baseline_mvo": w_base,
        "michaud_mean": Wdf.mean().values,
        "p05": Wdf.quantile(0.05).values,
        "p50": Wdf.quantile(0.50).values,
        "p95": Wdf.quantile(0.95).values,
        "std": Wdf.std(ddof=1).values,
    },
    index=assets,
)
summary["change_vs_baseline"] = summary["michaud_mean"] - summary["baseline_mvo"]

upper_map = {a: DETAIL_BOUNDS.get(a, (0.0, 1.0))[1] for a in assets}
summary["freq_at_zero"] = [float((Wdf[a] <= 1e-6).mean()) for a in assets]
summary["freq_at_upper"] = [
    float((Wdf[a] >= upper_map[a] - 1e-6).mean()) if upper_map[a] < 1.0 else np.nan
    for a in assets
]

w_mean = summary["michaud_mean"].values
baseline_metrics = metrics(w_base)
michaud_metrics = metrics(w_mean)
mean_groups = {
    "국내주식": w_mean[idx["국내주식"]],
    "해외주식": w_mean[idx["글로벌 선진국 DM"]] + w_mean[idx["글로벌 신흥국 EM"]],
    "국내채권": w_mean[idx["국내 국채(초장기)"]],
    "해외채권": w_mean[idx["글로벌 IG 크레딧"]],
    "대체투자": sum(w_mean[idx[a]] for a in ["사모주식 PE/VC", "실물 인프라", "사모대출 PD"]),
    "총주식": sum(w_mean[idx[a]] for a in ["국내주식", "글로벌 선진국 DM", "글로벌 신흥국 EM"]),
    "총채권": sum(w_mean[idx[a]] for a in ["국내 국채(초장기)", "글로벌 IG 크레딧"]),
}

Wdf.index.name = "iteration"
Wdf.to_csv(RESULTS / "step9_michaud_draws.csv", encoding="utf-8-sig")
summary.reset_index(names="asset").to_csv(
    RESULTS / "step9_michaud_summary.csv", index=False, encoding="utf-8-sig"
)
pd.DataFrame(mu_draws, columns=assets).to_csv(
    RESULTS / "step9_mu_draws.csv", index=False, encoding="utf-8-sig"
)
pd.DataFrame(
    {"iteration": np.arange(B), "lw_shrinkage": shrinkage}
).to_csv(RESULTS / "step9_lw_shrinkage.csv", index=False, encoding="utf-8-sig")

metadata = {
    "B": B,
    "seed": SEED,
    "T": T,
    "gamma": GAMMA,
    "method": "monthly bootstrap; recentered CMA mean shock; per-draw Ledoit-Wolf correlation; CMA marginal sigma; same policy constraints",
    "baseline_metrics": baseline_metrics,
    "michaud_mean_metrics": michaud_metrics,
    "average_lw_shrinkage": float(shrinkage.mean()),
    "lw_shrinkage_p05": float(np.quantile(shrinkage, 0.05)),
    "lw_shrinkage_p95": float(np.quantile(shrinkage, 0.95)),
    "mean_group_weights": {k: float(v) for k, v in mean_groups.items()},
    "sum_mean_weights": float(w_mean.sum()),
}
(RESULTS / "step9_michaud_meta.json").write_text(
    json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
)

print((summary * 100).round(2))
print(json.dumps(metadata, ensure_ascii=False, indent=2))
