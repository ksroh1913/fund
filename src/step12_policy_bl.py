# -*- coding: utf-8 -*-
"""STEP 12 - Policy Black-Litterman for Team 2 NPS assignment.

Key assignment choices
----------------------
- Prior portfolio: 2027 official mapped weights, NOT the Step 11 team target.
- Risky sleeve: 8 modeled assets; cash 0.1% is fixed outside the BL system.
- Prior risky weights are normalized to 100% before reverse optimization.
- delta = 2.5, tau = 0.025.
- Team CMA is total return; subtract rf=3.0% so Q is in excess-return units.
- P = I and Omega = Sigma / T, with T = 120 months / 12 = 10 years.
- Sigma_BL = Sigma for the baseline BL allocation.
- Final BL allocation constraints:
  * same policy constraints as earlier steps
  * |active weight| <= 3%p by detailed asset vs official mapped
  * ex-ante tracking error <= 1.0%
- If BL-MVO remains concentrated at bounds, apply Ellipsoidal Robust BL
  using V_BL as expected-return uncertainty, kappa at the Step 8 baseline.
"""

from pathlib import Path
import json
import math
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

RF = 0.03
DELTA = 2.5
TAU = 0.025
GAMMA = 4.0
CASH = 0.001
ACTIVE_CAP = 0.03
TE_CAP = 0.01
KAPPA = float(math.sqrt(chi2.ppf(0.50, df=8)))

OFFICIAL_TOTAL = {
    "국내주식": 0.208,
    "글로벌 선진국 DM": 0.356 * 0.92,
    "글로벌 신흥국 EM": 0.356 * 0.08,
    "국내 국채(초장기)": 0.218,
    "글로벌 IG 크레딧": 0.074,
    "사모주식 PE/VC": 0.143 / 3,
    "실물 인프라": 0.143 / 3,
    "사모대출 PD": 0.143 / 3,
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
mu_total = proxy.set_index("asset").loc[assets, "mu_cma"].astype(float).values
Sigma = cov_cma.loc[assets, assets].astype(float).values
R = rets[assets].astype(float).values
n = len(assets)

T_years = len(R) / 12.0
w_prior_total = np.array([OFFICIAL_TOTAL[a] for a in assets], dtype=float)
assert abs(w_prior_total.sum() - (1.0 - CASH)) < 1e-12

# BL reverse optimization is defined on the risky portfolio, so normalize 99.9% to 100%.
w_prior_risky = w_prior_total / w_prior_total.sum()

# Unit consistency: pi and Q are both excess returns.
mu_excess = mu_total - RF
P = np.eye(n)
pi = DELTA * Sigma @ w_prior_risky
Q = mu_excess.copy()
Omega = P @ (Sigma / T_years) @ P.T

prior_precision = np.linalg.inv(TAU * Sigma)
view_precision = P.T @ np.linalg.inv(Omega) @ P
V_BL = np.linalg.inv(prior_precision + view_precision)

mu_bl_excess = V_BL @ (
    prior_precision @ pi
    + P.T @ np.linalg.inv(Omega) @ Q
)
mu_bl_total = mu_bl_excess + RF

# Closed-form validation under P=I, Omega=Sigma/T.
prior_weight = 1.0 / (1.0 + TAU * T_years)
cma_weight = (TAU * T_years) / (1.0 + TAU * T_years)
mu_bl_closed = prior_weight * pi + cma_weight * Q
if np.max(np.abs(mu_bl_closed - mu_bl_excess)) > 1e-10:
    raise RuntimeError("BL full formula and closed form do not match")

def group_sum(names):
    ids = [idx[a] for a in names]
    return lambda w: w[ids].sum()

policy_cons = [
    {"type": "eq", "fun": lambda w: w.sum() - (1.0 - CASH)},
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

def solve(mu_vec, cov, use_bl_limits=False, robust_penalty=None, x0=None):
    cons = list(policy_cons)

    if use_bl_limits:
        for j in range(n):
            cons.append({
                "type": "ineq",
                "fun": lambda w, j=j: ACTIVE_CAP - (w[j] - w_prior_total[j]),
            })
            cons.append({
                "type": "ineq",
                "fun": lambda w, j=j: ACTIVE_CAP + (w[j] - w_prior_total[j]),
            })

        cons.append({
            "type": "ineq",
            "fun": lambda w: TE_CAP - np.sqrt(
                max(0.0, (w - w_prior_total) @ Sigma @ (w - w_prior_total))
            ),
        })

    if x0 is None:
        x0 = w_prior_total.copy()

    def objective(w):
        utility = mu_vec @ w - GAMMA / 2.0 * (w @ cov @ w)
        if robust_penalty is not None:
            utility -= robust_penalty(w)
        return -utility

    res = minimize(
        objective,
        x0,
        method="SLSQP",
        bounds=bounds,
        constraints=cons,
        options={"ftol": 1e-14, "maxiter": 5000},
    )
    if not res.success:
        raise RuntimeError(res.message)
    return res.x

def metrics(w, mu_eval):
    expected_return = float(mu_eval @ w + RF * CASH)
    volatility = float(np.sqrt(w @ Sigma @ w))
    active = w - w_prior_total
    te = float(np.sqrt(active @ Sigma @ active))
    turnover = float(0.5 * np.abs(active).sum())
    sharpe = float((expected_return - RF) / volatility)
    utility = float(expected_return - GAMMA / 2.0 * volatility**2)
    return {
        "expected_return": expected_return,
        "volatility": volatility,
        "sharpe": sharpe,
        "TE_vs_official": te,
        "turnover_vs_official": turnover,
        "utility": utility,
    }

# Baseline policy MVO for comparison.
w_base = solve(mu_total, Sigma, use_bl_limits=False)

# BL-MVO with assignment-mandated active and TE caps.
w_bl = solve(mu_bl_total, Sigma, use_bl_limits=True, x0=w_prior_total)

# Robust BL because BL-MVO still binds multiple detailed active caps.
robust_penalty = lambda w: KAPPA * np.sqrt(max(0.0, w @ V_BL @ w))
w_robust_bl = solve(
    mu_bl_total,
    Sigma,
    use_bl_limits=True,
    robust_penalty=robust_penalty,
    x0=w_bl,
)

table_d = pd.DataFrame({
    "asset": assets,
    "CMA_total": mu_total,
    "Q_excess": Q,
    "pi_excess": pi,
    "view_error_Q_minus_pi": Q - pi,
    "mu_BL_excess": mu_bl_excess,
    "mu_BL_total": mu_bl_total,
    "tauSigma_diag": np.diag(TAU * Sigma),
    "Omega_diag": np.diag(Omega),
    "V_BL_diag": np.diag(V_BL),
})
table_d.to_csv(RESULTS / "step12_tableD.csv", index=False, encoding="utf-8-sig")

pd.DataFrame(TAU * Sigma, index=assets, columns=assets).to_csv(
    RESULTS / "step12_tauSigma.csv", encoding="utf-8-sig"
)
pd.DataFrame(Omega, index=assets, columns=assets).to_csv(
    RESULTS / "step12_Omega.csv", encoding="utf-8-sig"
)
pd.DataFrame(V_BL, index=assets, columns=assets).to_csv(
    RESULTS / "step12_VBL.csv", encoding="utf-8-sig"
)

alloc = pd.DataFrame({
    "asset": assets,
    "official_mapped_total": w_prior_total,
    "baseline_policy_mvo": w_base,
    "BL_MVO": w_bl,
    "Robust_BL": w_robust_bl,
    "active_BL": w_bl - w_prior_total,
    "active_Robust_BL": w_robust_bl - w_prior_total,
})
alloc.to_csv(RESULTS / "step12_allocations.csv", index=False, encoding="utf-8-sig")

metric_rows = []
for name, w in [
    ("Baseline Policy MVO", w_base),
    ("BL-MVO", w_bl),
    ("Robust BL", w_robust_bl),
]:
    metric_rows.append({"method": name, **metrics(w, mu_bl_total)})
pd.DataFrame(metric_rows).to_csv(
    RESULTS / "step12_metrics.csv", index=False, encoding="utf-8-sig"
)

def reaggregate(w):
    s = pd.Series(w, index=assets)
    return {
        "국내주식": float(s["국내주식"]),
        "해외주식": float(s["글로벌 선진국 DM"] + s["글로벌 신흥국 EM"]),
        "국내채권": float(s["국내 국채(초장기)"]),
        "해외채권": float(s["글로벌 IG 크레딧"]),
        "대체투자": float(s[["사모주식 PE/VC", "실물 인프라", "사모대출 PD"]].sum()),
        "단기자금": CASH,
    }

common = pd.DataFrame({
    "Official 2027": reaggregate(w_prior_total),
    "BL-MVO": reaggregate(w_bl),
    "Robust BL": reaggregate(w_robust_bl),
})
common.to_csv(RESULTS / "step12_common_buckets.csv", encoding="utf-8-sig")

meta = {
    "delta": DELTA,
    "tau": TAU,
    "T_years": T_years,
    "rf": RF,
    "prior_weight": prior_weight,
    "cma_weight": cma_weight,
    "P": "identity",
    "Omega": "Sigma/T",
    "Sigma_BL_choice": "Sigma unchanged",
    "cash_treatment": "0.1% fixed outside 8-asset risky BL system",
    "active_cap": ACTIVE_CAP,
    "TE_cap": TE_CAP,
    "robust_BL_kappa": KAPPA,
    "final_BL_portfolio": "Robust BL",
    "notes": [
        "The Step 11 team target is not used as the BL prior.",
        "BL-MVO hit several +/-3% active caps and the 1.0% TE cap, so Ellipsoidal Robust BL was applied without manual weight edits.",
        "Robust BL still binds some caps; these are reported rather than manually overridden.",
    ],
}
(RESULTS / "step12_meta.json").write_text(
    json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8"
)

print("Prior/CMA weights:", prior_weight, cma_weight)
print("\nTable D")
print(table_d)
print("\nAllocations")
print(alloc)
print("\nMetrics")
print(pd.DataFrame(metric_rows))
