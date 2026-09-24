# -*- coding: utf-8 -*-
"""Reproduce Team 2 NPS assignment Steps 5-8 from CSV inputs."""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import chi2
from sklearn.covariance import LedoitWolf

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

SEED = 60080
RF = 0.03
GAMMAS = [2.0, 4.0, 6.0]
BASE_GAMMA = 4.0
BOOT_B = 5000

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
    "사모주식 PE/VC": 1/3,
    "실물 인프라": 1/3,
    "사모대출 PD": 1/3,
}

# Team analysis assumptions; not official NPS detailed limits.
DETAIL_BOUNDS = {
    "글로벌 신흥국 EM": (0.0, 0.10),
    "사모주식 PE/VC": (0.0, 0.10),
    "실물 인프라": (0.0, 0.10),
    "사모대출 PD": (0.0, 0.10),
}

BOX_D = {
    "국내주식": 0.010,
    "글로벌 선진국 DM": 0.010,
    "글로벌 신흥국 EM": 0.010,
    "국내 국채(초장기)": 0.005,
    "글로벌 IG 크레딧": 0.005,
    "사모주식 PE/VC": 0.015,
    "실물 인프라": 0.015,
    "사모대출 PD": 0.015,
}

proxy = pd.read_csv(DATA / "proxy_mapping.csv")
rets = pd.read_csv(DATA / "monthly_returns_krw.csv", index_col=0)
cov_cma = pd.read_csv(DATA / "covariance_cma.csv", index_col=0)

assets = proxy["asset"].tolist()
mu = proxy.set_index("asset").loc[assets, "mu_cma"].astype(float).values
sigma = proxy.set_index("asset").loc[assets, "sigma_cma"].astype(float).values
Sigma = cov_cma.loc[assets, assets].astype(float).values
R = rets[assets].astype(float).values
n = len(assets)
idx = {a:i for i,a in enumerate(assets)}

# Official target mapped to the 8-asset Team 2 structure.
scale = 1.0 / sum(OFFICIAL_TOP.values())
w_off = pd.Series(0.0, index=assets)
w_off["국내주식"] = OFFICIAL_TOP["국내주식"] * scale
w_off["글로벌 선진국 DM"] = OFFICIAL_TOP["해외주식"] * DM_SHARE * scale
w_off["글로벌 신흥국 EM"] = OFFICIAL_TOP["해외주식"] * EM_SHARE * scale
w_off["국내 국채(초장기)"] = OFFICIAL_TOP["국내채권"] * scale
w_off["글로벌 IG 크레딧"] = OFFICIAL_TOP["해외채권"] * scale
for a, share in ALT_SPLIT.items():
    w_off[a] = OFFICIAL_TOP["대체투자"] * share * scale


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

bounds = [(0.0,1.0)] * n
for a, b in DETAIL_BOUNDS.items():
    bounds[idx[a]] = b


def solve(mu_vec, cov, gamma, robust_penalty=None):
    def objective(w):
        utility = mu_vec @ w - gamma/2 * (w @ cov @ w)
        if robust_penalty is not None:
            utility -= robust_penalty(w)
        return -utility

    res = minimize(
        objective, w_off.values.copy(), method="SLSQP",
        bounds=bounds, constraints=policy_cons,
        options={"ftol":1e-12, "maxiter":5000}
    )
    if not res.success:
        raise RuntimeError(res.message)
    return res.x


def solve_simple(mu_vec, cov, gamma):
    cons = [{"type":"eq","fun":lambda w:w.sum()-1.0}]
    res = minimize(
        lambda w:-(mu_vec@w-gamma/2*(w@cov@w)),
        np.repeat(1/n,n), method="SLSQP",
        bounds=[(0,1)]*n, constraints=cons,
        options={"ftol":1e-12,"maxiter":5000}
    )
    if not res.success:
        raise RuntimeError(res.message)
    return res.x


def metrics(w, risk_cov=Sigma):
    er = float(mu @ w)
    vol = float(np.sqrt(w @ risk_cov @ w))
    d = w - w_off.values
    return {
        "expected_return": er,
        "volatility": vol,
        "sharpe": float((er-RF)/vol),
        "TE": float(np.sqrt(d @ Sigma @ d)),
        "turnover": float(0.5*np.abs(d).sum()),
    }


# STEP 5
rows = []
for g in GAMMAS:
    for label, w in [
        ("단순 MVO", solve_simple(mu,Sigma,g)),
        ("정책 MVO", solve(mu,Sigma,g)),
    ]:
        m = metrics(w)
        rows.append({
            "구분":label, "gamma":g,
            **dict(zip(assets,w)),
            "기대수익률":m["expected_return"],
            "변동성":m["volatility"],
            "Sharpe":m["sharpe"],
            "TE":m["TE"],
            "회전율":m["turnover"],
        })
pd.DataFrame(rows).to_csv(RESULTS/"step5_mvo.csv",index=False,encoding="utf-8-sig")

# STEP 6
lw = LedoitWolf().fit(R)
S_sample = np.cov(R,rowvar=False,ddof=0) * 12
S_lw = lw.covariance_ * 12
std_lw = np.sqrt(np.diag(S_lw))
rho_lw = S_lw / np.outer(std_lw,std_lw)
Sigma_lw_cma = np.outer(sigma,sigma) * rho_lw
w_lw = solve(mu,Sigma_lw_cma,BASE_GAMMA)

pd.DataFrame([
    ["표본 공분산",np.linalg.eigvalsh(S_sample).min(),np.linalg.cond(S_sample)],
    ["Ledoit-Wolf",np.linalg.eigvalsh(S_lw).min(),np.linalg.cond(S_lw)],
    ["CMA Sigma",np.linalg.eigvalsh(Sigma).min(),np.linalg.cond(Sigma)],
    ["LW 상관+CMA sigma",np.linalg.eigvalsh(Sigma_lw_cma).min(),np.linalg.cond(Sigma_lw_cma)],
],columns=["구분","최소고유값","조건수"]).to_csv(
    RESULTS/"step6_lw_diagnostics.csv",index=False,encoding="utf-8-sig"
)
print("LW shrinkage:",lw.shrinkage_)
print("LW weights:",dict(zip(assets,w_lw)))
print("LW metrics:",metrics(w_lw,Sigma_lw_cma))

# STEP 7
d = np.array([BOX_D[a] for a in assets])
mu_worst = mu-d
w_box = solve(mu_worst,Sigma,BASE_GAMMA)
pd.DataFrame({
    "자산":assets,
    "CMA_mu":mu,
    "Box폭":d,
    "최악_mu":mu_worst,
    "Box비중":w_box,
}).to_csv(RESULTS/"step7_box.csv",index=False,encoding="utf-8-sig")

# STEP 8
rng = np.random.default_rng(SEED)
T = len(R)
boot = np.empty((BOOT_B,n))
for b in range(BOOT_B):
    sample = R[rng.integers(0,T,T)]
    boot[b] = sample.mean(axis=0)*12

S_mu = np.cov(boot,rowvar=False,ddof=1)
pd.DataFrame(S_mu,index=assets,columns=assets).to_csv(
    RESULTS/"step8_S_mu.csv",encoding="utf-8-sig"
)

ell_rows = []
for p in [0.10,0.50,0.90]:
    kappa = float(np.sqrt(chi2.ppf(p,df=n)))
    penalty = lambda w,k=kappa:k*np.sqrt(max(0.0,w@S_mu@w))
    w = solve(mu,Sigma,BASE_GAMMA,penalty)
    m = metrics(w)
    mu_se = float(np.sqrt(w@S_mu@w))
    ell_rows.append({
        "p":p, "kappa":kappa,
        **dict(zip(assets,w)),
        "기대수익률":m["expected_return"],
        "변동성":m["volatility"],
        "Sharpe":m["sharpe"],
        "TE":m["TE"],
        "회전율":m["turnover"],
        "portfolio_mu_SE":mu_se,
        "worst_mu":float(mu@w-kappa*mu_se),
    })

pd.DataFrame(ell_rows).to_csv(
    RESULTS/"step8_ellipsoid.csv",index=False,encoding="utf-8-sig"
)
print("Saved results to",RESULTS)
