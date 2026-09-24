# -*- coding: utf-8 -*-
"""STEP 9 - Michaud resampling.

Baseline required for Table B:
- B=300, seed=60080.
- Bootstrap the frozen 120 monthly KRW proxy returns.
- Recenter the bootstrap mean shock around the forward-looking Team CMA.
- Re-estimate the *sample correlation* in every draw and re-scale with CMA sigma.
- Apply exactly the same gamma=4 policy constraints as the baseline MVO.

A separate LW+Michaud sensitivity is retained because it is a useful practical
combination, but Table B uses the pure Michaud result so the resampling effect
is not mixed with Ledoit-Wolf shrinkage.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

SEED, B, GAMMA, RF = 60080, 300, 4.0, 0.03

OFFICIAL_TOP = {
    "국내주식": 0.208, "해외주식": 0.356, "국내채권": 0.218,
    "해외채권": 0.074, "대체투자": 0.143,
}
DM_SHARE, EM_SHARE = 0.92, 0.08
ALT_SPLIT = {"사모주식 PE/VC":1/3, "실물 인프라":1/3, "사모대출 PD":1/3}
DETAIL_BOUNDS = {
    "글로벌 신흥국 EM":(0.0,0.10), "사모주식 PE/VC":(0.0,0.10),
    "실물 인프라":(0.0,0.10), "사모대출 PD":(0.0,0.10),
}

proxy = pd.read_csv(DATA/"proxy_mapping.csv")
rets = pd.read_csv(DATA/"monthly_returns_krw.csv", index_col=0)
cov_cma = pd.read_csv(DATA/"covariance_cma.csv", index_col=0)
assets = proxy["asset"].tolist()
mu = proxy.set_index("asset").loc[assets,"mu_cma"].astype(float).values
sigma = proxy.set_index("asset").loc[assets,"sigma_cma"].astype(float).values
Sigma = cov_cma.loc[assets,assets].astype(float).values
R = rets[assets].astype(float).values
T,n = R.shape
idx = {a:i for i,a in enumerate(assets)}

scale = 1.0/sum(OFFICIAL_TOP.values())
w_off = pd.Series(0.0,index=assets)
w_off["국내주식"] = OFFICIAL_TOP["국내주식"]*scale
w_off["글로벌 선진국 DM"] = OFFICIAL_TOP["해외주식"]*DM_SHARE*scale
w_off["글로벌 신흥국 EM"] = OFFICIAL_TOP["해외주식"]*EM_SHARE*scale
w_off["국내 국채(초장기)"] = OFFICIAL_TOP["국내채권"]*scale
w_off["글로벌 IG 크레딧"] = OFFICIAL_TOP["해외채권"]*scale
for a,s in ALT_SPLIT.items(): w_off[a] = OFFICIAL_TOP["대체투자"]*s*scale

def ids(names): return [idx[a] for a in names]
EQ=ids(["국내주식","글로벌 선진국 DM","글로벌 신흥국 EM"])
BD=ids(["국내 국채(초장기)","글로벌 IG 크레딧"])
ALT=ids(["사모주식 PE/VC","실물 인프라","사모대출 PD"])

def gsum(ii): return lambda w: w[ii].sum()
policy_cons=[
    {"type":"eq","fun":lambda w:w.sum()-1.0},
    {"type":"ineq","fun":lambda w:w[idx["국내주식"]]-.15},
    {"type":"ineq","fun":lambda w:.25-w[idx["국내주식"]]},
    {"type":"ineq","fun":lambda w:gsum(ids(["글로벌 선진국 DM","글로벌 신흥국 EM"]))(w)-.30},
    {"type":"ineq","fun":lambda w:.42-gsum(ids(["글로벌 선진국 DM","글로벌 신흥국 EM"]))(w)},
    {"type":"ineq","fun":lambda w:w[idx["국내 국채(초장기)"]]-.15},
    {"type":"ineq","fun":lambda w:.25-w[idx["국내 국채(초장기)"]]},
    {"type":"ineq","fun":lambda w:w[idx["글로벌 IG 크레딧"]]-.05},
    {"type":"ineq","fun":lambda w:.12-w[idx["글로벌 IG 크레딧"]]},
    {"type":"ineq","fun":lambda w:gsum(ALT)(w)-.12},
    {"type":"ineq","fun":lambda w:.18-gsum(ALT)(w)},
    {"type":"ineq","fun":lambda w:gsum(EQ)(w)-.52},
    {"type":"ineq","fun":lambda w:.60-gsum(EQ)(w)},
    {"type":"ineq","fun":lambda w:gsum(BD)(w)-.25},
    {"type":"ineq","fun":lambda w:.35-gsum(BD)(w)},
]
bounds=[(0.0,1.0)]*n
for a,b in DETAIL_BOUNDS.items(): bounds[idx[a]]=b

def solve(mu_vec,cov,x0=None,bounds_override=None):
    bb=bounds if bounds_override is None else bounds_override
    if x0 is None: x0=w_off.values.copy()
    res=minimize(lambda w:-(mu_vec@w-GAMMA/2*(w@cov@w)),x0,method="SLSQP",
                 bounds=bb,constraints=policy_cons,
                 options={"ftol":1e-12,"maxiter":5000})
    if not res.success: raise RuntimeError(res.message)
    return res.x

def metrics(w):
    er=float(mu@w); vol=float(np.sqrt(w@Sigma@w)); a=w-w_off.values
    return {"expected_return":er,"volatility":vol,"sharpe":float((er-RF)/vol),
            "TE":float(np.sqrt(a@Sigma@a)),"turnover":float(.5*np.abs(a).sum())}

def lw_corr(X):
    """Ledoit-Wolf 2004 scaled-identity correlation; sklearn convention."""
    Xc=X-X.mean(0); tt=len(X); Sm=Xc.T@Xc/tt
    m=np.trace(Sm)/n; F=m*np.eye(n); d2=np.sum((Sm-F)**2)
    b2=sum(np.sum((np.outer(x,x)-Sm)**2) for x in Xc)/tt**2
    b2=min(b2,d2); delta=b2/d2
    Sl=delta*F+(1-delta)*Sm; sd=np.sqrt(np.diag(Sl))
    return Sl/np.outer(sd,sd), float(delta)

def run_resampling(use_lw=False, bounds_override=None):
    rng=np.random.default_rng(SEED); full=R.mean(0); W=[]; M=[]; shrink=[]; prev=None
    for _ in range(B):
        rb=R[rng.integers(0,T,T)]
        mb=mu+12.0*(rb.mean(0)-full)
        if use_lw:
            corr,delta=lw_corr(rb); shrink.append(delta)
        else:
            corr=np.corrcoef(rb,rowvar=False)
        cov=np.outer(sigma,sigma)*corr
        try: wb=solve(mb,cov,prev,bounds_override)
        except RuntimeError: wb=solve(mb,cov,w_off.values.copy(),bounds_override)
        W.append(wb); M.append(mb); prev=wb
    return np.asarray(W),np.asarray(M),np.asarray(shrink)

w_base=solve(mu,Sigma)
W,mu_draws,_=run_resampling(False)
W_lw,_,shrink=run_resampling(True)

def summary_frame(W):
    d=pd.DataFrame(W,columns=assets)
    s=pd.DataFrame({
        "baseline_mvo":w_base,
        "michaud_mean":d.mean().values,
        "p05":d.quantile(.05).values,
        "p50":d.quantile(.50).values,
        "p95":d.quantile(.95).values,
        "std":d.std(ddof=1).values,
    },index=assets)
    s["change_vs_baseline"]=s["michaud_mean"]-s["baseline_mvo"]
    upper={a:DETAIL_BOUNDS.get(a,(0,1))[1] for a in assets}
    s["freq_at_zero"]=[float((d[a]<=1e-6).mean()) for a in assets]
    s["freq_at_upper"]=[float((d[a]>=upper[a]-1e-6).mean()) if upper[a]<1 else np.nan for a in assets]
    return s

pure=summary_frame(W); lw_sens=summary_frame(W_lw)
pure.reset_index(names="asset").to_csv(RESULTS/"step9_michaud_summary.csv",index=False,encoding="utf-8-sig")
lw_sens.reset_index(names="asset").to_csv(RESULTS/"step9_michaud_lw_sensitivity.csv",index=False,encoding="utf-8-sig")
pd.DataFrame(W,columns=assets).to_csv(RESULTS/"step9_michaud_draws.csv",index_label="iteration",encoding="utf-8-sig")
pd.DataFrame(mu_draws,columns=assets).to_csv(RESULTS/"step9_mu_draws.csv",index=False,encoding="utf-8-sig")

wm=pure["michaud_mean"].values
meta={
    "B":B,"seed":SEED,"T":T,"gamma":GAMMA,
    "baseline_method":"pure Michaud: bootstrap mean + per-draw sample correlation + CMA sigma",
    "sensitivity_method":"Michaud + per-draw Ledoit-Wolf correlation",
    "baseline_metrics":metrics(w_base),"michaud_mean_metrics":metrics(wm),
    "lw_sensitivity_mean_metrics":metrics(lw_sens["michaud_mean"].values),
    "average_lw_shrinkage":float(shrink.mean()),
}
(RESULTS/"step9_michaud_meta.json").write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding="utf-8")
print((pure*100).round(2))
print(json.dumps(meta,ensure_ascii=False,indent=2))
