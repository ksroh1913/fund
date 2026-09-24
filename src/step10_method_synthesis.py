# -*- coding: utf-8 -*-
"""STEP 10 - method synthesis.

Table B uses the six assignment methods separately:
Official mapped, Policy MVO, Ledoit-Wolf, Box, Ellipsoid, and *pure* Michaud.

For fair cross-method risk comparison, all portfolios are evaluated with the
common Team CMA covariance Sigma.  For Ledoit-Wolf, model-implied volatility
under the LW covariance is retained as an additional diagnostic.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data"; RESULTS=ROOT/"results"; RESULTS.mkdir(exist_ok=True)
GAMMA=4.0; RF=0.03
OFFICIAL_TOP={"국내주식":.208,"해외주식":.356,"국내채권":.218,"해외채권":.074,"대체투자":.143}
DM_SHARE=.92; EM_SHARE=.08
ALT_SPLIT={"사모주식 PE/VC":1/3,"실물 인프라":1/3,"사모대출 PD":1/3}
DETAIL_BOUNDS={"글로벌 신흥국 EM":(0,.10),"사모주식 PE/VC":(0,.10),"실물 인프라":(0,.10),"사모대출 PD":(0,.10)}

proxy=pd.read_csv(DATA/"proxy_mapping.csv")
rets=pd.read_csv(DATA/"monthly_returns_krw.csv",index_col=0)
cov_cma=pd.read_csv(DATA/"covariance_cma.csv",index_col=0)
assets=proxy["asset"].tolist(); idx={a:i for i,a in enumerate(assets)}
mu=proxy.set_index("asset").loc[assets,"mu_cma"].astype(float).values
sigma=proxy.set_index("asset").loc[assets,"sigma_cma"].astype(float).values
Sigma=cov_cma.loc[assets,assets].astype(float).values
R=rets[assets].astype(float).values; n=len(assets)

scale=1/sum(OFFICIAL_TOP.values())
w_off=pd.Series(0.,index=assets)
w_off["국내주식"]=OFFICIAL_TOP["국내주식"]*scale
w_off["글로벌 선진국 DM"]=OFFICIAL_TOP["해외주식"]*DM_SHARE*scale
w_off["글로벌 신흥국 EM"]=OFFICIAL_TOP["해외주식"]*EM_SHARE*scale
w_off["국내 국채(초장기)"]=OFFICIAL_TOP["국내채권"]*scale
w_off["글로벌 IG 크레딧"]=OFFICIAL_TOP["해외채권"]*scale
for a,s in ALT_SPLIT.items(): w_off[a]=OFFICIAL_TOP["대체투자"]*s*scale

def ids(names): return [idx[a] for a in names]
EQ=ids(["국내주식","글로벌 선진국 DM","글로벌 신흥국 EM"])
BD=ids(["국내 국채(초장기)","글로벌 IG 크레딧"])
ALT=ids(["사모주식 PE/VC","실물 인프라","사모대출 PD"])
def sm(ii): return lambda w:w[ii].sum()
policy=[
 {"type":"eq","fun":lambda w:w.sum()-1},
 {"type":"ineq","fun":lambda w:w[idx["국내주식"]]-.15},{"type":"ineq","fun":lambda w:.25-w[idx["국내주식"]]},
 {"type":"ineq","fun":lambda w:sm(ids(["글로벌 선진국 DM","글로벌 신흥국 EM"]))(w)-.30},
 {"type":"ineq","fun":lambda w:.42-sm(ids(["글로벌 선진국 DM","글로벌 신흥국 EM"]))(w)},
 {"type":"ineq","fun":lambda w:w[idx["국내 국채(초장기)"]]-.15},{"type":"ineq","fun":lambda w:.25-w[idx["국내 국채(초장기)"]]},
 {"type":"ineq","fun":lambda w:w[idx["글로벌 IG 크레딧"]]-.05},{"type":"ineq","fun":lambda w:.12-w[idx["글로벌 IG 크레딧"]]},
 {"type":"ineq","fun":lambda w:sm(ALT)(w)-.12},{"type":"ineq","fun":lambda w:.18-sm(ALT)(w)},
 {"type":"ineq","fun":lambda w:sm(EQ)(w)-.52},{"type":"ineq","fun":lambda w:.60-sm(EQ)(w)},
 {"type":"ineq","fun":lambda w:sm(BD)(w)-.25},{"type":"ineq","fun":lambda w:.35-sm(BD)(w)}
]
bounds=[(0,1)]*n
for a,b in DETAIL_BOUNDS.items(): bounds[idx[a]]=b
def solve(m,C,x0=None):
    if x0 is None:x0=w_off.values.copy()
    r=minimize(lambda w:-(m@w-GAMMA/2*(w@C@w)),x0,method="SLSQP",bounds=bounds,constraints=policy,
               options={"ftol":1e-12,"maxiter":5000})
    if not r.success: raise RuntimeError(r.message)
    return r.x

def lw_cov_corr(X):
    Xc=X-X.mean(0); T=len(X); S0=Xc.T@Xc/T
    target=np.trace(S0)/n*np.eye(n); d2=np.sum((S0-target)**2)
    b2=sum(np.sum((np.outer(x,x)-S0)**2) for x in Xc)/T**2
    b2=min(b2,d2); delta=b2/d2
    Sl=delta*target+(1-delta)*S0; sd=np.sqrt(np.diag(Sl))
    return Sl*12, Sl/np.outer(sd,sd), float(delta)

S_lw_hist,rho_lw,lw_delta=lw_cov_corr(R)
Sigma_lw=np.outer(sigma,sigma)*rho_lw
w_mvo=solve(mu,Sigma)
w_lw=solve(mu,Sigma_lw,w_mvo)
box=pd.read_csv(RESULTS/"step7_box.csv"); w_box=box.set_index("자산").loc[assets,"Box비중"].values
ell=pd.read_csv(RESULTS/"step8_ellipsoid.csv"); erow=ell.iloc[(ell["p"]-.5).abs().argsort()[:1]]
w_ell=erow[assets].iloc[0].values
mich=pd.read_csv(RESULTS/"step9_michaud_summary.csv"); w_mich=mich.set_index("asset").loc[assets,"michaud_mean"].values
methods={"Official mapped":w_off.values,"Policy MVO":w_mvo,"Ledoit-Wolf":w_lw,
         "Box Robust":w_box,"Ellipsoid Robust":w_ell,"Michaud mean":w_mich}
pd.DataFrame(methods,index=assets).to_csv(RESULTS/"step10_tableB_weights.csv",encoding="utf-8-sig")

def metric(w,C=Sigma):
    er=float(mu@w); vol=float(np.sqrt(w@C@w)); a=w-w_off.values
    return {"expected_return":er,"volatility":vol,"sharpe":float((er-RF)/vol),
            "TE":float(np.sqrt(a@Sigma@a)),"turnover":float(.5*np.abs(a).sum())}
rows=[]
for name,w in methods.items():
    common=metric(w,Sigma)
    row={"method":name,**common}
    if name=="Ledoit-Wolf":
        lm=metric(w,Sigma_lw)
        row["model_cov_volatility"]=lm["volatility"]; row["model_cov_sharpe"]=lm["sharpe"]
    else:
        row["model_cov_volatility"]=common["volatility"]; row["model_cov_sharpe"]=common["sharpe"]
    rows.append(row)
pd.DataFrame(rows).to_csv(RESULTS/"step10_tableB_metrics.csv",index=False,encoding="utf-8-sig")

re=[]
for name,w in methods.items():
    s=pd.Series(w,index=assets)
    re.append({"method":name,"국내주식":s["국내주식"],"해외주식":s["글로벌 선진국 DM"]+s["글로벌 신흥국 EM"],
               "국내채권":s["국내 국채(초장기)"],"해외채권":s["글로벌 IG 크레딧"],
               "대체투자":s[["사모주식 PE/VC","실물 인프라","사모대출 PD"]].sum()})
pd.DataFrame(re).to_csv(RESULTS/"step10_reaggregated.csv",index=False,encoding="utf-8-sig")

sens=[]
for j,a in enumerate(assets):
    for shock in (-.005,.005):
        m2=mu.copy();m2[j]+=shock;w=solve(m2,Sigma,w_mvo)
        sens.append({"asset_shocked":a,"shock":shock,"own_weight_base":w_mvo[j],"own_weight_new":w[j],
                     "own_weight_change":w[j]-w_mvo[j],"portfolio_L1_change":float(np.abs(w-w_mvo).sum()),
                     "portfolio_turnover_change":float(.5*np.abs(w-w_mvo).sum())})
sdf=pd.DataFrame(sens); sdf.to_csv(RESULTS/"step10_mu_50bp_sensitivity.csv",index=False,encoding="utf-8-sig")
summ=[]
for a,g in sdf.groupby("asset_shocked",sort=False):
    io=g["own_weight_change"].abs().idxmax(); il=g["portfolio_L1_change"].idxmax()
    summ.append({"asset":a,"max_abs_own_weight_change":abs(g.loc[io,"own_weight_change"]),
                 "shock_for_max_own":g.loc[io,"shock"],"max_portfolio_L1_change":g.loc[il,"portfolio_L1_change"],
                 "shock_for_max_L1":g.loc[il,"shock"]})
pd.DataFrame(summ).sort_values("max_abs_own_weight_change",ascending=False).to_csv(
    RESULTS/"step10_mu_50bp_sensitivity_summary.csv",index=False,encoding="utf-8-sig")

ev,V=np.linalg.eigh(Sigma); order=np.argsort(ev);ev=ev[order];V=V[:,order]
erows=[]
for k in range(min(4,n)):
    v=V[:,k]/np.max(np.abs(V[:,k]))
    erows.append({"rank_smallest":k+1,"eigenvalue":ev[k],**{assets[i]:v[i] for i in range(n)}})
pd.DataFrame(erows).to_csv(RESULTS/"step10_small_eigenvectors.csv",index=False,encoding="utf-8-sig")
v=V[:,0]/np.max(np.abs(V[:,0]))
meta={"baseline_gamma":GAMMA,"lw_shrinkage":lw_delta,"smallest_eigenvalue":float(ev[0]),
      "smallest_eigenvector":{assets[i]:float(v[i]) for i in range(n)},
      "primary_reference_method":"Pure Michaud resampling",
      "primary_reference_reason":"Separates resampling instability from Ledoit-Wolf covariance shrinkage; LW+Michaud is retained as a sensitivity.",
      "gamma_reason":"gamma=2 and gamma=4 produce nearly the same ~10.1% policy-MVO risk, close to the official-mapped ~10.7% risk; gamma=6 is a lower-risk sensitivity.",
      "risk_comparison":"Table B volatility and Sharpe use common Team CMA Sigma; LW model-covariance metrics are shown separately."}
(RESULTS/"step10_meta.json").write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding="utf-8")
print((pd.DataFrame(methods,index=assets)*100).round(2))
