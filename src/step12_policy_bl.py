# -*- coding: utf-8 -*-
"""STEP 12 - Policy Black-Litterman with assignment sensitivities.

Baseline
--------
- Prior = official 2027 mapped weights (never the Team target).
- 8 modeled risky assets; policy cash 0.1% fixed.
- delta=2.5, tau=0.025, rf=3%, gamma=4.
- P=I, Q=Team CMA excess returns.
- T=10 is a disclosed proxy assumption: the forward-looking W3 CMA is not a
  sample mean, so the 120-month risk/correlation sample is converted to a
  10-year equivalent solely to operationalize the assignment's Omega=Sigma/T.
- BL allocation: same policy constraints, detailed active <=3%p, TE<=1%.
- Robust BL: Ellipsoidal penalty around mu_BL using V_BL.

Sensitivities
-------------
1) T = 5, 10, 20, recomputing mu_BL, BL-MVO and Robust BL.
2) TE cap = 0.95% (baseline remains the assignment-required 1.0%).
3) Cash optimized in [0,2%] with mu_cash=rf instead of fixed 0.1%.
"""
from pathlib import Path
import json, math
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import chi2

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data"; RESULTS=ROOT/"results"; RESULTS.mkdir(exist_ok=True)
RF=.03; DELTA=2.5; TAU=.025; GAMMA=4.; CASH=.001; ACTIVE_CAP=.03; TE_CAP=.01
KAPPA=float(math.sqrt(chi2.ppf(.50,df=8)))
OFFICIAL_TOTAL={
 "국내주식":.208,"글로벌 선진국 DM":.356*.92,"글로벌 신흥국 EM":.356*.08,
 "국내 국채(초장기)":.218,"글로벌 IG 크레딧":.074,
 "사모주식 PE/VC":.143/3,"실물 인프라":.143/3,"사모대출 PD":.143/3,
}
DETAIL_BOUNDS={"글로벌 신흥국 EM":(0,.10),"사모주식 PE/VC":(0,.10),
               "실물 인프라":(0,.10),"사모대출 PD":(0,.10)}

proxy=pd.read_csv(DATA/"proxy_mapping.csv")
rets=pd.read_csv(DATA/"monthly_returns_krw.csv",index_col=0)
cov=pd.read_csv(DATA/"covariance_cma.csv",index_col=0)
assets=proxy["asset"].tolist();idx={a:i for i,a in enumerate(assets)};n=len(assets)
mu_total=proxy.set_index("asset").loc[assets,"mu_cma"].astype(float).values
Sigma=cov.loc[assets,assets].astype(float).values
R=rets[assets].astype(float).values
w_prior_total=np.array([OFFICIAL_TOTAL[a] for a in assets],float)
assert abs(w_prior_total.sum()-(1-CASH))<1e-12
w_prior_risky=w_prior_total/w_prior_total.sum()
mu_excess=mu_total-RF; pi=DELTA*Sigma@w_prior_risky; Q=mu_excess.copy()

def ids(names): return [idx[a] for a in names]
EQ=ids(["국내주식","글로벌 선진국 DM","글로벌 신흥국 EM"])
BD=ids(["국내 국채(초장기)","글로벌 IG 크레딧"])
ALT=ids(["사모주식 PE/VC","실물 인프라","사모대출 PD"])
FX=ids(["글로벌 선진국 DM","글로벌 신흥국 EM"])
def sm(ii): return lambda w:w[ii].sum()

bounds8=[(0,1)]*n
for a,b in DETAIL_BOUNDS.items(): bounds8[idx[a]]=b

def constraints8(total, bench=None, active_cap=None, te_cap=None):
    c=[
      {"type":"eq","fun":lambda w,t=total:w.sum()-t},
      {"type":"ineq","fun":lambda w:w[idx["국내주식"]]-.15},{"type":"ineq","fun":lambda w:.25-w[idx["국내주식"]]},
      {"type":"ineq","fun":lambda w:sm(FX)(w)-.30},{"type":"ineq","fun":lambda w:.42-sm(FX)(w)},
      {"type":"ineq","fun":lambda w:w[idx["국내 국채(초장기)"]]-.15},{"type":"ineq","fun":lambda w:.25-w[idx["국내 국채(초장기)"]]},
      {"type":"ineq","fun":lambda w:w[idx["글로벌 IG 크레딧"]]-.05},{"type":"ineq","fun":lambda w:.12-w[idx["글로벌 IG 크레딧"]]},
      {"type":"ineq","fun":lambda w:sm(ALT)(w)-.12},{"type":"ineq","fun":lambda w:.18-sm(ALT)(w)},
      {"type":"ineq","fun":lambda w:sm(EQ)(w)-.52},{"type":"ineq","fun":lambda w:.60-sm(EQ)(w)},
      {"type":"ineq","fun":lambda w:sm(BD)(w)-.25},{"type":"ineq","fun":lambda w:.35-sm(BD)(w)},
    ]
    if bench is not None and active_cap is not None:
        for j in range(n):
            c += [
              {"type":"ineq","fun":lambda w,j=j:active_cap-(w[j]-bench[j])},
              {"type":"ineq","fun":lambda w,j=j:active_cap+(w[j]-bench[j])},
            ]
    if bench is not None and te_cap is not None:
        c.append({"type":"ineq","fun":lambda w:te_cap-np.sqrt(max(0.,(w-bench)@Sigma@(w-bench)))})
    return c

def solve8(m,C,total=.999,bench=None,active_cap=None,te_cap=None,U=None,kappa=0.,x0=None):
    if x0 is None: x0=w_prior_total.copy()*(total/.999)
    def obj(w):
        u=m@w-GAMMA/2*(w@C@w)
        if U is not None and kappa: u-=kappa*np.sqrt(max(0.,w@U@w))
        return -u
    r=minimize(obj,x0,method="SLSQP",bounds=bounds8,
               constraints=constraints8(total,bench,active_cap,te_cap),
               options={"ftol":1e-14,"maxiter":10000})
    if not r.success: raise RuntimeError(r.message)
    return r.x

def bl_for_T(T):
    Omega=Sigma/T
    pp=np.linalg.inv(TAU*Sigma); vp=np.linalg.inv(Omega)
    V=np.linalg.inv(pp+vp)
    m_ex=V@(pp@pi+vp@Q)
    pwt=1/(1+TAU*T); cwt=TAU*T/(1+TAU*T)
    closed=pwt*pi+cwt*Q
    if np.max(np.abs(closed-m_ex))>1e-10: raise RuntimeError("BL closed-form mismatch")
    return Omega,V,m_ex+RF,pwt,cwt

def alloc_for_T(T,te_cap=.01):
    Omega,V,mbl,pwt,cwt=bl_for_T(T)
    wb=solve8(mbl,Sigma,bench=w_prior_total,active_cap=ACTIVE_CAP,te_cap=te_cap)
    wr=solve8(mbl,Sigma,bench=w_prior_total,active_cap=ACTIVE_CAP,te_cap=te_cap,U=V,kappa=KAPPA,x0=wb)
    return Omega,V,mbl,wb,wr,pwt,cwt

def metrics8(w,m,cash=CASH):
    er=float(m@w+RF*cash); vol=float(np.sqrt(w@Sigma@w)); a=w-w_prior_total
    return {"expected_return":er,"volatility":vol,"sharpe":float((er-RF)/vol),
            "TE_vs_official":float(np.sqrt(a@Sigma@a)),"turnover_vs_official":float(.5*np.abs(a).sum()+.5*abs(cash-CASH)),
            "utility":float(er-GAMMA/2*vol**2)}

Omega,V_BL,mu_bl_total,w_bl,w_rbl,prior_weight,cma_weight=alloc_for_T(10.)
w_base=solve8(mu_total,Sigma)

# Baseline assignment outputs.
pd.DataFrame({"asset":assets,"CMA_total":mu_total,"Q_excess":Q,"pi_excess":pi,
 "view_error_Q_minus_pi":Q-pi,"mu_BL_excess":mu_bl_total-RF,"mu_BL_total":mu_bl_total,
 "tauSigma_diag":np.diag(TAU*Sigma),"Omega_diag":np.diag(Omega),"V_BL_diag":np.diag(V_BL)}
).to_csv(RESULTS/"step12_tableD.csv",index=False,encoding="utf-8-sig")
pd.DataFrame(TAU*Sigma,index=assets,columns=assets).to_csv(RESULTS/"step12_tauSigma.csv",encoding="utf-8-sig")
pd.DataFrame(Omega,index=assets,columns=assets).to_csv(RESULTS/"step12_Omega.csv",encoding="utf-8-sig")
pd.DataFrame(V_BL,index=assets,columns=assets).to_csv(RESULTS/"step12_VBL.csv",encoding="utf-8-sig")
pd.DataFrame({"asset":assets,"official_mapped_total":w_prior_total,"baseline_policy_mvo":w_base,
 "BL_MVO":w_bl,"Robust_BL":w_rbl,"active_BL":w_bl-w_prior_total,"active_Robust_BL":w_rbl-w_prior_total}
).to_csv(RESULTS/"step12_allocations.csv",index=False,encoding="utf-8-sig")
pd.DataFrame([{"method":"Baseline Policy MVO",**metrics8(w_base,mu_bl_total)},
              {"method":"BL-MVO",**metrics8(w_bl,mu_bl_total)},
              {"method":"Robust BL",**metrics8(w_rbl,mu_bl_total)}]
).to_csv(RESULTS/"step12_metrics.csv",index=False,encoding="utf-8-sig")

def agg(w,cash=CASH):
    s=pd.Series(w,index=assets)
    return {"국내주식":float(s["국내주식"]),"해외주식":float(s["글로벌 선진국 DM"]+s["글로벌 신흥국 EM"]),
            "국내채권":float(s["국내 국채(초장기)"]),"해외채권":float(s["글로벌 IG 크레딧"]),
            "대체투자":float(s[["사모주식 PE/VC","실물 인프라","사모대출 PD"]].sum()),"단기자금":cash}
pd.DataFrame({"Official 2027":agg(w_prior_total),"BL-MVO":agg(w_bl),"Robust BL":agg(w_rbl)}
).to_csv(RESULTS/"step12_common_buckets.csv",encoding="utf-8-sig")

# T sensitivity: every T recomputes posterior and both allocations.
trows=[]
for T in (5.,10.,20.):
    Om,V,mbl,wb,wr,pwt,cwt=alloc_for_T(T)
    for nm,w in [("BL-MVO",wb),("Robust BL",wr)]:
        trows.append({"T_years":T,"prior_weight":pwt,"cma_weight":cwt,"method":nm,
                      **{assets[i]:w[i] for i in range(n)},**metrics8(w,mbl)})
pd.DataFrame(trows).to_csv(RESULTS/"step12_T_sensitivity.csv",index=False,encoding="utf-8-sig")

# TE 0.95% sensitivity, while 1.0% remains the assignment baseline.
_,_,m95,wb95,wr95,_,_=alloc_for_T(10.,te_cap=.0095)
pd.DataFrame([{"method":"BL-MVO_TE095",**{assets[i]:wb95[i] for i in range(n)},**metrics8(wb95,m95)},
              {"method":"Robust_BL_TE095",**{assets[i]:wr95[i] for i in range(n)},**metrics8(wr95,m95)}]
).to_csv(RESULTS/"step12_TE095_sensitivity.csv",index=False,encoding="utf-8-sig")

# Cash sensitivity: optimize a ninth cash asset with mu=rf, sigma=0, cash in [0,2%].
def solve_cash(mbl,V=None,kappa=0.,te_cap=.01,use_bl_limits=True):
    x0=np.r_[w_prior_total,CASH]
    b9=bounds8+[(0,.02)]
    def obj(x):
        w=x[:n];cash=x[-1]
        u=mbl@w+RF*cash-GAMMA/2*(w@Sigma@w)
        if V is not None and kappa:u-=kappa*np.sqrt(max(0.,w@V@w))
        return -u
    c=[
      {"type":"eq","fun":lambda x:x.sum()-1},
      {"type":"ineq","fun":lambda x:x[idx["국내주식"]]-.15},{"type":"ineq","fun":lambda x:.25-x[idx["국내주식"]]},
      {"type":"ineq","fun":lambda x:x[FX].sum()-.30},{"type":"ineq","fun":lambda x:.42-x[FX].sum()},
      {"type":"ineq","fun":lambda x:x[idx["국내 국채(초장기)"]]-.15},{"type":"ineq","fun":lambda x:.25-x[idx["국내 국채(초장기)"]]},
      {"type":"ineq","fun":lambda x:x[idx["글로벌 IG 크레딧"]]-.05},{"type":"ineq","fun":lambda x:.12-x[idx["글로벌 IG 크레딧"]]},
      {"type":"ineq","fun":lambda x:x[ALT].sum()-.12},{"type":"ineq","fun":lambda x:.18-x[ALT].sum()},
      {"type":"ineq","fun":lambda x:x[EQ].sum()-.52},{"type":"ineq","fun":lambda x:.60-x[EQ].sum()},
      {"type":"ineq","fun":lambda x:x[BD].sum()-.25},{"type":"ineq","fun":lambda x:.35-x[BD].sum()},
    ]
    if use_bl_limits:
        for j in range(n):
            c += [{"type":"ineq","fun":lambda x,j=j:ACTIVE_CAP-(x[j]-w_prior_total[j])},
                  {"type":"ineq","fun":lambda x,j=j:ACTIVE_CAP+(x[j]-w_prior_total[j])}]
        c.append({"type":"ineq","fun":lambda x:te_cap-np.sqrt(max(0.,(x[:n]-w_prior_total)@Sigma@(x[:n]-w_prior_total)))})
    r=minimize(obj,x0,method="SLSQP",bounds=b9,constraints=c,options={"ftol":1e-14,"maxiter":10000})
    if not r.success: raise RuntimeError(r.message)
    return r.x
cm=solve_cash(mu_total,use_bl_limits=False)
cb=solve_cash(mu_bl_total); cr=solve_cash(mu_bl_total,V_BL,KAPPA)
cashrows=[]
for nm,x,m_eval in [("Policy_MVO_cash_optimized",cm,mu_total),("BL-MVO_cash_optimized",cb,mu_bl_total),("Robust_BL_cash_optimized",cr,mu_bl_total)]:
    cashrows.append({"method":nm,**{assets[i]:x[i] for i in range(n)},"단기자금":x[-1],**metrics8(x[:n],m_eval,x[-1])})
pd.DataFrame(cashrows).to_csv(RESULTS/"step12_cash_sensitivity.csv",index=False,encoding="utf-8-sig")

meta={"delta":DELTA,"tau":TAU,"T_years_baseline":10.,"rf":RF,"gamma":GAMMA,
 "prior_weight":prior_weight,"cma_weight":cma_weight,"P":"identity","Omega":"Sigma/T",
 "T_interpretation":"Proxy assumption: W3 CMA is forward-looking, so the 120-month risk sample is converted to a 10-year equivalent solely to operationalize the assignment Omega formula.",
 "delta_vs_gamma":"delta=2.5 reverse-engineers the policy prior; gamma=4 is the portfolio utility risk-aversion coefficient. They play different roles.",
 "cash_baseline":"Policy cash fixed at official 0.1%; 0-2% optimization is reported as sensitivity.",
 "active_cap":ACTIVE_CAP,"TE_cap":TE_CAP,"robust_BL_kappa":KAPPA,
 "notes":["Step 11 Team target is not used as the BL prior.",
          "Baseline TE limit remains the assignment-required 1.0%; 0.95% is only an implementation-buffer sensitivity.",
          "BL-MVO and Robust BL are reported without manual weight edits."]}
(RESULTS/"step12_meta.json").write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding="utf-8")
print("Prior/CMA weights",prior_weight,cma_weight)
print("Robust BL",w_rbl)
