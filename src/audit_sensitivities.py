# -*- coding: utf-8 -*-
"""Supplemental sensitivity checks required by the independent audit.

These are not alternative baselines. They document how conclusions change when:
- Box widths use bootstrap standard errors,
- Ellipsoid uses Sigma/T instead of proxy-bootstrap S_mu,
- the Team's 10% detailed caps are relaxed.

Baseline assignment results remain in Steps 5-13.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import chi2

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data"; RESULTS=ROOT/"results"; RESULTS.mkdir(exist_ok=True)
RF=.03;G=4.;K=float(np.sqrt(chi2.ppf(.5,8)))
proxy=pd.read_csv(DATA/"proxy_mapping.csv")
rets=pd.read_csv(DATA/"monthly_returns_krw.csv",index_col=0)
assets=proxy["asset"].tolist();idx={a:i for i,a in enumerate(assets)};n=len(assets)
mu=proxy.set_index("asset").loc[assets,"mu_cma"].astype(float).values
sig=proxy.set_index("asset").loc[assets,"sigma_cma"].astype(float).values
S=pd.read_csv(DATA/"covariance_cma.csv",index_col=0).loc[assets,assets].astype(float).values
R=rets[assets].astype(float).values;T=len(R)
w_off=np.array([.208,.356*.92,.356*.08,.218,.074,.143/3,.143/3,.143/3]);w_slv=w_off/w_off.sum()
iKR,iDM,iEM,iGB,iIG,iPE,iIN,iPD=range(8);EQ=[iKR,iDM,iEM];BD=[iGB,iIG];ALT=[iPE,iIN,iPD];FX=[iDM,iEM]
base_bounds=[(0,1),(0,1),(0,.10),(0,1),(0,1),(0,.10),(0,.10),(0,.10)]
relaxed=[(0,1)]*n
def constraints():
 return [
 {"type":"eq","fun":lambda w:w.sum()-1},
 {"type":"ineq","fun":lambda w:w[iKR]-.15},{"type":"ineq","fun":lambda w:.25-w[iKR]},
 {"type":"ineq","fun":lambda w:w[FX].sum()-.30},{"type":"ineq","fun":lambda w:.42-w[FX].sum()},
 {"type":"ineq","fun":lambda w:w[iGB]-.15},{"type":"ineq","fun":lambda w:.25-w[iGB]},
 {"type":"ineq","fun":lambda w:w[iIG]-.05},{"type":"ineq","fun":lambda w:.12-w[iIG]},
 {"type":"ineq","fun":lambda w:w[ALT].sum()-.12},{"type":"ineq","fun":lambda w:.18-w[ALT].sum()},
 {"type":"ineq","fun":lambda w:w[EQ].sum()-.52},{"type":"ineq","fun":lambda w:.60-w[EQ].sum()},
 {"type":"ineq","fun":lambda w:w[BD].sum()-.25},{"type":"ineq","fun":lambda w:.35-w[BD].sum()}]
def solve(m,C,bounds=base_bounds,U=None,kappa=0,x0=None):
 if x0 is None:x0=w_slv.copy()
 def obj(w):
  u=m@w-G/2*(w@C@w)
  if U is not None:u-=kappa*np.sqrt(max(0,w@U@w))
  return -u
 r=minimize(obj,x0,method="SLSQP",bounds=bounds,constraints=constraints(),options={"ftol":1e-12,"maxiter":10000})
 if not r.success:raise RuntimeError(r.message)
 return r.x

# Bootstrap S_mu used by baseline Ellipsoid.
rng=np.random.default_rng(60080)
boot=np.array([R[rng.integers(0,T,T)].mean(0)*12 for _ in range(5000)])
Smu=np.cov(boot,rowvar=False,ddof=1);se=np.sqrt(np.diag(Smu))
D=np.array([.01,.01,.01,.005,.005,.015,.015,.015])

# Manual Ledoit-Wolf correlation.
X=R-R.mean(0);Sm=X.T@X/T;tr=np.trace(Sm)/n;F=tr*np.eye(n);d2=np.sum((Sm-F)**2)
b2=min(sum(np.sum((np.outer(x,x)-Sm)**2) for x in X)/T**2,d2);delta=b2/d2
Sl=delta*F+(1-delta)*Sm;sd=np.sqrt(np.diag(Sl));rho_lw=Sl/np.outer(sd,sd);S_lw=np.outer(sig,sig)*rho_lw

rows=[]
cases=[
 ("MVO",mu,S,None,0),
 ("Ledoit-Wolf",mu,S_lw,None,0),
 ("Box default",mu-D,S,None,0),
 ("Box bootstrap-SE",mu-se,S,None,0),
 ("Ellipsoid bootstrap",mu,S,Smu,K),
 ("Ellipsoid Sigma/T",mu,S,S/(T/12),K),
]
for name,m,C,U,k in cases:
 for cap_label,bounds in [("10% detailed caps",base_bounds),("relaxed detailed caps",relaxed)]:
  w=solve(m,C,bounds,U,k)
  rows.append({"case":name,"detail_bounds":cap_label,**{assets[i]:w[i] for i in range(n)}})
pd.DataFrame(rows).to_csv(RESULTS/"audit_method_sensitivities.csv",index=False,encoding="utf-8-sig")

# Pure Michaud with relaxed detailed caps (baseline pure Michaud is Step 9).
full=R.mean(0);rng=np.random.default_rng(60080);W=[]
for _ in range(300):
 rb=R[rng.integers(0,T,T)];mb=mu+12*(rb.mean(0)-full)
 corr=np.corrcoef(rb,rowvar=False);C=np.outer(sig,sig)*corr
 W.append(solve(mb,C,relaxed))
W=np.asarray(W)
pd.DataFrame({"asset":assets,"mean_relaxed":W.mean(0),"p05_relaxed":np.quantile(W,.05,0),
              "p50_relaxed":np.quantile(W,.50,0),"p95_relaxed":np.quantile(W,.95,0)}
).to_csv(RESULTS/"audit_michaud_relaxed_caps.csv",index=False,encoding="utf-8-sig")

notes=pd.DataFrame([
 {"item":"Box bootstrap-SE","interpretation":"Sensitivity only. Much wider uncertainty widths push the solution to policy corners; this demonstrates that calibration size, not only set geometry, drives Box-vs-Ellipsoid differences."},
 {"item":"Ellipsoid Sigma/T","interpretation":"Sensitivity only. Baseline remains the assignment-recommended proxy bootstrap S_mu; Sigma/T shows dependence on uncertainty-scale choice."},
 {"item":"10% detail caps","interpretation":"Exact 10% levels are Team assumptions, not official NPS limits. Relaxed-cap outputs show how much Part-2 results depend on them."},
])
notes.to_csv(RESULTS/"audit_sensitivity_notes.csv",index=False,encoding="utf-8-sig")
print(pd.DataFrame(rows))
