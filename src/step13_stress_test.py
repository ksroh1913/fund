# -*- coding: utf-8 -*-
"""STEP 13 - assignment stress tests after independent audit.

Baseline portfolios are read from Step 11/12 result files; no hard-coded Team
or BL weights remain.

Stress 1
--------
Baseline follows the assignment literally: CMA misses by 1 standard error in
the direction adverse to the Robust-BL active position.  Fixed portfolios are
evaluated directly under stressed CMA returns (not passed through BL again).
SE_i = sqrt(Omega_ii).  Absolute return, active return versus Official, and
utility are reported.

A separate *self-active-direction* sensitivity applies the same 1SE rule
against each candidate's own active weights.  This is not a replacement for
the assignment stress; it is a fairness diagnostic showing how each candidate
behaves when its own bets are wrong.

Stress 2
--------
Baseline policy-equity definition = domestic equity, DM, EM.  A sensitivity
also treats PE as an equity asset.  If a shocked correlation matrix is not PSD,
the correction is disclosed.

Stress 3
--------
Alternative vol x1.5 and correlation with DM/EM +0.20.  Both eigenvalue-clipping
and Higham nearest-correlation corrections are reported.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data"; RESULTS=ROOT/"results"; RESULTS.mkdir(exist_ok=True)
RF=.03; GAMMA=4.; CASH=.001

proxy=pd.read_csv(DATA/"proxy_mapping.csv")
assets=proxy["asset"].tolist();idx={a:i for i,a in enumerate(assets)}
mu=proxy.set_index("asset").loc[assets,"mu_cma"].astype(float).values
sigma=proxy.set_index("asset").loc[assets,"sigma_cma"].astype(float).values
Sigma=pd.read_csv(DATA/"covariance_cma.csv",index_col=0).loc[assets,assets].astype(float).values
corr=pd.read_csv(DATA/"correlation.csv",index_col=0).loc[assets,assets].astype(float).values
Omega=pd.read_csv(RESULTS/"step12_Omega.csv",index_col=0).loc[assets,assets].astype(float).values
se=np.sqrt(np.diag(Omega))

# Read portfolio weights from the immediately preceding steps.
table_c=pd.read_csv(RESULTS/"step11_tableC_detailed.csv",index_col=0)
alloc=pd.read_csv(RESULTS/"step12_allocations.csv").set_index("asset")
w_off=alloc.loc[assets,"official_mapped_total"].values
w_team=table_c.loc[assets,"Team scenario"].values
w_bl=alloc.loc[assets,"Robust_BL"].values
for w in (w_off,w_team,w_bl):
    if abs(w.sum()-.999)>1e-8: raise ValueError("Expected risky total 99.9%")

portfolios={"Official":w_off,"Team":w_team,"Robust_BL":w_bl}

def pm(w,m,C,cash=CASH):
    er=float(m@w+RF*cash);vol=float(np.sqrt(w@C@w));u=float(er-GAMMA/2*vol**2)
    return er,vol,u

def clip_corr(C,eps=1e-8):
    C=(C+C.T)/2;v,V=np.linalg.eigh(C);X=V@np.diag(np.maximum(v,eps))@V.T
    d=np.sqrt(np.diag(X));X=X/np.outer(d,d);np.fill_diagonal(X,1)
    return (X+X.T)/2

def higham_corr(C,it=2000,tol=1e-12):
    Y=(C+C.T)/2; dS=np.zeros_like(Y)
    for _ in range(it):
        R=Y-dS;v,V=np.linalg.eigh((R+R.T)/2);X=V@np.diag(np.maximum(v,0))@V.T
        dS=X-R;Y0=Y;Y=X.copy();np.fill_diagonal(Y,1)
        if np.max(np.abs(Y-Y0))<tol: break
    return (Y+Y.T)/2

def cov_from(C,s): return np.outer(s,s)*C

# ---------------- Stress 1: direct CMA miss, fixed weights ----------------
active_bl=w_bl-w_off
shock=np.where(active_bl>0,-se,np.where(active_bl<0,se,0.))
mu_s1=mu+shock
base_cma={k:pm(w,mu,Sigma) for k,w in portfolios.items()}
s1={k:pm(w,mu_s1,Sigma) for k,w in portfolios.items()}
s1rows=[]
for k,w in portfolios.items():
    er0,v0,u0=base_cma[k];er1,v1,u1=s1[k]
    off0=base_cma["Official"][0];off1=s1["Official"][0]
    s1rows.append({"portfolio":k,"baseline_expected_return":er0,"stress_expected_return":er1,
                   "expected_return_change":er1-er0,
                   "baseline_active_return_vs_official":er0-off0,
                   "stress_active_return_vs_official":er1-off1,
                   "active_return_change":(er1-off1)-(er0-off0),
                   "baseline_volatility":v0,"stress_volatility":v1,
                   "baseline_utility":u0,"stress_utility":u1,"utility_change":u1-u0})
pd.DataFrame(s1rows).to_csv(RESULTS/"step13_stress1_direct_cma.csv",index=False,encoding="utf-8-sig")
pd.DataFrame({"asset":assets,"CMA_total":mu,"Omega_SE":se,"BL_active":active_bl,
              "CMA_shock":shock,"stressed_CMA_total":mu_s1}
).to_csv(RESULTS/"step13_stress1_view_miss.csv",index=False,encoding="utf-8-sig")

# ---- Stress 1 supplemental sensitivity: shock each portfolio's own active direction ----
w_blend=.5*w_team+.5*w_bl
self_ports={"Team":w_team,"Robust_BL":w_bl,"Blend_50_RBL":w_blend}
self_rows=[]
for k,w in self_ports.items():
    own_active=w-w_off
    own_shock=np.where(own_active>0,-se,np.where(own_active<0,se,0.))
    m_self=mu+own_shock
    er0,v0,u0=pm(w,mu,Sigma); er1,v1,u1=pm(w,m_self,Sigma)
    off0=pm(w_off,mu,Sigma)[0]; off1=pm(w_off,m_self,Sigma)[0]
    self_rows.append({
        "portfolio":k,
        "baseline_expected_return":er0,
        "stress_expected_return":er1,
        "expected_return_change":er1-er0,
        "baseline_active_return_vs_official":er0-off0,
        "stress_active_return_vs_official":er1-off1,
        "active_return_change":(er1-off1)-(er0-off0),
        "baseline_volatility":v0,
        "stress_volatility":v1,
        "baseline_utility":u0,
        "stress_utility":u1,
        "utility_change":u1-u0,
    })
pd.DataFrame(self_rows).to_csv(
    RESULTS/"step13_stress1_self_active_sensitivity.csv",
    index=False,encoding="utf-8-sig"
)

# ---------------- Stress 2 ----------------
def stress2(equity_assets):
    eq=[idx[a] for a in equity_assets]
    m=mu.copy();m[eq]-=.02
    C=corr.copy()
    for z,i in enumerate(eq):
        for j in eq[z+1:]: C[i,j]=C[j,i]=min(.99,C[i,j]+.15)
    raw=float(np.linalg.eigvalsh(C).min());fixed=False
    if raw<1e-10: C=higham_corr(C);fixed=True
    cov=cov_from(C,sigma)
    return m,C,cov,raw,fixed
eq_base=["국내주식","글로벌 선진국 DM","글로벌 신흥국 EM"]
eq_pe=eq_base+["사모주식 PE/VC"]
m2,c2,C2,e2,fix2=stress2(eq_base)
m2p,c2p,C2p,e2p,fix2p=stress2(eq_pe)

# ---------------- Stress 3 ----------------
alt=["사모주식 PE/VC","실물 인프라","사모대출 PD"]; ge=["글로벌 선진국 DM","글로벌 신흥국 EM"]
c3=corr.copy(); s3=sigma.copy()
for a in alt:s3[idx[a]]*=1.5
for a in alt:
    for b in ge:
        i,j=idx[a],idx[b];c3[i,j]=c3[j,i]=min(.99,c3[i,j]+.20)
raw3=float(np.linalg.eigvalsh(c3).min())
c3clip=clip_corr(c3);c3hi=higham_corr(c3)
C3clip=cov_from(c3clip,s3); C3hi=cov_from(c3hi,s3)

def scenario_rows(label,m,C,baseline_mu=mu,baseline_C=Sigma):
    out=[]
    for k,w in portfolios.items():
        er,v,u=pm(w,m,C);er0,v0,u0=pm(w,baseline_mu,baseline_C)
        out.append({"scenario":label,"portfolio":k,"expected_return":er,"volatility":v,"utility":u,
                    "expected_return_change":er-er0,"volatility_change":v-v0,"utility_change":u-u0})
    return out
rows=[]
rows += scenario_rows("Stress2_policy_equity",m2,C2)
rows += scenario_rows("Stress2_include_PE_sensitivity",m2p,C2p)
rows += scenario_rows("Stress3_eigen_clip",mu,C3clip)
rows += scenario_rows("Stress3_Higham",mu,C3hi)
pd.DataFrame(rows).to_csv(RESULTS/"step13_stress_results.csv",index=False,encoding="utf-8-sig")

# T sensitivity: T changes posterior allocation and therefore BL active direction;
# direct Stress-1 still evaluates the candidate under CMA +/- sqrt(Sigma/T).
Ts=pd.read_csv(RESULTS/"step12_T_sensitivity.csv")
trows=[]
for T in sorted(Ts["T_years"].unique()):
    r=Ts[(Ts["T_years"]==T)&(Ts["method"]=="Robust BL")].iloc[0]
    wr=r[assets].astype(float).values
    seT=np.sqrt(np.diag(Sigma/float(T)))
    act=wr-w_off; sh=np.where(act>0,-seT,np.where(act<0,seT,0.))
    mT=mu+sh
    for name,w in [("Team",w_team),("Robust_BL",wr)]:
        er0=pm(w,mu,Sigma)[0];er1=pm(w,mT,Sigma)[0]
        off0=pm(w_off,mu,Sigma)[0];off1=pm(w_off,mT,Sigma)[0]
        trows.append({"T_years":T,"portfolio":name,"baseline_expected_return":er0,
                      "stress_expected_return":er1,"expected_return_change":er1-er0,
                      "baseline_active_return_vs_official":er0-off0,
                      "stress_active_return_vs_official":er1-off1,
                      "active_return_change":(er1-off1)-(er0-off0)})
pd.DataFrame(trows).to_csv(RESULTS/"step13_T_stress1_sensitivity.csv",index=False,encoding="utf-8-sig")

# Diagnostics and realized shock disclosure.
diag=[
 {"scenario":"Stress2_policy_equity","raw_min_corr_eigenvalue":e2,"correction":"Higham" if fix2 else "none",
  "final_min_corr_eigenvalue":float(np.linalg.eigvalsh(c2).min())},
 {"scenario":"Stress2_include_PE","raw_min_corr_eigenvalue":e2p,"correction":"Higham" if fix2p else "none",
  "final_min_corr_eigenvalue":float(np.linalg.eigvalsh(c2p).min())},
 {"scenario":"Stress3","raw_min_corr_eigenvalue":raw3,"correction":"eigenvalue clipping",
  "final_min_corr_eigenvalue":float(np.linalg.eigvalsh(c3clip).min())},
 {"scenario":"Stress3","raw_min_corr_eigenvalue":raw3,"correction":"Higham",
  "final_min_corr_eigenvalue":float(np.linalg.eigvalsh(c3hi).min())},
]
pd.DataFrame(diag).to_csv(RESULTS/"step13_covariance_diagnostics.csv",index=False,encoding="utf-8-sig")
pd.DataFrame({"asset":alt,
 "intended_DM_shock":[.20]*3,
 "realized_DM_shock_clip":[c3clip[idx[a],idx["글로벌 선진국 DM"]]-corr[idx[a],idx["글로벌 선진국 DM"]] for a in alt],
 "realized_DM_shock_Higham":[c3hi[idx[a],idx["글로벌 선진국 DM"]]-corr[idx[a],idx["글로벌 선진국 DM"]] for a in alt]}
).to_csv(RESULTS/"step13_stress3_realized_shocks.csv",index=False,encoding="utf-8-sig")

meta={
 "stress1":"Assignment baseline: direct CMA 1SE miss; SE=sqrt(diag(Omega)); adverse direction defined by Robust BL active sign; fixed weights; no posterior re-filtering.",
 "stress1_self_active_sensitivity":"Supplemental only: each candidate is shocked against its own active direction using the same SE=sqrt(diag(Omega)). This does not replace the assignment-defined Stress 1.",
 "stress2_baseline_equities":eq_base,
 "stress2_PE_sensitivity":"PE is mapped as an alternative in the baseline; a PE-as-equity sensitivity is separately reported.",
 "stress3":"Both eigenvalue clipping and Higham nearest-correlation corrections are reported because the raw shocked matrix is not PSD.",
 "evaluation_basis":"CMA total returns and common CMA Sigma, so all portfolios are compared on one basis.",
}
(RESULTS/"step13_meta.json").write_text(json.dumps(meta,ensure_ascii=False,indent=2),encoding="utf-8")
print(pd.DataFrame(s1rows))
print(pd.DataFrame(self_rows))
print(pd.DataFrame(rows))
