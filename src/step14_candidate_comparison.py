# -*- coding: utf-8 -*-
"""STEP 14 candidate comparison before the final IC decision.

This file intentionally does not hard-code the final decision first.  It
constructs a common-basis comparison of:
- Team target
- Robust BL (assignment TE cap 1.0%)
- Robust BL with 0.95% TE implementation buffer
- 25/50/75% convex blends between Team and baseline Robust BL

Stress 1 uses the baseline Robust-BL active direction, as required by the
assignment's BL-active adverse-view-miss definition.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/"data";RESULTS=ROOT/"results";RESULTS.mkdir(exist_ok=True)
RF=.03;G=4.;CASH=.001;AUM=1865.6;ACTUAL_KR=.291
proxy=pd.read_csv(DATA/"proxy_mapping.csv");assets=proxy["asset"].tolist();idx={a:i for i,a in enumerate(assets)}
mu=proxy.set_index("asset").loc[assets,"mu_cma"].astype(float).values
sig=proxy.set_index("asset").loc[assets,"sigma_cma"].astype(float).values
S=pd.read_csv(DATA/"covariance_cma.csv",index_col=0).loc[assets,assets].astype(float).values
rho=pd.read_csv(DATA/"correlation.csv",index_col=0).loc[assets,assets].astype(float).values
alloc=pd.read_csv(RESULTS/"step12_allocations.csv").set_index("asset")
tc=pd.read_csv(RESULTS/"step11_tableC_detailed.csv",index_col=0)
w_off=alloc.loc[assets,"official_mapped_total"].values
w_team=tc.loc[assets,"Team scenario"].values
w_rbl=alloc.loc[assets,"Robust_BL"].values
te95=pd.read_csv(RESULTS/"step12_TE095_sensitivity.csv")
w_rbl95=te95.loc[te95["method"]=="Robust_BL_TE095",assets].iloc[0].astype(float).values
Omega=pd.read_csv(RESULTS/"step12_Omega.csv",index_col=0).loc[assets,assets].astype(float).values

def met(w,m=mu,C=S):
    er=float(m@w+RF*CASH);v=float(np.sqrt(w@C@w));a=w-w_off
    return er,v,float(er-G/2*v*v),float(np.sqrt(a@S@a)),float(.5*np.abs(a).sum())

# Stress 1
se=np.sqrt(np.diag(Omega));active=w_rbl-w_off;shock=np.where(active>0,-se,np.where(active<0,se,0.));m1=mu+shock
# Stress 2
EQ=[idx[a] for a in ["국내주식","글로벌 선진국 DM","글로벌 신흥국 EM"]]
m2=mu.copy();m2[EQ]-=.02;c2=rho.copy()
for z,i in enumerate(EQ):
    for j in EQ[z+1:]:c2[i,j]=c2[j,i]=min(.99,c2[i,j]+.15)
C2=np.outer(sig,sig)*c2
# Stress 3 Higham
ALT=[idx[a] for a in ["사모주식 PE/VC","실물 인프라","사모대출 PD"]]
FX=[idx[a] for a in ["글로벌 선진국 DM","글로벌 신흥국 EM"]]
c3=rho.copy();s3=sig.copy();s3[ALT]*=1.5
for i in ALT:
    for j in FX:c3[i,j]=c3[j,i]=min(.99,c3[i,j]+.20)
def higham(C,it=2000,tol=1e-12):
    Y=(C+C.T)/2;dS=np.zeros_like(C)
    for _ in range(it):
        R=Y-dS;v,V=np.linalg.eigh((R+R.T)/2);X=V@np.diag(np.maximum(v,0))@V.T
        dS=X-R;Y0=Y;Y=X.copy();np.fill_diagonal(Y,1)
        if np.max(np.abs(Y-Y0))<tol:break
    return (Y+Y.T)/2
C3=np.outer(s3,s3)*higham(c3)

candidates={"Team":w_team,"Robust BL TE1.00":w_rbl,"Robust BL TE0.95":w_rbl95}
for alpha in (.25,.50,.75):
    candidates[f"Blend {int(alpha*100)}% RBL"]=(1-alpha)*w_team+alpha*w_rbl

off_base=met(w_off);off_s1=met(w_off,m1,S)
rows=[]
for name,w in candidates.items():
    er,v,u,te,to=met(w); er1,v1,u1,_,_=met(w,m1,S); er2,v2,u2,_,_=met(w,m2,C2); er3,v3,u3,_,_=met(w,mu,C3)
    active0=er-off_base[0];active1=er1-off_s1[0]
    rows.append({"candidate":name,**{assets[i]:w[i] for i in range(len(assets))},
                 "expected_return":er,"volatility_common_Sigma":v,"sharpe_common_Sigma":(er-RF)/v,"utility":u,
                 "TE":te,"turnover":to,
                 "stress1_active_loss":active1-active0,
                 "stress2_utility_change":u2-u,"stress3_Higham_utility_change":u3-u,
                 "domestic_equity_full_sell_trn":(ACTUAL_KR-w[idx["국내주식"]])*AUM,
                 "domestic_equity_incremental_vs_official_trn":(.208-w[idx["국내주식"]])*AUM})
df=pd.DataFrame(rows)
df.to_csv(RESULTS/"step14_candidate_comparison.csv",index=False,encoding="utf-8-sig")
print(df[["candidate","expected_return","volatility_common_Sigma","utility","TE","turnover","stress1_active_loss","stress2_utility_change","stress3_Higham_utility_change","domestic_equity_incremental_vs_official_trn"]])
