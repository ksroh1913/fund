# -*- coding: utf-8 -*-
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'data'
RESULTS = ROOT / 'results'
RESULTS.mkdir(exist_ok=True)

RF=0.03
DELTA=2.5
TAU=0.025
GAMMA=4.0
CASH=0.001
T_YEARS=10.0

proxy=pd.read_csv(DATA/'proxy_mapping.csv')
assets=proxy['asset'].tolist()
mu=proxy.set_index('asset').loc[assets,'mu_cma'].astype(float).values
sigma=proxy.set_index('asset').loc[assets,'sigma_cma'].astype(float).values
Sigma_df=pd.read_csv(DATA/'covariance_cma.csv',index_col=0)
Sigma=Sigma_df.loc[assets,assets].astype(float).values
corr_df=pd.read_csv(DATA/'correlation.csv',index_col=0)
corr=corr_df.loc[assets,assets].astype(float).values
Smu_df=pd.read_csv(RESULTS/'step8_S_mu.csv',index_col=0)
Smu=Smu_df.loc[assets,assets].astype(float).values
se=np.sqrt(np.diag(Smu))
idx={a:i for i,a in enumerate(assets)}

w_off=np.array([
    0.208,
    0.356*0.92,
    0.356*0.08,
    0.218,
    0.074,
    0.143/3,
    0.143/3,
    0.143/3,
])
assert abs(w_off.sum()-0.999)<1e-12

w_team=np.array([0.220,0.275,0.065,0.205,0.085,0.030,0.065,0.054])
assert abs(w_team.sum()-0.999)<1e-12

w_bl=np.array([
0.1916614561381548,
0.30087353641282205,
0.027465007449023242,
0.24800000000000005,
0.09246857679294879,
0.01766666666666664,
0.043198089873717926,
0.07766666666666668,
])
assert abs(w_bl.sum()-0.999)<1e-12

w_prior_risky=w_off/w_off.sum()
pi=DELTA*Sigma@w_prior_risky
Q=mu-RF
prior_weight=1/(1+TAU*T_YEARS)
cma_weight=TAU*T_YEARS/(1+TAU*T_YEARS)
mu_bl_excess=prior_weight*pi+cma_weight*Q
mu_bl=mu_bl_excess+RF

portfolios={'Official':w_off,'Team':w_team,'Robust_BL':w_bl}

def nearest_corr_psd(C, eps=1e-8):
    C=(C+C.T)/2
    vals, vecs=np.linalg.eigh(C)
    clipped=np.maximum(vals,eps)
    Cp=vecs@np.diag(clipped)@vecs.T
    d=np.sqrt(np.diag(Cp))
    Cp=Cp/np.outer(d,d)
    np.fill_diagonal(Cp,1.0)
    return (Cp+Cp.T)/2

def cov_from_corr(C, sig):
    return np.outer(sig,sig)*C

def pm(w, mu_eval, cov_eval):
    er=float(mu_eval@w + RF*CASH)
    vol=float(np.sqrt(w@cov_eval@w))
    util=float(er - GAMMA/2*vol**2)
    return er,vol,util

base={name:pm(w,mu_bl,Sigma) for name,w in portfolios.items()}

# Stress 1: CMA miss by 1 bootstrap SE adverse to final BL active sign.
active_bl=w_bl-w_off
shock_cma=np.where(active_bl>0,-se,np.where(active_bl<0,se,0.0))
Q_s1=Q+shock_cma
mu_s1=prior_weight*pi+cma_weight*Q_s1+RF
Sigma_s1=Sigma.copy()

# Stress 2: public-equity expected return -2pp and pairwise equity correlation +0.15.
equity=['국내주식','글로벌 선진국 DM','글로벌 신흥국 EM']
mu_s2=mu_bl.copy()
for a in equity:
    mu_s2[idx[a]]-=0.02
corr_s2=corr.copy()
for i,a in enumerate(equity):
    for b in equity[i+1:]:
        ia,ib=idx[a],idx[b]
        corr_s2[ia,ib]=corr_s2[ib,ia]=min(0.99,corr_s2[ia,ib]+0.15)
min_eig_s2_raw=float(np.linalg.eigvalsh(corr_s2).min())
psd_fix_s2=False
if min_eig_s2_raw < 1e-10:
    corr_s2=nearest_corr_psd(corr_s2)
    psd_fix_s2=True
Sigma_s2=cov_from_corr(corr_s2,sigma)

# Stress 3: alternatives vol +50%; correlation with DM/EM +0.20.
alt=['사모주식 PE/VC','실물 인프라','사모대출 PD']
global_eq=['글로벌 선진국 DM','글로벌 신흥국 EM']
mu_s3=mu_bl.copy()
sig_s3=sigma.copy()
for a in alt:
    sig_s3[idx[a]]*=1.5
corr_s3=corr.copy()
for a in alt:
    for b in global_eq:
        ia,ib=idx[a],idx[b]
        corr_s3[ia,ib]=corr_s3[ib,ia]=min(0.99,corr_s3[ia,ib]+0.20)
min_eig_s3_raw=float(np.linalg.eigvalsh(corr_s3).min())
psd_fix_s3=False
if min_eig_s3_raw < 1e-10:
    corr_s3=nearest_corr_psd(corr_s3)
    psd_fix_s3=True
Sigma_s3=cov_from_corr(corr_s3,sig_s3)

scenarios={
    'Baseline':(mu_bl,Sigma),
    'Stress1_CMA_1SE':(mu_s1,Sigma_s1),
    'Stress2_Equity':(mu_s2,Sigma_s2),
    'Stress3_Alternatives':(mu_s3,Sigma_s3),
}

rows=[]
for scen,(m,c) in scenarios.items():
    for pname,w in portfolios.items():
        er,vol,util=pm(w,m,c)
        ber,bvol,butil=base[pname]
        rows.append({
            'scenario':scen,
            'portfolio':pname,
            'expected_return':er,
            'volatility':vol,
            'utility':util,
            'expected_return_change':er-ber,
            'volatility_change':vol-bvol,
            'utility_change':util-butil,
        })
res=pd.DataFrame(rows)
res.to_csv(RESULTS/'step13_stress_results.csv',index=False,encoding='utf-8-sig')

s1_detail=pd.DataFrame({
    'asset':assets,
    'CMA_total':mu,
    'CMA_SE':se,
    'BL_active':active_bl,
    'CMA_shock':shock_cma,
    'stressed_CMA_total':mu+shock_cma,
    'baseline_BL_total':mu_bl,
    'stressed_BL_total':mu_s1,
    'BL_total_change':mu_s1-mu_bl,
})
s1_detail.to_csv(RESULTS/'step13_stress1_view_miss.csv',index=False,encoding='utf-8-sig')

diag=pd.DataFrame([
    ['Stress2', 'raw_min_corr_eigenvalue', min_eig_s2_raw],
    ['Stress2', 'PSD_fix_applied', psd_fix_s2],
    ['Stress2', 'final_min_corr_eigenvalue', float(np.linalg.eigvalsh(corr_s2).min())],
    ['Stress2', 'final_min_cov_eigenvalue', float(np.linalg.eigvalsh(Sigma_s2).min())],
    ['Stress3', 'raw_min_corr_eigenvalue', min_eig_s3_raw],
    ['Stress3', 'PSD_fix_applied', psd_fix_s3],
    ['Stress3', 'final_min_corr_eigenvalue', float(np.linalg.eigvalsh(corr_s3).min())],
    ['Stress3', 'final_min_cov_eigenvalue', float(np.linalg.eigvalsh(Sigma_s3).min())],
],columns=['scenario','check','value'])
diag.to_csv(RESULTS/'step13_covariance_diagnostics.csv',index=False,encoding='utf-8-sig')

pd.DataFrame(corr_s2,index=assets,columns=assets).to_csv(RESULTS/'step13_corr_stress2.csv',encoding='utf-8-sig')
pd.DataFrame(Sigma_s2,index=assets,columns=assets).to_csv(RESULTS/'step13_cov_stress2.csv',encoding='utf-8-sig')
pd.DataFrame(corr_s3,index=assets,columns=assets).to_csv(RESULTS/'step13_corr_stress3.csv',encoding='utf-8-sig')
pd.DataFrame(Sigma_s3,index=assets,columns=assets).to_csv(RESULTS/'step13_cov_stress3.csv',encoding='utf-8-sig')

rel=[]
for scen,(m,c) in scenarios.items():
    er_off=pm(w_off,m,c)[0]
    for pname,w in {'Team':w_team,'Robust_BL':w_bl}.items():
        er=pm(w,m,c)[0]
        rel.append({'scenario':scen,'portfolio':pname,'expected_return_vs_official':er-er_off})
rel_df=pd.DataFrame(rel)
rel_df.to_csv(RESULTS/'step13_relative_to_official.csv',index=False,encoding='utf-8-sig')

baseline_alpha = rel_df[rel_df['scenario']=='Baseline'].set_index('portfolio')['expected_return_vs_official']
stress1_alpha = rel_df[rel_df['scenario']=='Stress1_CMA_1SE'].set_index('portfolio')['expected_return_vs_official']
active_loss_rows=[]
for pname in ['Team','Robust_BL']:
    active_loss_rows.append({
        'portfolio': pname,
        'baseline_alpha_vs_official': baseline_alpha[pname],
        'stress1_alpha_vs_official': stress1_alpha[pname],
        'stress1_active_loss': stress1_alpha[pname] - baseline_alpha[pname],
    })
pd.DataFrame(active_loss_rows).to_csv(
    RESULTS/'step13_stress1_active_loss.csv', index=False, encoding='utf-8-sig'
)

meta={
    'baseline_expected_return':'Policy BL posterior total return',
    'baseline_covariance':'Team CMA covariance Sigma',
    'stress1':{
        'rule':'CMA total return +/- 1 bootstrap SE adverse to Robust BL active sign; then propagated through Policy BL posterior (20% CMA weight)',
        'CMA_SE_source':'Step 8 bootstrap S_mu diagonal',
    },
    'stress2':{
        'equity_assets':equity,
        'expected_return_shock_pp':-2.0,
        'pairwise_equity_correlation_shock':0.15,
        'PSD_fix_applied':psd_fix_s2,
    },
    'stress3':{
        'alternative_assets':alt,
        'global_equity_assets':global_eq,
        'alternative_vol_multiplier':1.5,
        'alt_global_equity_correlation_shock':0.20,
        'PSD_fix_applied':psd_fix_s3,
    },
    'gamma_for_utility':GAMMA,
    'cash_weight':CASH,
}
(RESULTS/'step13_meta.json').write_text(
    json.dumps(meta,ensure_ascii=False,indent=2),encoding='utf-8'
)

print(res.to_string(index=False))
print('\nRelative to official:')
print(rel_df.to_string(index=False))
print('\nCov diag:')
print(diag.to_string(index=False))
