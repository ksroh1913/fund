# -*- coding: utf-8 -*-
"""STEP 2 - Team 2 NPS 8-asset correlation/covariance builder.

- Keeps Team 2 CMA mu and sigma.
- Estimates rho from 120 complete monthly observations, 2016-09 to 2026-08.
- Converts USD assets to unhedged KRW returns.
- Writes CSV inputs used by later optimization steps.
"""

from pathlib import Path
import numpy as np
import pandas as pd
import yfinance as yf

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

START = "2016-08-01"
END = "2026-09-01"   # 2026-08 month-end is the last complete month
MIN_T = 120

PROXIES = {
    "국내주식": ("069500.KS", "KODEX 200", True),
    "글로벌 선진국 DM": ("URTH", "iShares MSCI World ETF", False),
    "글로벌 신흥국 EM": ("EEM", "iShares MSCI Emerging Markets ETF", False),
    "국내 국채(초장기)": ("152380.KS", "KODEX 국채선물10년", True),
    "글로벌 IG 크레딧": ("LQD", "iShares iBoxx $ Investment Grade Corporate Bond ETF", False),
    "사모주식 PE/VC": ("PSP", "Invesco Global Listed Private Equity ETF", False),
    "실물 인프라": ("IGF", "iShares Global Infrastructure ETF", False),
    "사모대출 PD": ("BIZD", "VanEck BDC Income ETF", False),
}
FX_TICKER = "KRW=X"

CMA_MU = pd.Series({
    "국내주식": 0.062,
    "글로벌 선진국 DM": 0.056,
    "글로벌 신흥국 EM": 0.084,
    "국내 국채(초장기)": 0.036,
    "글로벌 IG 크레딧": 0.053,
    "사모주식 PE/VC": 0.068,
    "실물 인프라": 0.068,
    "사모대출 PD": 0.078,
})

CMA_SIGMA = pd.Series({
    "국내주식": 0.195,
    "글로벌 선진국 DM": 0.162,
    "글로벌 신흥국 EM": 0.210,
    "국내 국채(초장기)": 0.065,
    "글로벌 IG 크레딧": 0.078,
    "사모주식 PE/VC": 0.220,
    "실물 인프라": 0.110,
    "사모대출 PD": 0.095,
})


def monthly_level(ticker: str) -> pd.Series:
    d = yf.download(
        ticker, start=START, end=END,
        auto_adjust=False, progress=False,
        actions=False, threads=False
    )
    if d.empty:
        raise RuntimeError(f"{ticker}: no price data")

    if isinstance(d.columns, pd.MultiIndex):
        lvl0 = d.columns.get_level_values(0)
        col = "Adj Close" if "Adj Close" in lvl0 else "Close"
        s = d[col].iloc[:, 0]
    else:
        col = "Adj Close" if "Adj Close" in d.columns else "Close"
        s = d[col]

    s = pd.to_numeric(s, errors="coerce").dropna()
    s.index = pd.to_datetime(s.index).tz_localize(None)
    return s.resample("ME").last()


levels = {}
for asset, (ticker, _, _) in PROXIES.items():
    levels[asset] = monthly_level(ticker)
fx = monthly_level(FX_TICKER).rename("USDKRW")

levels = pd.concat(levels, axis=1)
krw = levels.copy()

for asset, (_, _, domestic_krw) in PROXIES.items():
    if not domestic_krw:
        krw[asset] = levels[asset] * fx

krw = krw.dropna(how="any")

# Need 121 prices for 120 monthly returns.
if len(krw) < MIN_T + 1:
    raise RuntimeError(f"only {len(krw)} common month-end prices")

krw = krw.tail(MIN_T + 1)
rets = krw.pct_change().dropna()

expected_start = pd.Timestamp("2016-09-30")
expected_end = pd.Timestamp("2026-08-31")
if len(rets) != 120 or rets.index[0] != expected_start or rets.index[-1] != expected_end:
    raise RuntimeError(
        f"unexpected sample: {rets.index[0]} ~ {rets.index[-1]}, T={len(rets)}"
    )

corr = rets.corr()
sig = CMA_SIGMA.loc[corr.columns].values
cov_cma = pd.DataFrame(
    np.outer(sig, sig) * corr.values,
    index=corr.index, columns=corr.columns
)

proxy_rows = []
for asset, (ticker, proxy, domestic_krw) in PROXIES.items():
    proxy_rows.append({
        "asset": asset,
        "ticker": ticker,
        "proxy": proxy,
        "currency_basis": "KRW" if domestic_krw else "KRW unhedged",
        "mu_cma": CMA_MU[asset],
        "sigma_cma": CMA_SIGMA[asset],
    })
proxy = pd.DataFrame(proxy_rows)

diagnostics = pd.DataFrame([
    ["T_monthly", len(rets)],
    ["start_month", rets.index.min().strftime("%Y-%m-%d")],
    ["end_month", rets.index.max().strftime("%Y-%m-%d")],
    ["corr_symmetry_max_abs_error", float(np.abs(corr.values-corr.values.T).max())],
    ["corr_diag_max_abs_error", float(np.abs(np.diag(corr.values)-1).max())],
    ["corr_min_eigenvalue", float(np.linalg.eigvalsh(corr.values).min())],
    ["corr_condition_number", float(np.linalg.cond(corr.values))],
    ["cov_min_eigenvalue", float(np.linalg.eigvalsh(cov_cma.values).min())],
    ["cov_condition_number", float(np.linalg.cond(cov_cma.values))],
], columns=["check", "value"])

proxy.to_csv(DATA / "proxy_mapping.csv", index=False, encoding="utf-8-sig")
rets.to_csv(DATA / "monthly_returns_krw.csv", encoding="utf-8-sig")
corr.to_csv(DATA / "correlation.csv", encoding="utf-8-sig")
cov_cma.to_csv(DATA / "covariance_cma.csv", encoding="utf-8-sig")
diagnostics.to_csv(DATA / "diagnostics.csv", index=False, encoding="utf-8-sig")

print(diagnostics.to_string(index=False))
print("\nSaved to", DATA)
