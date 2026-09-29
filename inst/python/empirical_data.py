"""Raw-file readers and paper sample selection. Inputs stay outside the package."""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
DATA_DIR = None
OUT = None
CACHE = None

def dump(value,path):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,default=str),encoding='utf-8')

def configure(data_dir,output_dir):
    global DATA_DIR, OUT, CACHE
    DATA_DIR=Path(data_dir).resolve()
    OUT=Path(output_dir).resolve()
    OUT.mkdir(parents=True,exist_ok=True)
    CACHE=OUT/'data_cache'

def get_data():
    (OUT/'panels').mkdir(parents=True,exist_ok=True)
    cache=OUT/'source_panels.pkl'  # Rebuild from caller-supplied files; never reuse stale panels.
    chunks=[]; audit=[]
    for f in sorted((DATA_DIR/'China/CSMAR').glob('TRDNEW_Mnth_*_CSMAR.csv')):
        x=pd.read_csv(f,dtype={'Stkcd':str})
        audit.append(dict(file=f.name,rows=len(x),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
        chunks.append(x)
    if not chunks:
        raise FileNotFoundError('No China/CSMAR/TRDNEW_Mnth_*_CSMAR.csv files found under data_dir')
    x=pd.concat(chunks,ignore_index=True)
    x['Stkcd']=x.Stkcd.str.zfill(6)
    dup=x.duplicated(['Stkcd','Trdmnt'],keep=False)
    if dup.any():
        conflicts=x[dup].groupby(['Stkcd','Trdmnt'])[['Mretwd','Msmvttl3','Markettype']].nunique().gt(1).any(axis=1)
        if conflicts.any(): raise ValueError('Conflicting stock-month duplicates')
    x=x.drop_duplicates(['Stkcd','Trdmnt'])
    # Exact match to P9706: Shanghai/Shenzhen main-board A shares only.
    x=x[x.Markettype.isin([1,4])].copy()
    x['month']=pd.PeriodIndex(x.Trdmnt,freq='M')
    first=x.groupby('Stkcd').month.min()
    r=x.pivot(index='month',columns='Stkcd',values='Mretwd').sort_index()
    cap=x.pivot(index='month',columns='Stkcd',values='Msmvttl3').sort_index()
    factfile=DATA_DIR/'China/CSMAR/STK_MKT_THRFACMONTH_Three_factor_2020501_202604.csv'
    f=pd.read_csv(factfile);f=f[f.MarkettypeID=='P9706'].copy()
    f.index=pd.PeriodIndex(f.TradingMonth,freq='M');f=f[['RiskPremium2','SMB2','HML2']].sort_index()
    assert not f.index.duplicated().any()
    assert (r.stack()>=-1).all()
    # Published benchmark annual rates; source A through 2015-05, source B afterwards.
    rate_entries=[('2004-10-29',1.71),('2006-08-19',1.80),('2007-03-18',1.98),
      ('2007-05-19',2.07),('2007-07-21',2.34),('2007-08-22',2.61),('2007-09-15',2.88),
      ('2007-12-21',3.33),('2008-10-09',3.15),('2008-10-30',2.88),('2008-11-27',1.98),
      ('2008-12-23',1.71),('2010-10-20',1.91),('2010-12-26',2.25),('2011-02-09',2.60),
      ('2011-04-06',2.85),('2011-07-07',3.10),('2012-06-08',2.85),('2012-07-06',2.60),
      ('2014-11-22',2.35),('2015-03-01',2.10),('2015-05-11',1.85),('2015-06-28',1.60),
      ('2015-08-26',1.35),('2015-10-24',1.10)]
    rates=pd.Series(dict(rate_entries),name='annual_3m_deposit_pct');rates.index=pd.to_datetime(rates.index)
    monthly=rates.reindex(r.index.to_timestamp(how='end'),method='ffill')
    rf=pd.Series(np.expm1(np.log1p(monthly.to_numpy()/100)/12),index=r.index,name='rf_3m_reconstructed')
    # Exact source monthlyization is unavailable: retain explicit approximation and sensitivities.
    rf1_path=DATA_DIR/'China/CSMAR/BND_Exchange_月度化无风险利率.csv'
    if rf1_path.exists():
        rf1=pd.read_csv(rf1_path)
        rf1['month']=pd.to_datetime(rf1.Clsdt).dt.to_period('M')
        rf1=rf1.groupby('month').Nrrmtdt.last()/100
    else:
        rf1=pd.Series(np.nan,index=r.index)
    pd.DataFrame({'effective_date':rates.index,'annual_rate_pct':rates.values}).to_csv(OUT/'rf_rate_history.csv',index=False)
    pd.DataFrame({'rf_3m_reconstructed':rf,'rf_1y_supplied':rf1,'rf_3m_simple':pd.Series(monthly.to_numpy()/1200,index=r.index)}).to_csv(OUT/'rf_monthly_audit.csv')
    audit.append(dict(note='RF reconstructed from public three-month benchmark rates; effective monthly conversion and month-end convention are analyst choices, not verified CSMAR export.',
      source_a='https://www.fjnx.com.cn/xpls/dtjrxx/rmbckllb/index.html',
      source_b='https://pdf.dfcfw.com/pdf/H2_AN201906281336805446_1.pdf',
      source_c='https://www.beijing.gov.cn/fwcj/jiage/zffw1/65dc557711a82834a863972e.html',
      factor_market='P9706',asset_markettypes=[1,4],duplicate_rows=int(dup.sum()),asset_count=r.shape[1]))
    dump(audit,OUT/'data_audit.json')
    obj=dict(returns=r,cap=cap,factors=f,rf=rf,rf1=rf1,first=first)
    return obj

def panel(data,s):
    idx=pd.period_range(s['start'],s['end'],freq='M'); start=idx[0]
    r=data['returns'].reindex(idx)
    cap=data['cap'].loc[start-1].dropna(); cap=cap[cap>0].sort_values(ascending=False,kind='mergesort')
    eligible=cap.index.intersection(r.columns,sort=False)
    chosen=list(eligible if s['pool']=='full' else eligible[:int(s['pool'][3:])])
    mat=r.reindex(columns=chosen)
    good=mat.notna().all() & mat.var().gt(1e-14)
    kept=list(mat.columns[good]); f=data['factors'].reindex(idx)
    assert not f.isna().any().any()
    rr=mat[kept];y=rr.sub(data['rf'].reindex(idx),axis=0)
    attr=dict(N=len(kept),N_candidate=len(chosen),retention=len(kept)/len(chosen),
      missing_rate=float(mat.isna().to_numpy().mean()),zero_return_rate=float((rr==0).to_numpy().mean()),
      extreme_return_rate=float((rr.abs()>.5).to_numpy().mean()),cap_date=str(start-1),
      market='China main-board A',rf='3m benchmark reconstructed')
    return y,f,attr,mat.notna().sum()

def read_french_factors() -> pd.DataFrame:
    path = DATA_DIR / "America" / "F-F_Research_Data_Factors.CSV"
    lines = path.read_text(encoding="latin1").splitlines()
    header = next(i for i, line in enumerate(lines) if line.startswith(",Mkt-RF"))
    rows = []
    for line in lines[header + 1 :]:
        if not line.strip():
            break
        cells = [x.strip() for x in line.split(",")]
        if len(cells) >= 5 and len(cells[0]) == 6 and cells[0].isdigit():
            rows.append(cells[:5])
    f = pd.DataFrame(rows, columns=["yyyymm", "MKT", "SMB", "HML", "RF"])
    f.index = pd.PeriodIndex(f.pop("yyyymm"), freq="M")
    return f.apply(pd.to_numeric, errors="coerce") / 100.0

def read_us_monthly() -> dict:
    cache = CACHE / "us_adjusted_monthly_combined.pkl"
    daily_path = DATA_DIR / "America" / "sp500_stocks.csv"
    d = pd.read_csv(daily_path, usecols=["date", "close", "symbol"])
    d["date"] = pd.to_datetime(d["date"], errors="coerce")
    d = d.dropna(subset=["date", "close", "symbol"])
    d["symbol"] = d["symbol"].astype(str).str.replace(".", "-", regex=False)
    # The source was produced with yfinance auto_adjust=True.  Its OHLC fields
    # are adjusted; the supplied raw-price workbook is retained only for audit.
    d["month"] = d["date"].dt.to_period("M")
    last = d.sort_values("date").groupby(["month", "symbol"], sort=False).tail(1)
    current_prices = last.pivot(index="month", columns="symbol", values="close").sort_index()

    # Supplement the current-constituent file with point-in-time constituents
    # that later left the index.  The downloaded CSV has a malformed first
    # header row; only the first eight documented fields are read, and the
    # malformed 1972 observation is outside every analysis window.
    pit_path = DATA_DIR / "America" / "sp500_pit_data.csv"
    pit_parts = []
    if pit_path.exists():
        pit_names = ["Date", "Ticker", "Adj Close", "Close", "High", "Low", "Open", "Volume"]
        for chunk in pd.read_csv(pit_path, skiprows=1, names=pit_names, usecols=range(8), chunksize=500_000):
            chunk["Date"] = pd.to_datetime(chunk["Date"], errors="coerce")
            chunk["Adj Close"] = pd.to_numeric(chunk["Adj Close"], errors="coerce")
            chunk = chunk.dropna(subset=["Date", "Ticker", "Adj Close"])
            chunk["Ticker"] = chunk["Ticker"].str.replace(".", "-", regex=False)
            chunk["month"] = chunk["Date"].dt.to_period("M")
            pit_parts.append(chunk[["month", "Ticker", "Date", "Adj Close"]])
        pit_daily = pd.concat(pit_parts, ignore_index=True)
        pit_last = pit_daily.sort_values("Date").groupby(["month", "Ticker"], sort=False).tail(1)
        pit_prices = pit_last.pivot(index="month", columns="Ticker", values="Adj Close").sort_index()
        prices = current_prices.combine_first(pit_prices)
    else:
        raise FileNotFoundError("Paper reproduction requires America/sp500_pit_data.csv")
    returns = prices.pct_change(fill_method=None)

    hist = pd.read_csv(DATA_DIR / "America" / "sp500_historical_components.csv")
    hist["date"] = pd.to_datetime(hist["date"])
    if not hist["date"].is_monotonic_increasing or hist["date"].duplicated().any():
        raise ValueError('Membership snapshot dates must be unique and sorted chronologically')
    hist["members"] = hist["tickers"].str.replace(".", "-", regex=False).str.split(",").map(set)
    obj = {"prices": prices, "returns": returns, "membership": hist, "factors": read_french_factors()}
    return obj

def continuous_members(hist: pd.DataFrame, start: str, end: str) -> set[str]:
    first = pd.Period(start, "M").start_time
    last = pd.Period(end, "M").end_time.normalize()
    prior = hist.loc[hist["date"] <= first]
    during = hist.loc[(hist["date"] >= first) & (hist["date"] <= last)]
    sets = list(during["members"])
    if not prior.empty:
        sets.insert(0, prior.iloc[-1]["members"])
    if not sets:
        raise ValueError(f"No S&P membership snapshots for {start}--{end}")
    return set.intersection(*sets)

def us_panel(data: dict, start: str, end: str, universe: str = "continuous"):
    months = pd.period_range(start, end, freq="M")
    ret = data["returns"].reindex(months)
    complete = set(ret.columns[ret.notna().all(axis=0) & ret.var(axis=0).gt(1e-14)])
    if universe == "continuous":
        eligible = continuous_members(data["membership"], start, end)
        assets = sorted(complete & eligible)
    elif universe == "current_complete":
        eligible = set(ret.columns)
        assets = sorted(complete)
    else:
        raise ValueError(universe)
    f = data["factors"].reindex(months)
    if f.isna().any().any():
        raise ValueError(f"Missing French factors in {start}--{end}")
    y = ret[assets].sub(f["RF"], axis=0)
    attr = {
        "market": "United States",
        "start": start,
        "end": end,
        "T": len(months),
        "N": len(assets),
        "N_membership_intersection": len(eligible),
        "N_complete_prices": len(complete),
        "universe": universe,
        "N_over_T": len(assets) / len(months),
    }
    return y, f[["MKT", "SMB", "HML"]], attr, assets

def french(path, five=False):
    lines = path.read_text(encoding="latin1").splitlines()
    header = next(i for i, line in enumerate(lines) if line.startswith(",Mkt-RF"))
    names = ["yyyymm", "MKT", "SMB", "HML"] + (["RMW", "CMA"] if five else []) + ["RF"]
    rows = []
    for line in lines[header + 1:]:
        if not line.strip():
            break
        cells = [cell.strip() for cell in line.split(",")]
        if len(cells) >= len(names) and len(cells[0]) == 6 and cells[0].isdigit():
            rows.append(cells[:len(names)])
    result = pd.DataFrame(rows, columns=names)
    result.index = pd.PeriodIndex(result.pop("yyyymm"), freq="M")
    return result.apply(pd.to_numeric, errors="raise") / 100.0
