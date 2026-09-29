# External data for empirical reproduction

No market observations, prepared panels, cached arrays, or database exports are
included in LFPM. Supply a directory outside the package and the repository.
The methods also accept your own complete return/factor matrices through
`lfpm_empirical_test()`.

## Sources and required files

| Input | Source | Local filename under `data_dir` |
|---|---|---|
| China monthly stock returns, capitalization, market type | [CSMAR](https://data.csmar.com/) | `China/CSMAR/TRDNEW_Mnth_*_CSMAR.csv` (one or more files) |
| China value-weighted FF3, market P9706 | [CSMAR](https://data.csmar.com/) | `China/CSMAR/STK_MKT_THRFACMONTH_Three_factor_2020501_202604.csv` |
| Optional China one-year RF audit series | [CSMAR](https://data.csmar.com/) | `China/CSMAR/BND_Exchange_月度化无风险利率.csv` |
| U.S. adjusted daily OHLC archive | [Kaggle: S&P 500 Stocks, 25 Years](https://www.kaggle.com/datasets/darkmatternet/s-and-p-500-stocks-25-years-of-data-updated-daily) | `America/sp500_stocks.csv` |
| U.S. point-in-time constituent adjusted daily archive | [Kaggle: S&P500 index stocks](https://www.kaggle.com/datasets/joebeachcapital/s-and-p500-index-stocks-daily-updated) | `America/sp500_pit_data.csv` |
| Historical S&P 500 constituents | [fja05680/sp500](https://github.com/fja05680/sp500) | `America/sp500_historical_components.csv` |
| U.S. monthly FF3 and RF | [Kenneth French Data Library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html) | `America/F-F_Research_Data_Factors.CSV` |
| U.S. monthly FF5 and RF | [Kenneth French Data Library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html) | `America/F-F_Research_Data_5_Factors_2x3.csv` |

The unusual `2020501` spelling in the China factor filename is retained from
the study's export name. Rename an equivalent export to that name after checking
its contents. CSMAR is an access-controlled database: a reader needs the relevant
institutional subscription or other authorized access. Its [official access
guide](https://www.csmar.com/contents/9/861.html) describes institutional permissions
and purchases. Kaggle may require an account. Download and use each source under
its own terms; the package does not automate authentication or redistribute data.

The two U.S. price archives and component history are third-party sources, not
CRSP. They lack CRSP permanent identifiers and delisting returns. Files whose
values are raw rather than adjusted prices are not interchangeable with the
specified inputs. The first price archive's `close` field was produced with
`yfinance` auto-adjustment in the study's source vintage. Verify that convention
when acquiring a replacement vintage.

## Field contract

* CSMAR stock CSVs: `Stkcd`, `Trdmnt` (`YYYY-MM`), `Mretwd`, `Msmvttl3`,
  `Markettype`. `Mretwd` is the **decimal** monthly return including cash-dividend
  reinvestment. Stock codes are padded to six digits. The loader rejects
  conflicting duplicate stock-month rows and keeps `Markettype` 1 or 4.
* CSMAR FF3: `MarkettypeID`, `TradingMonth`, `RiskPremium2`, `SMB2`, `HML2`.
  Keep P9706, with factors already in decimal monthly return units.
* Optional RF audit file: `Clsdt`, `Nrrmtdt`; the original parser takes the last
  observation each month and divides by 100. This file does not define the main
  reconstructed three-month RF and is not required for main results.
* First U.S. daily file: `date`, `close`, `symbol`; use adjusted closing prices.
* Point-in-time daily file: skip its first header row and read the first eight
  fields as `Date,Ticker,Adj Close,Close,High,Low,Open,Volume`, matching the archived
  export. `Adj Close` must include the source's split/dividend adjustments.
* Membership file: `date`, `tickers`, with comma-separated tickers in a quoted
  field, unique dates, and rows sorted chronologically. Ticker dots are normalized
  to hyphens in the prices and membership sources.
* French files: preserve the official monthly CSV header beginning `,Mkt-RF`.
  Monthly rows use six-digit `YYYYMM`. FF3 columns are `Mkt-RF,SMB,HML,RF`; FF5
  adds `RMW,CMA` before `RF`. The loader divides percent returns by 100 once.

## Selection and timing

China uses January 2021--December 2025, the top 1,000 firms ranked by **December
2020** total market capitalization, then complete nonconstant monthly returns.
The paper has 980 assets and 60 months. Supply at least December 2020--December
2025 stock records, and factors for January 2021--December 2025.

The main Chinese RF is reconstructed from public **three-month deposit** annual
benchmarks using month-end effective rates and `(1+r)^(1/12)-1`. Its annual rate
is 1.10% throughout the main window. This is an explicitly documented benchmark
reconstruction, not an exact CSMAR internal RF export; CSMAR timing and rounding
could not be independently verified. Sources recorded in the original research
include [Fujian Rural Credit history](https://www.fjnx.com.cn/xpls/dtjrxx/rmbckllb/index.html),
[postal-bank disclosure](https://pdf.dfcfw.com/pdf/H2_AN201906281336805446_1.pdf), and
[the 2015 rate-change notice](https://www.beijing.gov.cn/fwcj/jiage/zffw1/65dc557711a82834a863972e.html).

U.S. uses January 2012--December 2019. Daily adjusted prices are sampled at the
last available trading day each month before computing simple monthly returns.
Supply prices beginning no later than December 2011. Membership is intersected
across snapshots through the window, including the last snapshot at/before its
start, and complete nonconstant returns are required. The paper has 288 assets
and 96 months. This is a long-tenure constituent sample with survivorship-related
selection, not a CRSP universe. Both U.S. price files are required.

The U.S. FF3 model includes factor-loading interactions from 2017-12; FF5 uses
2015-07. Each has one full-window intercept. FF5 excess returns are reconstructed
from the stock returns with **FF5's own RF column**.

## Reproducibility boundary

Public sources are revised and the Kaggle archives update daily. A current
download can differ from the study's input snapshot even within historical
months. Exact reproduction needs the same source vintages, adjustment rules,
and fields; links alone do not guarantee identical market observations. The
paper's FF5 export header identified the 202608 CRSP database vintage.

Each run writes input SHA-256 hashes, source filenames, software versions,
random seeds, draw counts, and method conventions to the requested output
directory. A paper workflow stops if its selected panel dimensions differ from
980 by 60 or 288 by 96. This catches sample mismatch but does not prove the input
values match the paper. No original market observations are bundled to bypass
the data-access or vintage requirements.
