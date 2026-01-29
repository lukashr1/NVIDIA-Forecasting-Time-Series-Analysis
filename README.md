# NVIDIA Stock Price Forecasting - Time Series Analysis

This project performs comprehensive time series analysis and forecasting of NVIDIA (NVDA) stock prices using ARIMA and GARCH models in R.

## Overview

The analysis includes:
- **ARIMA modeling** for mean return forecasting
- **GARCH modeling** for volatility forecasting
- Train/Validation/Test split evaluation
- Expanding window forecasting methodology
- Comprehensive model diagnostics and comparison

## Project Structure

```
.
├── nvidia.R          # Main R script with complete analysis
├── nvda_data.rds     # Cached NVIDIA stock data
├── nvidia_4165.pdf   # Analysis report/documentation
└── README.md         # This file
```

## Requirements

### R Packages
The following R packages are required:

```r
library(quantmod)      # Financial data and charting
library(forecast)      # Time series forecasting
library(xts)          # Time series objects
library(tseries)      # Time series analysis
library(urca)         # Unit root and cointegration tests
library(xtable)       # Table formatting
library(ggplot2)      # Data visualization
library(dplyr)        # Data manipulation
library(rugarch)      # GARCH modeling
library(zoo)          # Time series infrastructure
library(FinTS)        # Financial time series tools
```

### Installation

Install required packages in R:

```r
install.packages(c("quantmod", "forecast", "xts", "tseries", "urca", 
                   "xtable", "ggplot2", "dplyr", "rugarch", "zoo", "FinTS"))
```

## Data

The script downloads NVIDIA stock data from Yahoo Finance:
- **Ticker**: NVDA
- **Date Range**: 2020-03-03 to 2025-03-03
- **Data Source**: Yahoo Finance via `quantmod` package

Data is automatically cached in `nvda_data.rds` to avoid repeated downloads.

## Analysis Workflow

### 1. Data Preparation
- Downloads/loads NVIDIA stock data
- Extracts adjusted closing prices
- Splits data into:
  - Training set (60%)
  - Validation set (20%)
  - Test set (20%)

### 2. ARIMA Model Analysis

#### Stationarity Testing
- Augmented Dickey-Fuller (ADF) tests on price levels and returns
- Log transformation for variance stabilization
- First differencing to achieve stationarity

#### Model Selection
Tests multiple ARIMA specifications:
- ARIMA(1,0,0) - Best model by BIC
- ARIMA(0,0,1)
- ARIMA(0,0,2)
- ARIMA(2,0,0)
- ARIMA(2,0,1) - Selected by auto.arima

#### Model Diagnostics
- Ljung-Box tests for residual autocorrelation
- ACF/PACF plots
- Normality tests (Q-Q plots, histograms)
- Residual analysis

#### Forecasting Strategy
- **Validation Set**: 60-day expanding window forecasting
- **Test Set**: One-step-ahead forecasts with ARIMA(1,0,0)
- Hit ratio analysis (directional accuracy)
- Binomial test for forecast significance

### 3. GARCH Model Analysis

#### ARCH Effects Detection
- Tests for volatility clustering in ARIMA residuals
- Ljung-Box tests on squared residuals
- ARCH LM tests

#### GARCH Model Selection
Tests GARCH specifications:
- GARCH(1,1) - Selected model (best by BIC)
- GARCH(1,2)
- GARCH(2,1)
- GARCH(2,2)

Uses Student's t-distribution for error terms.

#### Volatility Forecasting
- Expanding window approach on test set
- One-step-ahead conditional volatility forecasts
- Comparison with realized volatility (30-day rolling SD)

#### Model Diagnostics
- Standardized residual tests
- ARCH effect removal verification
- ACF/PACF of standardized residuals

## Running the Analysis

1. **Set Working Directory** (if using RStudio):
   ```r
   setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
   ```

2. **Run the Script**:
   ```r
   source("nvidia.R")
   ```

The script will:
- Download data (or load from cache)
- Perform all analyses
- Generate diagnostic plots
- Output forecast accuracy metrics

## Key Results

### ARIMA Model Performance
- **Best Model**: ARIMA(1,0,0) based on BIC
- **Test Set Metrics**: RMSE, MAE, MAPE
- **Hit Ratio**: Percentage of correct directional forecasts
- **Statistical Significance**: Binomial test results

### GARCH Model Performance
- **Best Model**: GARCH(1,1) with Student's t-distribution
- **Volatility Forecasting**: RMSE and MAE on test set
- **ARCH Effects**: Successfully removed from residuals

## Visualizations

The script generates multiple plots:
- Stock price time series with train/val/test splits
- Log returns and differenced series
- ACF and PACF plots
- Forecast vs actual comparisons
- Squared residuals (ARCH effects)
- Hit ratio scatter plots
- Volatility comparison plots

## Model Evaluation Metrics

- **RMSE**: Root Mean Squared Error
- **MAE**: Mean Absolute Error  
- **MAPE**: Mean Absolute Percentage Error
- **Hit Ratio**: Directional accuracy of forecasts
- **Information Criteria**: AIC, BIC for model selection

## Notes

- The analysis uses expanding window forecasting for realistic out-of-sample evaluation
- Models are refitted at each forecast step to incorporate new information
- The two-step approach (ARIMA + GARCH) separately models mean and volatility
- Results depend on the data period and market conditions

## References

- **ARIMA**: Autoregressive Integrated Moving Average models
- **GARCH**: Generalized Autoregressive Conditional Heteroskedasticity models
- **ADF Test**: Augmented Dickey-Fuller test for unit roots
- **Ljung-Box Test**: Test for autocorrelation in residuals

## License

This project is for educational and research purposes.

## Author

Time series analysis of NVIDIA stock prices using statistical forecasting methods.
