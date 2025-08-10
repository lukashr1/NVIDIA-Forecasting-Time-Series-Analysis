setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Packages
library(quantmod)
library(forecast)
library(xts)
library(tseries)
library(urca)
library(xtable)
library(ggplot2)
library(dplyr)
library(rugarch)
library(zoo)
library(FinTS)

# Downloading the data 

# File path
file_path <- "nvda_data.rds"

# Checking if the file already exists
if (file.exists(file_path)) {
  # Loading saved data
  nvda_df <- readRDS(file_path)
  print("Loaded saved NVDA data.")
} else {
  # Downloading data and saving it
  nvda_df <- getSymbols("NVDA", src = "yahoo", auto.assign = FALSE, from="2020-03-03", to="2025-03-03")
  saveRDS(nvda_df, file = file_path)
  print("Downloaded and saved NVDA data.")
}

### Analysis ARIMA model ###

# Extract adjusted closing price as a time series
nvda_ts <- xts(nvda_df$NVDA.Adjusted, order.by = index(nvda_df))
plot(nvda_ts, main = "NVIDIA Adjusted Stock Price")

# Convert xts to data frame for ggplot
nvda_df <- data.frame(
  date  = index(nvda_ts),
  price = as.numeric(nvda_ts)
)

# Determine split indices for training, validation, and test sets
total_obs   <- nrow(nvda_df)
train_end   <- floor(0.6 * total_obs)
val_end     <- floor(0.8 * total_obs)

# Extract split dates for plotting
train_end_date <- nvda_df$date[train_end]
val_end_date   <- nvda_df$date[val_end]

# Create data frames for each set
train_df <- nvda_df[1:train_end, ]
val_df   <- nvda_df[(train_end + 1):val_end, ]
test_df  <- nvda_df[(val_end + 1):total_obs, ]

# Create xts objects for modeling 
train_ts <- xts(train_df$price, order.by = train_df$date)
val_ts   <- xts(val_df$price, order.by = val_df$date)
test_ts  <- xts(test_df$price, order.by = test_df$date)

# Plotting the full series with train/val/test split markers
ggplot(nvda_df, aes(x = date, y = price)) +
  geom_line(color = "black") +
  geom_vline(xintercept = as.numeric(train_end_date), color = "blue", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = as.numeric(val_end_date), color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = train_end_date, y = max(nvda_df$price), label = "Val/Train Split", hjust = -0.2, color = "blue") +
  annotate("text", x = val_end_date, y = max(nvda_df$price), label = "Val/Test Split", hjust = -0.2, color = "red") +
  labs(
    title = "Train, Validation, and Test Split of NVIDIA Prices",
    x = "Date",
    y = "Price"
  ) +
  scale_x_date(date_labels = "%Y", date_breaks = "1 year") +
  theme_minimal()

#---------Work on training set/ARIMA-----------

# Plotting training set
plot(train_ts, main = "Training Set: NVIDIA Stock Prices", ylab = "Price", xlab = "Date")

# Dickey-Fuller test to confirm non-stationarity, using urca package, using AIC according to Stock&Watson
ur01_AIC <- ur.df(train_ts, type = "none", lags = 20, selectlags = c("AIC"))

# Value of test-statistic is: 0.16, comparing with critical values from Dickey-Fuller,
# we find that there is atleast one unit root, making the time series non-stationary
summary(ur01_AIC)

# Taking log of the closing price, stabilizing variance
train_log <- log(train_ts)  

# Looks non-stationary, check with Dickey-Fuller test
plot(train_log, main="Log of NVIDIA Stock Price")

# Dickey-Fuller test on log transformed time series
ur02_AIC <- ur.df(train_log, type = "none", lags = 20, selectlags = c("AIC"))

# Value of test-statistic is: 1.18, still atleast one unit root present
summary(ur02_AIC)

# Trying to eliminate the unit root by doing first differencing, also eliminating NAs
train_returns <- na.omit(diff(train_log))

# Now exhibits stationary properties, lets check with ADF
plot(train_returns, main = "Training Set Returns")

# Now with type = none, we don´t observe a clear trend in the timeseries
ur03_AIC <- ur.df(train_returns, type="none", lags = 20, selectlags = c("AIC"))

# Value of test-statistic: -6.8 < DF critical values, reject null hypothesis of existing unit root
summary(ur03_AIC)

# Plot ACF and PACF for NVIDIA log returns

# ACF yields one significant lag 1, and cuts off after, suggesting MA(1) process
acf(train_returns, main="ACF of NVIDIA Returns", lag.max = 15)

# PACF yields one significant lag 1, and cuts off after, suggesting AR(1) process
pacf(train_returns, main="PACF of NVIDIA Returns", lag.max = 15)

# Looking at the ACF and PACF, it suggests an ARIMA(1,0,1) is a decent fit

# Fitting ARIMA(1,0,1)
arima_101 <- Arima(train_returns, order=c(1,0,1))

# ARMA(1,0,1)
summary(arima_101)

# Lets create a loop walking us through different ARIMA models, 
# letting the loop do the heavy work of picking our models

# First we define the range of our AR(p) and MA(q) terms, with no differencing 
p_values <- 0:8 # AR terms
q_values <- 0:8 # MA terms

# Creating a empty data frame to store the results
results <- data.frame(Model = character(), 
                      BIC = numeric(),
                      stringsAsFactors = FALSE
                      )

# Looping over all combinations of (p,0,q)
for (p in p_values) {
    for (q in q_values) {
      
      # Fit ARIMA model
      model <- Arima(train_returns, order = c(p, 0, q)) # 0 to reflect no differencing
      # Store BIC
      results <- rbind(results, data.frame(
        Model = paste0("ARIMA(", p, ",0,", q, ")"),
        BIC = BIC(model)
      ))
    }
}

# Sorting by BIC (best models first)
results_BIC <- results[order(results$BIC), ]

# Printing top models with the least AIC/BIC
print(results_BIC) # Suggests ARIMA(1,0,0), -6419.279

# ARIMA(1,0,0), best model according to BIC
arima_100 <- Arima(train_returns, order = c(1,0,0))
summary(arima_100)

# Lets test with auto.arima
auto_model <- auto.arima(train_returns)

# Auto.arima chooses ARIMA(2,0,1) as the best model
summary(auto_model)

# We will test 5 different ARIMA models on forecasting:
arima_100 <- Arima(train_returns, order = c(1,0,0)) 
arima_001 <- Arima(train_returns, order = c(0,0,1)) 
arima_002 <- Arima(train_returns, order = c(0,0,2)) 
arima_200 <- Arima(train_returns, order = c(2,0,0))
arima_201 <- Arima(train_returns, order = c(2,0,1))


summary(arima_100)
summary(arima_200)
summary(arima_001)
summary(arima_002)
summary(arima_201)

# Table for comparison of ICs

# Extracting BIC values for each model
model_names <- c("ARIMA(1,0,0)", "ARIMA(2,0,0)", "ARIMA(0,0,1)", "ARIMA(0,0,2)", "ARIMA(2,0,1")

bic_values <- c(
  BIC(arima_100),
  BIC(arima_200),
  BIC(arima_001),
  BIC(arima_002),
  BIC(arima_201)
)

# Creating a data frame for comparison on ICs
comparison_table <- data.frame(
  Model = model_names,
  BIC = bic_values
)

# Printing the table
print(comparison_table)

### DIAGNOSTICS ###

# Autocorrelation tests #

# Plotting the residuals from all models
plot(resid(arima_100), main="Residuals of ARIMA(1,0,0)") 
plot(resid(arima_200), main="Residuals of ARIMA(2,0,0)") 
plot(resid(arima_001), main="Residuals of ARIMA(0,0,1)")
plot(resid(arima_002), main="Residuals of ARIMA(0,0,2)") 
plot(resid(arima_201), main="Residuals of ARIMA(2,0,1)") 

# Ljung-Box test to check for autocorrelation in residuals, H0: white noise residuals - no autocorrelation
Box.test(resid(arima_100), type="Ljung-Box", lag = 30) 
Box.test(resid(arima_200), type="Ljung-Box", lag = 30) 
Box.test(resid(arima_001), type="Ljung-Box", lag = 30) # All models well specified
Box.test(resid(arima_002), type="Ljung-Box", lag = 30) 
Box.test(resid(arima_201), type="Ljung-Box", lag = 30) 

# ACF and PACF for all models
acf(resid(arima_100), main = "ACF of Residuals, ARIMA(1,0,0)")
pacf(resid(arima_100), main = "PACF of Residuals, ARIMA(1,0,0)")

acf(resid(arima_200), main = "ACF of Residuals, ARIMA(2,0,0)")
pacf(resid(arima_200), main = "PACF of Residuals, ARIMA(2,0,0)")

acf(resid(arima_001), main = "ACF of Residuals, ARIMA(0,0,1)")
pacf(resid(arima_001), main = "PACF of Residuals, ARIMA(0,0,1)")

acf(resid(arima_002), main = "ACF of Residuals, ARIMA(0,0,2)")
pacf(resid(arima_002), main = "PACF of Residuals, ARIMA(0,0,2)")

acf(resid(arima_201), main = "ACF of Residuals, ARIMA(2,0,1)")
pacf(resid(arima_201), main = "PACF of Residuals, ARIMA(2,0,1)")

# Histograms of residuals
hist(resid(arima_100), main = "Histogram of Residuals (ARIMA(1,0,0))", breaks = 30)
hist(resid(arima_200), main = "Histogram of Residuals (ARIMA(2,0,0))", breaks = 30)
hist(resid(arima_101), main = "Histogram of Residuals (ARIMA(1,0,1))", breaks = 30)
hist(resid(arima_001), main = "Histogram of Residuals (ARIMA(0,0,1))", breaks = 30)
hist(resid(arima_002), main = "Histogram of Residuals (ARIMA(0,0,2))", breaks = 30)
hist(resid(arima_201), main = "Histogram of Residuals (ARIMA(0,0,2))", breaks = 30)

# Q-Q plots of residuals
qqnorm(resid(arima_100), main = "Q-Q Plot of Residuals ARIMA(1,0,0)"); qqline(resid(arima_100), col = "red")
qqnorm(resid(arima_200), main = "Q-Q Plot of Residuals ARIMA(2,0,0)"); qqline(resid(arima_200), col = "red")
qqnorm(resid(arima_001), main = "Q-Q Plot of Residuals ARIMA(0,0,1)"); qqline(resid(arima_001), col = "red")
qqnorm(resid(arima_002), main = "Q-Q Plot of Residuals ARIMA(0,0,2)"); qqline(resid(arima_002), col = "red")
qqnorm(resid(arima_201), main = "Q-Q Plot of Residuals ARIMA(0,0,2)"); qqline(resid(arima_201), col = "red")


#---------------Validation set/ARIMA-----------------------------

# FORECASTING: Expanding-window ARIMA forecasts on validation set

# Compute log returns
full_returns_log <- log(nvda_ts)
full_returns <- na.omit(diff(full_returns_log))  # Log returns

# Setup window size
window_size <- 60

# Forecasting range
total_obs <- length(full_returns)
train_end <- floor(0.6 * total_obs)
val_end   <- floor(0.8 * total_obs)

# Validation starts immediately after training ends
val_start <- train_end + 1

# Validation ends at val_end, inclusive
validation_length <- val_end - train_end

# Create storage
forecasts_100 <- numeric(validation_length)
forecasts_001 <- numeric(validation_length)
forecasts_002 <- numeric(validation_length)
forecasts_200 <- numeric(validation_length)
forecasts_201 <- numeric(validation_length)

# Expanding-window loop: starting 60 days before validation start, then growing
for (i in 1:validation_length) {
  # First point in validation is train_end + 1, so:
  end_idx <- train_end + i - 1
  start_idx <- end_idx - window_size + 1  # ensures window starts 60 days back
  window_data <- full_returns[start_idx:end_idx]
  
  # Forecast with each ARIMA model
  model_100 <- Arima(window_data, order = c(1, 0, 0), method = "ML")
  forecasts_100[i] <- forecast(model_100, h = 1)$mean
  
  model_001 <- Arima(window_data, order = c(0, 0, 1), method = "ML")
  forecasts_001[i] <- forecast(model_001, h = 1)$mean
  
  model_002 <- Arima(window_data, order = c(0, 0, 2), method = "ML")
  forecasts_002[i] <- forecast(model_002, h = 1)$mean
  
  model_200 <- Arima(window_data, order = c(2, 0, 0), method = "ML")
  forecasts_200[i] <- forecast(model_200, h = 1)$mean
  
  model_201 <- Arima(window_data, order = c(2, 0, 1), method = "ML")
  forecasts_201[i] <- forecast(model_201, h = 1)$mean
}

# Actual values to compare
actual_values <- full_returns[(train_end + 1):(train_end + validation_length)]

# Evaluate forecast accuracy
accuracy_100 <- accuracy(forecasts_100, actual_values)
accuracy_001 <- accuracy(forecasts_001, actual_values)
accuracy_002 <- accuracy(forecasts_002, actual_values)
accuracy_200 <- accuracy(forecasts_200, actual_values)
accuracy_201 <- accuracy(forecasts_201, actual_values)

# Print results
print(accuracy_100)
print(accuracy_001)
print(accuracy_002)
print(accuracy_200)
print(accuracy_201)

## Visualizing ARIMA forecasts 

# Creating a data frame for plotting forecasts on validation set
forecast_df <- data.frame(
  Time     = index(full_returns)[(train_end + 1):(val_end)],
  Actual   = actual_values,
  ARIMA100 = forecasts_100,
  ARIMA001 = forecasts_001,
  ARIMA002 = forecasts_002,
  ARIMA200 = forecasts_200,
  ARIMA201 = forecasts_201
)

# Column names
colnames(forecast_df) <- c("Time", "Actual", 
                           "ARIMA(1,0,0)", "ARIMA(0,0,1)", "ARIMA(0,0,2)", 
                           "ARIMA(2,0,0)", "ARIMA(2,0,1)")

# Plot forecasts vs actual on validation period
ggplot(forecast_df, aes(x = Time)) +
  geom_line(aes(y = Actual,          color = "Actual Returns"), linewidth = 0.4) +
  geom_line(aes(y = `ARIMA(1,0,0)`,  color = "ARIMA(1,0,0) Forecast"), linewidth = 0.4) +
  geom_line(aes(y = `ARIMA(0,0,1)`,  color = "ARIMA(0,0,1) Forecast"), linewidth = 0.4) +
  geom_line(aes(y = `ARIMA(0,0,2)`,  color = "ARIMA(0,0,2) Forecast"), linewidth = 0.4) +
  geom_line(aes(y = `ARIMA(2,0,0)`,  color = "ARIMA(2,0,0) Forecast"), linewidth = 0.4) +
  geom_line(aes(y = `ARIMA(2,0,1)`,  color = "ARIMA(2,0,1) Forecast"), linewidth = 0.4) +
  labs(
    title = "Validation Set: Actual Returns vs ARIMA Forecasts",
    x = "Time", y = "Returns"
  ) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  scale_color_manual(values = c(
    "Actual Returns"           = "black",
    "ARIMA(1,0,0) Forecast"    = "red",
    "ARIMA(0,0,1) Forecast"    = "blue",
    "ARIMA(0,0,2) Forecast"    = "green",
    "ARIMA(2,0,0) Forecast"    = "purple",
    "ARIMA(2,0,1) Forecast"    = "orange"
  )) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

#----------------------Test set/ARIMA-------------------------

### FORECASTING: Expanding-window ARIMA(1,0,0) on test set ###

# Compute returns
full_returns_log <- log(nvda_ts)
full_returns <- na.omit(diff(full_returns_log))

# Set minimum window size
window_size <- 60

# Define range for forecasting over test set
returns_length <- length(full_returns)
val_end_returns <- floor(0.8 * returns_length)
test_length <- returns_length - val_end_returns 

# Create storage vector
forecasts_test_100 <- numeric(test_length)

# Expanding window forecast loop
for (i in 1:test_length) {
  end_idx <- val_end_returns + i - 1
  start_idx <- end_idx - window_size + 1  # Start 60 obs before first test point
  window_data <- full_returns[start_idx:end_idx]
  
  model_100 <- Arima(window_data, order = c(1, 0, 0), method = "ML")
  forecasts_test_100[i] <- forecast(model_100, h = 1)$mean
}

# Actual test values
actual_test_values <- full_returns[(val_end_returns + 1):returns_length]

# Evaluate forecast accuracy
accuracy_test_100 <- accuracy(forecasts_test_100, actual_test_values)
print(accuracy_test_100)

# Calculate Hit Ratio
hit_vector <- sign(forecasts_test_100) == sign(actual_test_values)
hit_ratio <- mean(hit_vector, na.rm = TRUE)

# Print Hit Ratio
cat("Hit Ratio (ARIMA(1,0,0) on Test Set):", round(hit_ratio * 100, 2), "%\n")

# Create time index for test set
test_index <- index(full_returns)[(val_end_returns + 1):returns_length]

hit_df <- data.frame(
  Date = test_index,
  Hit = hit_vector
)

# Plot of HIT ratio
colnames(hit_df) <- c("Date", "Hit")

ggplot(hit_df, aes(x = Date, y = as.numeric(Hit))) +
  geom_point(color = "steelblue", alpha = 0.5) +
  scale_y_continuous(breaks = c(0, 1), labels = c("Miss", "Hit")) +
  labs(
    title = "Hit Ratio of ARIMA(1,0,0) Forecasts on Test Set",
    x = "Date", y = "Forecast Direction Accuracy"
  ) +
  theme_minimal()

# Binomial test: is hit ratio significantly greater than 0.5?
hit_successes <- sum(hit_vector, na.rm = TRUE)
n_hits <- sum(!is.na(hit_vector))
binom_result <- binom.test(hit_successes, n_hits, p = 0.5, alternative = "greater")

# Print result
cat("Binomial Test Result:\n")
print(binom_result)

# Build forecast dataframe for plotting
forecast_df <- data.frame(
  Time = test_index,
  Actual = actual_test_values,
  ARIMA100 = forecasts_test_100
)

# Rename columns to match legend
colnames(forecast_df) <- c("Time", "Actual", "ARIMA(1,0,0) Forecast")

# Plot actual vs forecast
ggplot(forecast_df, aes(x = Time)) +
  geom_line(aes(y = Actual, color = "Actual Returns"), linewidth = 0.4) +
  geom_line(aes(y = `ARIMA(1,0,0) Forecast`, color = "ARIMA(1,0,0) Forecast"), linewidth = 0.4) +
  labs(
    title = "Test Set: Actual Returns vs ARIMA(1,0,0) Forecast",
    x = "Time",
    y = "Returns"
  ) +
  scale_color_manual(values = c(
    "Actual Returns" = "black",
    "ARIMA(1,0,0) Forecast" = "red"
  )) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

# Forecast errors on the test set
errors_100 <- actual_test_values - forecasts_test_100

# Ljung-Box test for autocorrelation in errors
Box.test(errors_100, lag = 30, type = "Ljung-Box")

# ACF and PACF plots
acf(errors_100, main = "ACF of Forecast Errors (Test Set)")
pacf(errors_100, main = "PACF of Forecast Errors (Test Set)")

# Fit ARIMA model on training data
arima_100 <- Arima(train_returns, order = c(1,0,0), method = "ML")
resid_100 <- resid(arima_100)

# ARCH effects
plot(resid_100^2)

# Convert to a regular data frame
df <- data.frame(
  Date = index(train_returns),
  SquaredResiduals = as.numeric(resid_100^2)
)

# Plot with years on x-axis
ggplot(df, aes(x = Date, y = SquaredResiduals)) +
  geom_line(color = "black") +
  labs(
    title = "Squared Residuals from ARIMA Model on Training Set",
    subtitle = "Inspection for ARCH effects over time",
    x = "Year",
    y = expression(hat(epsilon)^2)
  ) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(margin = margin(b = 10)),
    axis.title = element_text(face = "bold")
  )


# # ACF and PACF is indicative of ARCH effects being present
acf(resid_100^2)
pacf(resid_100^2)

# We observe ARCH effects
Box.test(resid_100^2, lag = 30, type = "Ljung-Box") # p-value < 0.01
ArchTest(resid_100, lags = 10) # p-value < 0.01


# Ljung-Box test of the squared residuals shows volatility clustering, this suggests
# a GARCH model for volatility forecasting. We will use the ARIMA(1,0,0) model as a foundation of our
# GARCH analyis.

#---------------------Training set/GARCH--------------------------------------

### Analysis GARCH model ###

# Compute returns
full_returns_log <- log(nvda_ts)
full_returns <- na.omit(diff(full_returns_log))

# Define training set (80% of the data)
train_end <- floor(0.8 * length(full_returns))
train_returns <- full_returns[1:train_end]

# Fit ARIMA model on training data
arima_100 <- Arima(train_returns, order = c(1,0,0), method = "ML")
resid_100 <- resid(arima_100)

# Fitting a simple GARCH(1,1) model
garch_11_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0)), # Skip mean modeling here since ARIMA(1,0,0) already modeled the mean – GARCH should only model volatility
  distribution.model = "std" 
)

# Fitting the GARCH(1,1) model
garch_11 <- ugarchfit(spec = garch_11_spec, data = resid(arima_100))  # Using ARIMA residuals

# Printing results
print(garch_11)

# Plot diagnostics
plot(garch_11, which = 1)

# Extract standardized residuals
resid_std <- residuals(garch_11, standardize = TRUE)

# Ljung-Box
Box.test(resid_std, lag = 30, type = "Ljung-Box")
Box.test(resid_std^2, lag = 30, type = "Ljung-Box")

# ARCH LM test
ArchTest(resid_std, lags = 10)

# Extract conditional sigma (volatility forecasts) from GARCH fit
fitted_volatility <- sigma(garch_11)  # conditional standard deviations

# Compute realized volatility proxy (squared returns)
realized_volatility <- residuals(garch_11)^2  # squared residuals as proxy

# Compute RMSE and MAE between predicted and realized volatility
rmse_garch_train <- sqrt(mean((fitted_volatility^2 - realized_volatility)^2))
mae_garch_train <- mean(abs(fitted_volatility^2 - realized_volatility))

# Print the results
cat("Training RMSE (volatility):", rmse_garch_train, "\n")
cat("Training MAE (volatility):", mae_garch_train, "\n")


# Lets make a loop here to find our most suitable GARCH-model according to BIC. 

# Define possible orders for GARCH(p, q)
p_values <- 1:2
q_values <- 1:2

# Create empty dataframe to store results
results <- data.frame(Model = character(), 
                      BIC = numeric(), 
                      stringsAsFactors = FALSE)

# Loop through combinations
for (p in p_values) {
  for (q in q_values) {
    
    # Define GARCH specification
    garch_spec <- ugarchspec(
      variance.model = list(model = "sGARCH", garchOrder = c(p, q)),
      mean.model     = list(armaOrder = c(0, 0), include.mean = FALSE),  # mean is already modeled in ARIMA
      distribution.model = "std"
    )
    
    # Fit the model on ARIMA residuals
    garch_fit <- ugarchfit(spec = garch_spec, data = resid(arima_100))
    
    # Extract AIC and BIC
    bic_value <- infocriteria(garch_fit)[1]
    
    # Store results
    results <- rbind(results, data.frame(
      Model = paste0("GARCH(", p, ",", q, ")"),
      BIC = bic_value
    ))
  }
}

# Sorting models by best BIC
results_BIC <- results[order(results$BIC), ]
print(results_BIC)

# BIC suggests GARCH(1,1)

#--------------------------Test set/GARCH-----------------------------------------

# Expanding window GARCH(1,1) forecasting, forecasting with ARIMA(1,0,0) residuals #

# Define train/test split (80% train, 20% test)
total_length  <- length(full_returns)
split_index   <- floor(0.8 * total_length)      # last index of training set
test_start    <- split_index + 1                # first index of test set
forecast_horizon <- total_length - split_index  # number of out-of-sample points (20% of data)

# Initialize vector for volatility forecasts over test set
vol_forecasts <- numeric(forecast_horizon)

# Expanding-window one-step-ahead volatility forecasts
for (i in 1:forecast_horizon) {
  # Training data indices from start of series up to current point (expanding window)
  train_end   <- split_index + i - 1            # end index of training data for this iteration
  train_data  <- full_returns[1:train_end]      # use all data up to train_end (expanding window)
  
  # Fit ARIMA(1,0,0) to model the mean (AR(1) on returns)
  arima_fit   <- Arima(train_data, order = c(1, 0, 0), include.mean = TRUE, method = "ML")
  arima_resid <- residuals(arima_fit)           # residuals (return series with AR(1) mean removed)
  
  # Specify GARCH(1,1) model for the residuals (mean model has no ARMA terms since mean is filtered by ARIMA)
  garch_spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model     = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "std"
  )
  
  # Fit GARCH(1,1) model to the ARIMA residuals
  garch_fit <- ugarchfit(spec = garch_spec, data = arima_resid, solver = "hybrid")
  
  # Forecast one-step-ahead volatility (conditional standard deviation)
  garch_forecast <- ugarchforecast(garch_fit, n.ahead = 1)
  vol_forecasts[i] <- as.numeric(sigma(garch_forecast))[1]  # extract forecasted sigma and store
}

# Prepare actual vs forecast volatility for the test period
test_dates   <- index(full_returns)[test_start:total_length] # dates for test set
realized_vol <- rollapply(full_returns, width = 30, FUN = sd, align = "right", fill = NA)
realized_vol <- realized_vol[test_start:total_length] # realized volatility (e.g., 30-day rolling SD) over test
vol_compare  <- data.frame(Date = test_dates,
                           ActualVol = realized_vol,
                           ForecastVol = vol_forecasts)

# Compute accuracy metrics 
accuracy_metrics <- accuracy(vol_forecasts, realized_vol)
print(accuracy_metrics)

# Ensure column names are consistent
colnames(vol_compare) <- c("Date", "ActualVol", "ForecastVol")

# Drop any NA rows 
vol_compare_clean <- na.omit(vol_compare)

# Plot actual vs forecasted volatility
ggplot(vol_compare_clean, aes(x = Date)) +
  geom_line(aes(y = ActualVol, color = "Actual Volatility"), linewidth = 0.4) +
  geom_line(aes(y = ForecastVol, color = "Forecasted Volatility"), linewidth = 0.4) +
  labs(
    title = "Test Set: Actual vs GARCH(1,1) Forecasted Volatility",
    x = "Date",
    y = "Volatility (Standard Deviation)"
  ) +
  scale_color_manual(values = c(
    "Actual Volatility" = "black",
    "Forecasted Volatility" = "red"
  )) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

### Diagnostics of GARCH(1,1) Forecast ###

# Check if standardized residuals have autocorrelation

# Get standardized residuals from the fitted GARCH model
std_resid <- residuals(garch_fit, standardize = TRUE) # Extract standardized residuals: removes volatility effects so we can test for autocorrelation

# Perform the Ljung-Box test on standardized residuals
Box.test(std_resid, lag = 30, type = "Ljung-Box") # No autocorrelation - good!

# Plot ACF
acf(std_resid, main = "ACF of Standardized Residuals (GARCH Model)")

# Plot PACF
pacf(std_resid, main = "PACF of Standardized Residuals (GARCH Model)")

ArchTest(std_resid, lags = 10) # GARCH removed ARCH effects 
