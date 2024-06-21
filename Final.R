library(forecast)
library(zoo)
library(tidyverse)
library(lubridate)

# Import, Prep, & Plot Data ---------------------------------------------------------------------

path <- here::here("Data/HDTGPDUSQ163N.csv")  
data <- readr::read_csv(path)

data$DATE <- as.Date(data$DATE, format = "%m/%d/%Y")

#View(data)

#Add year and quarter columns
data$Year <- year(data$DATE)
data$Quarter <- quarter(data$DATE)
data <- data %>% select(-DATE)

#rename columns
data <- data[, c("Year", "Quarter", setdiff(names(data), c("Year", "Quarter")))]

#Convert to ratio
data$DEBTGDP <- data$DEBTGDP / 100

#convert to time series
debt.ts <- ts(data$DEBTGDP, start = c(2005,1), freq = 4)

#plot data
plot(debt.ts,
     main = "US Household Debt to GDP (2005-2023)",
     ylab = "Ratio")

#create seasonplot of data
seasonplot(window(debt.ts, start = 2018),
           col = rainbow(5),
           year.labels.left = TRUE, 
           ylab = "Ratio", 
           main = "US Household Debt to GDP (2018-2023)")

# split data into estimation & holdout ~ 75/25 split
estimation <- window(debt.ts, end = c(2018,4))
holdout <- window(debt.ts, start = c(2019,1))

n <- length(debt.ts)
m <- length(holdout)

# Moving Average (MA) Model  ------------------------------------------------------------

ma_order <- c(1,3,5,7)

#Create matrix to store forecast
ma_frc <- matrix(NA, nrow = m, ncol = 4)
colnames(ma_frc) <- paste0("MA ", ma_order)
rownames(ma_frc) <- time(holdout)

#Create forecast
#Loop estimation sample and receive the rolling origin
#For each new origin, calculate the mean for 3 and 7 years
for (i in 1:m) {
  ma_estimate <- debt.ts[1:(n - m + i - 1)]
  for (j in 1:length(ma_order)) {
    ma_frc[i,j] <- mean(tail(ma_estimate, ma_order[j]))
  }
}

#Convert forecast to time series object
ma_frc.ts <- ts(ma_frc, start = start(holdout), freq = frequency(holdout))

#Plot Forecasting Results
plot(debt.ts,ylab = "Debt-GDP")
lines(ma_frc.ts[,1], col = 2)
lines(ma_frc.ts[,2], col = 3)
lines(ma_frc.ts[,3], col = 4)
lines(ma_frc.ts[,4], col = 5)
legend("topright",c("Data", colnames(ma_frc)), col = c("black",2,3),lty = 1)


#Calculate Error
ma_error <- matrix(rep(holdout, 4), ncol = 4) - ma_frc

# MSE & MAE
ma_RMSE <- sqrt(apply(ma_error^2,2,mean))
ma_RMSE

ma_MAE <- apply(abs(ma_error),2,mean)
ma_MAE

# Create matrix to store AIC values
ma_ic <- matrix(NA, nrow = 2, ncol = length(ma_order))
rownames(ma_ic) <- c("AIC","BIC")
colnames(ma_ic) <- paste0("MA ", ma_order)

# Loop over each MA model order
for (j in 1:length(ma_order)) {
  # Calculate error for the j-th MA model
  ma_error <- holdout - ma_frc[,j]
  # Calculate MSE
  ma_mse <- mean(ma_error^2)
  # Calculate AIC using the formula
  ma_ic[1,j] <- length(holdout) * log(ma_mse) + 2 * ma_order[j]
  ma_ic[2,j] <- length(holdout) * log(ma_mse) + log(length(holdout)) * ma_order[j]
}

# Print AIC values for each MA model
ma_ic

# Simple Exponential Smoothing (SES) ----------------------------------------------

# Specify alpha values
alpha_lst <- c(0.1,0.5,0.9)

# Create matrix to store AIC values
ses_ic_values <- matrix(NA, nrow = (length(alpha_lst)+1), ncol = 2)

# Loop over alpha values and fit models
for (j in 1:length(alpha_lst)){
  fit <- ets(debt.ts, model = "ANN", alpha = alpha_lst[j])
  fit$initstate[1] <- window(debt.ts, start = 2005, end = 2005)
  fit <- ets(debt.ts, model = fit, use.initial.values = TRUE)

  # Compute AIC
  ses_ic_values[j,1] <- AIC(fit)
  ses_ic_values[j,2] <- BIC(fit)
}

# let R pick alpha
ses_fit = ets(debt.ts, model="ANN")
alpha_hat = ses_fit$par[1]

# fit model using selected alpha (alpha_hat)
ses_fit_alpha_hat = ets(debt.ts, alpha = alpha_hat, model = "ANN")

# compute AIC for alpha hat
aic_alpha_hat <- AIC(ses_fit_alpha_hat)
bic_alpha_hat <- BIC(ses_fit_alpha_hat)

# store value in AIC matrix
ses_ic_values[4,1] <- aic_alpha_hat
ses_ic_values[4,2] <- bic_alpha_hat

# Rename columns and rows
rownames(ses_ic_values) <- paste("a =", c(alpha_lst, round(alpha_hat,2)))
colnames(ses_ic_values) <- c("AIC","BIC")

# Display AIC for each smoothing parameter (alpha)
ses_ic_values

#Best SES model is with alpha = 0.74

#### Computing SES error using holdout ####
# Specify alpha values
alpha_lst <- c(0.1,0.5,0.9)

ses_frc <- matrix(NA, nrow = m, ncol = 3)
colnames(ses_frc) <- paste("a =", alpha_lst)
rownames(ses_frc) <- time(holdout)
ses_frc

for (i in 1:m){
  ses_estimate <- debt.ts[1:(n-m+i-1)]
  for (j in 1:length(alpha_lst)){
    fit <- ets(ses_estimate, model = "ANN", alpha = alpha_lst[j])
    fit$initstate[1] <- window(debt.ts, start = 2005, end = 2005)
    fit <- ets(ses_estimate, model = fit, use.initial.values = TRUE)
    ses_frc[i,j] <- forecast(fit, h = 1)$mean
  }
}
ses_frc

#Convert SES forecast to time series
ses.ts <- ts(ses_frc,frequency=frequency(holdout),start=start(holdout))

#Plot SES forecast
plot(debt.ts, ylab="Debt-GDP")
lines(ses.ts[,1], col=2)
lines(ses.ts[,2], col=3)
lines(ses.ts[,3], col=4)
legend("right",c("Data",colnames(ses_frc)),col=c("black",2,3,4),lty=2)

#Compute error
ses_error <- matrix(rep(holdout,3), ncol = 3) - ses_frc
ses_error

#Compute RMSE
ses_RMSE <- sqrt(apply(ses_error^2, 2, mean))
ses_RMSE

#Compute MAE
ses_MAE <- apply(abs(ses_error), 2, mean)
ses_MAE

# let R pick alpha

ses_fit = ets(estimation, model="ANN")
alpha_hat = ses_fit$par[1]
alpha_hat

ses_fit_alpha_hat = ets(debt.ts, alpha = alpha_hat, model = "AAN")

ses_ah_error = holdout - window(ses_fit_alpha_hat$fitted, start = c(2019,1))

mse_alpha_hat = sqrt(mean(ses_ah_error^2))
mse_alpha_hat

mae_alpha_hat = mean(abs(ses_ah_error))
mae_alpha_hat

# De-Seasonalized LES -----------------------------------------------------------

# Method 1: Forecasting Decomposition

# Decompose data into seasonal and trend components
components <- decompose(debt.ts)

plot(components$seasonal, main = "Seasonal Component Plot", ylab = "Change")
plot(components$trend, main = "Trend Component Plot", ylab = "Debt-GDP")

# Remove seasonal component from data
deseasonal_1 <- debt.ts - components$seasonal

# Fit deseasonal LES model
ds_fit_1 <- ets(deseasonal_1, model = "AAN") 

# Compute AIC
ds_aic_1 <- AIC(ds_fit_1)
ds_bic_1 <- BIC(ds_fit_1)

#### Forecasting w/ deasonal fit 1 ####
# Error
ds_error_1 <- holdout - window(ds_fit_1$fitted + components$seasonal, start = c(2019,1))
ds_error_1

rmse_ds_1 <- sqrt(mean(ds_error_1^2))

mae_ds_1 <- mean(abs(ds_error_1))

c(rmse_ds_1,mae_ds_1)

frc_ds_1 <- forecast(ds_fit_1, h = 12)$mean + components$seasonal[1:4]

# Plot the original data
plot(debt.ts, xlim = c(2005, 2028), ylim = range(c(debt.ts, frc_ds_1)),
     xlab = "Year", ylab = "Debt-GDP", main = "Forecasting Decomposition Forecast (H=12)")

# Add the forecasted values
lines(frc_ds_1, col = 2)

# Add legend
legend("topright", legend = c("Original Data", "Forecast"), col = c("black", "blue"), lty = 1)

#### ####

# Method 2: Pure Decomposition (Centered Moving Average)

# Create centered moving average by quarters
CMA <- ma(debt.ts, order = 4)

# Detrending
detrended <- debt.ts / CMA

# Deseasonalization
seasonal <- sapply(1:4, function(q) mean(detrended[seq(q, length(detrended), 4)], na.rm = TRUE))
adj_seasonal <- seasonal / mean(seasonal)
rep_adj_seasonal <- rep(adj_seasonal, length.out = length(debt.ts))
deseasonal_2 <- debt.ts / rep_adj_seasonal

# Fit a LES model on the deseasonalized series
ds_fit_2 <- ets(deseasonal_2, model = "AAN")

# Compute AIC for deseasonal fit
ds_aic_2 <- AIC(ds_fit_2)
ds_bic_2 <- BIC(ds_fit_2)

ds_ic <- matrix(c(ds_aic_1,ds_bic_1,ds_aic_2,ds_bic_2), nrow = 2, ncol = 2)
colnames(ds_ic) <- c("Forecast", "Pure")
rownames(ds_ic) <- c("AIC", "BIC")

ds_ic

# Best model is with forecasting decomposition

# Holt-Winters --------------------------------------------------------------------

# let R choose the model 

hw_auto_fit <- hw(debt.ts)
hw_auto_fit$model$aic
hw_auto_fit$method

hw_auto_fit$model$par[1:3]

a = hw_auto_fit$model$par[1]
b = hw_auto_fit$model$par[2]
g = hw_auto_fit$model$par[3]

# g is so small it almost does nothing

hw_auto_fit <- ets(debt.ts, model = "AAA", alpha = a, beta = b, gamma = g)
hw_auto_fit_aic <- AIC(hw_auto_fit)
hw_auto_fit_bic <- BIC(hw_auto_fit)
  
# Additive no seasonal

hw_additive_fit <- ets(debt.ts, model = "AAN", damped = F)

#what parameters are chosen
hw_additive_fit$par[1:2]

hw_add_aic <- AIC(hw_additive_fit)
hw_add_bic <- BIC(hw_additive_fit)

# Damped Additive

hw_damped_add <- ets(debt.ts, model = "AAA", damped = T)

hw_damped_add$par[1:4]

hw_damped_add_aic <- AIC(hw_damped_add)
hw_damped_add_bic <- BIC(hw_damped_add)

# Multiplicative
hw_mult_fit <- ets(debt.ts, model = "MAM", damped = F)

hw_mult_fit$par[1:3]

hw_mult_aic <- AIC(hw_mult_fit)
hw_mult_bic <- BIC(hw_mult_fit)

#Best model is additive holt-winters

hw_ic <- matrix(c(hw_auto_fit_aic, hw_add_aic, hw_damped_add_aic, hw_mult_aic,
                         hw_auto_fit_bic, hw_add_bic, hw_damped_add_bic, hw_mult_bic),
                       nrow = 2, ncol = 4, byrow = TRUE,
                       dimnames = list(c("AIC", "BIC"),
                                       c("Auto", "Additive", "Damped Add", "Multiplicative")))
hw_ic
# ARIMA ---------------------------------------------------------------------------

#shows seasonality and trend
plot(debt.ts)

Acf(debt.ts)
#autocorrelation plot confirms non-stationarity of data

#Must make data stationary by differencing
#First order differencing
diff_1st_debt_ts <- diff(debt.ts, lag = 4)
plot(diff_1st_debt_ts)
Acf(diff_1st_debt_ts)

#Still exhibits non-stationarity, use 2nd order differencing
diff_2nd_debt_ts <- diff(diff(debt.ts, lag = 4))
plot(diff_2nd_debt_ts)
#data appears stationary now, can use for ARIMA

#now we look at acf & pacf plots to determine p,d,q for our model
Acf(diff_2nd_debt_ts)
#tails off as damped wave pattern
pacf(diff_2nd_debt_ts)
#tails off linearly

# p is for AR scheme, d for differencing, q for MA scheme
#AR scheme of 1 from plots
p <- 1
#second order differencing, Differencing scheme of 2
d <- 2
#MA scheme is 1
q <- 1

#the combination of parameters above provides the best log likelihood/lowest AIC
arima_1 <- arima(debt.ts, order = c(p,d,q))
arima_1

par(mar = c(5, 5, 4, 2) + 0.1)

#now we should analyze the residuals
plot(x = 1:length(arima_1$residuals), y = arima_1$residuals, 
     main = "Scatterplot of ARIMA 1 Residuals", xlab = "Indexed Quarters", ylab = "Residuals", 
     pch = 16, col = "black")

#significant spike at lag 4, 8, 12 ... indicates seasonality still uncaptured
Acf(arima_1$residuals)

#now we need to include the seasonal components

#seasonal parameters
P <- 1
D <- 0
Q <- 1

#combination of seasonal parameters above provides best fit 
arima_2 <- arima(debt.ts, order = c(p,d,q), seasonal = c(P,D,Q))
arima_2

AIC(arima_2)

plot(x = 1:length(arima_2$residuals), y = arima_2$residuals, 
     main = "Scatterplot of ARIMA 2 Residuals", xlab = "Indexed Quarters", ylab = "Residuals", 
     pch = 16, col = "black")

#alternatively, we could use the auto.arima function to choose parameters for us

arima_3 <- auto.arima(debt.ts)
arima_3

plot(x = 1:length(arima_3$residuals), y = arima_3$residuals, 
     main = "Scatterplot of ARIMA 3 Residuals", xlab = "Indexed Quarters", ylab = "Residuals", 
     pch = 16, col = "black")

#this model uses p = 0, d = 1, q = 0, and P = 0, D = 1, Q = 1

#Calculate AIC for each model to choose the best model
arima_1_aic <- AIC(arima_1)
arima_1_bic <- BIC(arima_1)
arima_2_aic <- AIC(arima_2)
arima_2_bic <- BIC(arima_2)
arima_3_aic <- AIC(arima_3)
arima_3_bic <- BIC(arima_3)

arima_ic <- matrix(c(arima_1_aic,arima_2_aic,arima_3_aic,arima_1_bic,arima_2_bic,
                     arima_3_bic), nrow = 2, ncol = 3, byrow = TRUE,
                   dimnames = list(c("AIC", "BIC"),
                                   c("Model 1", "Model 2", "Model 3")))

arima_ic
# Calculate error for ARIMA model 2

arima_refit <- arima(estimation, order = c(1,2,1), seasonal = c(1,1,1))
arima_frc <- forecast(arima_refit, h = length(holdout))
arima_frc_point <- arima_frc$mean

arima_frc.ts <- ts(arima_frc_point, start = start(holdout), end = end(holdout), frequency = 4)

plot(debt.ts)
lines(arima_frc.ts, col = 2)

arima_error <- holdout - arima_frc.ts

rmse_arima <- sqrt(mean(arima_error^2))
mae_arima <- mean(abs(arima_error))

c(rmse_arima_3,mae_arima_3)

#Model 3 diagnostics
# Decrease plot margins
par(mar = c(3, 3, 1, 1))  # Set bottom, left, top, and right margins
tsdiag(arima_2)


#Check normality of the residuals
qqnorm(arima_2$residuals)
qqline(arima_2$residuals)
hist(arima_2$residuals)

# Model Selection -----------------------------------------------------------------

#Compare AIC from the best model using each method

# Create a matrix to store AIC values
model_names <- c("MA", "SES", "Deseasonal LES", "Holt-Winters", "ARIMA")
best_AIC_values <- c(ma_ic[1,1], ses_ic_values[4,1], ds_aic_1, hw_damped_add_aic, arima_2_aic)
best_BIC_values <- c(ma_ic[2,1], ses_ic_values[4,2], ds_bic_1, hw_damped_add_bic, arima_2_bic)
comparison_matrix <- matrix(c(best_AIC_values,best_BIC_values), nrow = length(model_names), ncol = 2)

# Add row and column names
rownames(comparison_matrix) <- model_names
colnames(comparison_matrix) <- c("AIC","BIC")

# Print the comparison matrix
print(comparison_matrix)



# Forecasting  -----------------------------------------------------

h4_forecast <- forecast(arima_2, h = 4)

h4_forecast

plot(h4_forecast, main = "Forecast Next 4 Quarters",
     ylab = "Debt-GDP Ratio",
     xlab = "Time",
     type = "l")

h8_forecast <- forecast(arima_3, h = 8)

plot(h8_forecast, main = "Forecast Next 8 Quarters",
     ylab = "Debt-GDP Ratio",
     xlab = "Time",
     type = "l")

h12_forecast <- forecast(arima_2, h = 12)

plot(h12_forecast, main = "Forecast Next 12 Quarters",
     ylab = "Debt-GDP Ratio",
     xlab = "Time",
     type = "l")

h12_forecast


debt.ts
