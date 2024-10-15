import numpy as np
import pandas as pd
from statsmodels.tsa.arima.model import ARIMA
from sklearn.metrics import root_mean_squared_error
import warnings

# Suppress all warnings
warnings.filterwarnings("ignore")


df = pd.read_csv("C:/Users/matte/Documents/STO MOD AND SIM/Progetto Esame/Data_wolves_moose_Isle_Royale_June2019.csv", header=None, sep = ";")

# Select the relevant columns
df = df[[0, 1, 2]]

# Rename specific columns
df.rename(columns={0: 'year', 1: 'pred_abs', 2: 'prey_abs'}, inplace=True)


# Calculate the 'Sum' column
df['Sum'] = df['pred_abs'] + df['prey_abs']

# Calculate the 'prey' and 'pred' columns
df['prey'] = df['prey_abs'] / df['Sum']
df['pred'] = df['pred_abs'] / df['Sum']

# Select the desired columns
df_selected = df[['prey', 'pred']]
print(df_selected.head())


train_arima = df_selected.iloc[:40]  
val_arima = df_selected.iloc[40:50]  
full_train = df_selected.iloc[:50] 




train_arima.reset_index(drop=True, inplace=True)
val_arima.reset_index(drop=True, inplace=True)
full_train.reset_index(drop=True, inplace=True)

best_rmse_prey = float('inf')
best_params_prey = (0, 0, 0)
best_rmse_pred = float('inf')
best_params_pred = (0, 0, 0)


for p in range(6):  # Try p from 0 to 5
    for d in range(2):  # Try d from 0 to 1
        for q in range(6):  # Try q from 0 to 5
            try:
                # Fit ARIMA(p, d, q)
                model = ARIMA(train_arima.prey, order=(p, d, q))
                model_fit = model.fit()
                
                # Forecast on the test set
                forecast = model_fit.forecast(steps=len(val_arima))
                
                # Compute the Root Mean Squared Error (RMSE)
                rmse = root_mean_squared_error(val_arima.prey, forecast)
                
                # Update best parameters if this is the lowest RMSE
                if rmse < best_rmse_prey:
                    best_rmse_prey = rmse
                    best_params_prey = (p, d, q)
            except Exception as e:
                # Handle any fitting errors (e.g., non-invertible models)
                print(f"Error fitting prey ARIMA({p}, {d}, {q}): {e}")
                continue


model_prey = ARIMA(full_train.prey, order=(best_params_prey[0], best_params_prey[1], best_params_prey[2]))  
fitted_prey = model_prey.fit() 
fc_prey = fitted_prey.forecast(11)
fc_prey_series = pd.Series(fc_prey)


for p in range(6):  # Try p from 0 to 5
    for d in range(2):  # Try d from 0 to 1
        for q in range(6):  # Try q from 0 to 5
            try:
                # Fit ARIMA(p, d, q)
                model = ARIMA(train_arima.pred, order=(p, d, q))
                model_fit = model.fit()
                
                # Forecast on the test set
                forecast = model_fit.forecast(steps=len(val_arima))
                
                # Compute the Root Mean Squared Error (RMSE)
                rmse = root_mean_squared_error(val_arima.pred, forecast)
                
                # Update best parameters if this is the lowest RMSE
                if rmse < best_rmse_pred:
                    best_rmse_pred = rmse
                    best_params_pred = (p, d, q)
            except Exception as e:
                # Handle any fitting errors (e.g., non-invertible models)
                print(f"Error fitting pred ARIMA({p}, {d}, {q}): {e}")
                continue


model_pred = ARIMA(full_train.pred, order=(best_params_pred[0], best_params_pred[1], best_params_pred[2]))  
fitted_pred = model_pred.fit() 
fc_pred = fitted_pred.forecast(11)
fc_pred_series = pd.Series(fc_pred)

df = pd.DataFrame({'prey fore': fc_prey_series, 'pred fore': fc_pred_series})
df.to_csv('ARIMA.csv', index=False)
