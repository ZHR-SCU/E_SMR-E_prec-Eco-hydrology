import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import shap
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error

# === 1. Load example dataset ===
# Replace with your actual data file path
file_path = 'your_data.xlsx'
df = pd.read_excel(file_path, sheet_name='YourSheetName', engine='openpyxl')

# === 2. Select variables ===
# Replace with actual column names in your dataset
ALT = df['ALT_column']
Prec = df['Prec_column']
Target1 = df['Target_SM_bottom']
Target2 = df['Target_SM_root']

# === 3. Filter missing data ===
valid_idx = ~ALT.isna() & ~Prec.isna() & ~Target1.isna() & ~Target2.isna()
ALT = ALT[valid_idx].values
Prec = Prec[valid_idx].values
Target1 = Target1[valid_idx].values
Target2 = Target2[valid_idx].values

# === 4. Standardize predictors ===
scaler = StandardScaler()
Prec_scaled = scaler.fit_transform(Prec.reshape(-1, 1))
ALT_scaled = scaler.fit_transform(ALT.reshape(-1, 1))
interaction_scaled = Prec_scaled * ALT_scaled

# Combine predictors
X = np.hstack([Prec_scaled, ALT_scaled, interaction_scaled])
feature_names = ['Precipitation (scaled)', 'ALT (scaled)', 'Precip × ALT (scaled)']

# === 5. Standardize targets ===
scaler_y1 = StandardScaler()
scaler_y2 = StandardScaler()
Target1_scaled = scaler_y1.fit_transform(Target1.reshape(-1, 1)).flatten()
Target2_scaled = scaler_y2.fit_transform(Target2.reshape(-1, 1)).flatten()

# === 6. Fit Random Forest for Target1 ===
rf1 = RandomForestRegressor(n_estimators=10000, random_state=0)
rf1.fit(X, Target1_scaled)
pred1 = rf1.predict(X)
mse1 = mean_squared_error(Target1_scaled, pred1)
print(f"\n=== OOB MSE (Target1) === {mse1:.4f} ===\n")

# === 7. Fit Random Forest for Target2 ===
rf2 = RandomForestRegressor(n_estimators=10000, random_state=0)
rf2.fit(X, Target2_scaled)
pred2 = rf2.predict(X)
mse2 = mean_squared_error(Target2_scaled, pred2)
print(f"\n=== OOB MSE (Target2) === {mse2:.4f} ===\n")

# === 8. SHAP analysis for Target1 ===
explainer1 = shap.TreeExplainer(rf1)
shap_values1 = explainer1.shap_values(X)

mean_shap1 = np.mean(shap_values1, axis=0)
abs_mean_shap1 = np.mean(np.abs(shap_values1), axis=0)

print("\n=== SHAP Contribution (Target1) ===\n")
for i, name in enumerate(feature_names):
    print(f"{name}: Mean SHAP = {mean_shap1[i]:.4f}, Mean |SHAP| = {abs_mean_shap1[i]:.4f}")

# Feature importance from Random Forest
importances = rf1.feature_importances_
importances_pct = 100 * importances / np.sum(importances)

print("\n=== Feature Importances (Target1) ===\n")
for i, name in enumerate(feature_names):
    print(f"{name}: {importances_pct[i]:.2f}%")

# SHAP summary plot
shap.summary_plot(shap_values1, X, feature_names=feature_names)

# === 9. SHAP analysis for Target2 ===
explainer2 = shap.TreeExplainer(rf2)
shap_values2 = explainer2.shap_values(X)

mean_shap2 = np.mean(shap_values2, axis=0)
abs_mean_shap2 = np.mean(np.abs(shap_values2), axis=0)

print("\n=== SHAP Contribution (Target2) ===\n")
for i, name in enumerate(feature_names):
    print(f"{name}: Mean SHAP = {mean_shap2[i]:.4f}, Mean |SHAP| = {abs_mean_shap2[i]:.4f}")

# Feature importance from Random Forest
importances2 = rf2.feature_importances_
importances2_pct = 100 * importances2 / np.sum(importances2)

print("\n=== Feature Importances (Target2) ===\n")
for i, name in enumerate(feature_names):
    print(f"{name}: {importances2_pct[i]:.2f}%")

# SHAP summary plot
shap.summary_plot(shap_values2, X, feature_names=feature_names)
