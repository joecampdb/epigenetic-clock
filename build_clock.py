import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import cross_val_predict, LeaveOneOut
from sklearn.preprocessing import StandardScaler
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

print("="*60, flush=True)
print("MUSCLE EPIGENETIC CLOCK BUILDER", flush=True)
print("="*60, flush=True)

# 1. Load Data
print("\n[1] LOADING DATA", flush=True)
meth = pd.read_csv('GTEx_Muscle.meth.csv', sep='\t', index_col=0)
print(f"    Methylation: {meth.shape[0]:,} CpGs x {meth.shape[1]} samples", flush=True)

anno = pd.read_csv('GTEx_Muscle.anno.csv', sep='\t', index_col=0).T
print(f"    Annotation: {anno.shape[0]} samples x {anno.shape[1]} variables", flush=True)

def age_midpoint(age_str):
    lo, hi = map(int, age_str.split('-'))
    return (lo + hi) / 2

anno['age_numeric'] = anno['age'].apply(age_midpoint)
print(f"\n    Age: mean={anno['age_numeric'].mean():.1f}, SD={anno['age_numeric'].std():.1f}, range={anno['age_numeric'].min():.0f}-{anno['age_numeric'].max():.0f}", flush=True)

X = meth.T.values.astype(np.float32)
y = anno.loc[meth.columns, 'age_numeric'].values.astype(np.float32)
cpg_names = meth.index.values
del meth

# 2. Correlations
print("\n[2] COMPUTING CpG-AGE CORRELATIONS", flush=True)
X = np.nan_to_num(X, nan=0.5, posinf=1.0, neginf=0.0)
X_mean = X.mean(axis=0)
X_std = X.std(axis=0)
X_std[X_std == 0] = 1
X_norm = (X - X_mean) / X_std
y_norm = (y - y.mean()) / y.std()
correlations = (X_norm.T @ y_norm) / len(y)

print(f"    Correlation range: [{correlations.min():.4f}, {correlations.max():.4f}]", flush=True)

top_pos_idx = np.argsort(correlations)[-5:][::-1]
top_neg_idx = np.argsort(correlations)[:5]

print("\n    Top 5 positive:", flush=True)
for i in top_pos_idx:
    print(f"      {cpg_names[i]}: r = {correlations[i]:.4f}", flush=True)
print("    Top 5 negative:", flush=True)
for i in top_neg_idx:
    print(f"      {cpg_names[i]}: r = {correlations[i]:.4f}", flush=True)

# 3. Feature Selection - MORE STRINGENT to reduce overfitting
print("\n[3] FEATURE SELECTION (stringent: |r| > 0.5)", flush=True)
sig_mask = np.abs(correlations) > 0.5
n_sig = sig_mask.sum()
print(f"    CpGs with |r| > 0.5: {n_sig:,} ({100*n_sig/len(correlations):.2f}%)", flush=True)

X_filtered = X[:, sig_mask]
cpg_filtered = cpg_names[sig_mask]
corr_filtered = correlations[sig_mask]

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_filtered)

# 4. Build Clock with STRONGER regularization
print("\n[4] BUILDING CLOCK (high regularization)", flush=True)
model = ElasticNetCV(
    l1_ratio=[0.5, 0.7, 0.9, 0.95, 0.99, 1.0],  # Favor Lasso for sparsity
    alphas=np.logspace(-2, 2, 100),  # Higher alpha range
    cv=5,
    max_iter=20000,
    random_state=42,
    n_jobs=-1
)

print("    Running LOO-CV (47 folds)...", flush=True)
loo = LeaveOneOut()
y_pred_loo = cross_val_predict(model, X_scaled, y, cv=loo, n_jobs=-1)
model.fit(X_scaled, y)

print(f"    Selected alpha: {model.alpha_:.4f}", flush=True)
print(f"    Selected l1_ratio: {model.l1_ratio_:.2f}", flush=True)
n_coefs = int((model.coef_ != 0).sum())
print(f"    Non-zero coefficients: {n_coefs} / {len(model.coef_)}", flush=True)

# 5. Metrics
print("\n[5] MODEL EVALUATION", flush=True)
r, p = pearsonr(y, y_pred_loo)
mae = np.mean(np.abs(y - y_pred_loo))
rmse = np.sqrt(np.mean((y - y_pred_loo)**2))
r2 = 1 - np.sum((y - y_pred_loo)**2) / np.sum((y - y.mean())**2)

print("="*50, flush=True)
print("MUSCLE CLOCK PERFORMANCE (Leave-One-Out CV)", flush=True)
print("="*50, flush=True)
print(f"Pearson r:        {r:.4f}", flush=True)
print(f"p-value:          {p:.2e}", flush=True)
print(f"R-squared:        {r2:.4f}", flush=True)
print(f"MAE:              {mae:.2f} years", flush=True)
print(f"RMSE:             {rmse:.2f} years", flush=True)
print(f"Clock CpGs:       {n_coefs}", flush=True)
print("="*50, flush=True)

# 6. Overfitting check
print("\n[6] OVERFITTING DIAGNOSTIC", flush=True)
y_pred_train = model.predict(X_scaled)
r_train, _ = pearsonr(y, y_pred_train)
mae_train = np.mean(np.abs(y - y_pred_train))

print(f"    {'Metric':<15} {'Training':<12} {'LOO-CV':<12} {'Gap':<10}", flush=True)
print(f"    {'-'*49}", flush=True)
print(f"    {'Pearson r':<15} {r_train:<12.4f} {r:<12.4f} {r_train-r:<10.4f}", flush=True)
print(f"    {'MAE (years)':<15} {mae_train:<12.2f} {mae:<12.2f} {mae-mae_train:<10.2f}", flush=True)

gap = r_train - r
if gap < 0.1:
    print("    Status: GOOD - Minimal overfitting", flush=True)
elif gap < 0.2:
    print("    Status: ACCEPTABLE - Moderate overfitting for small n", flush=True)
else:
    print("    Status: HIGH - Consider more regularization", flush=True)

# 7. Visualization
print("\n[7] GENERATING PLOTS", flush=True)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

ax = axes[0]
ax.scatter(y, y_pred_loo, alpha=0.7, s=80, edgecolor='black', c='steelblue')
lims = [30, 80]
ax.plot(lims, lims, 'k--', lw=2, label='Perfect prediction')
ax.set_xlim(lims)
ax.set_ylim(lims)
ax.set_xlabel('Chronological Age (years)', fontsize=12)
ax.set_ylabel('Predicted Age (years)', fontsize=12)
ax.set_title(f'Muscle Clock: LOO-CV\nr = {r:.3f}, MAE = {mae:.1f} years', fontsize=14)
ax.legend()
ax.set_aspect('equal')

ax = axes[1]
residuals = y_pred_loo - y
ax.scatter(y, residuals, alpha=0.7, s=80, edgecolor='black', c='coral')
ax.axhline(0, color='black', linestyle='--', lw=2)
ax.set_xlabel('Chronological Age (years)', fontsize=12)
ax.set_ylabel('Residual (Predicted - Actual)', fontsize=12)
ax.set_title(f'Age Acceleration\nMean = {residuals.mean():.2f}, SD = {residuals.std():.2f}', fontsize=14)

plt.tight_layout()
plt.savefig('muscle_clock_performance.png', dpi=150, bbox_inches='tight')
print("    Saved: muscle_clock_performance.png", flush=True)

# Top CpGs
fig, axes = plt.subplots(2, 3, figsize=(12, 8))
top_cpgs_idx = list(top_pos_idx[:3]) + list(top_neg_idx[:3])
for idx, ax in zip(top_cpgs_idx, axes.flat):
    ax.scatter(y, X[:, idx], alpha=0.7, s=60, edgecolor='black')
    ax.set_xlabel('Age')
    ax.set_ylabel('Beta')
    ax.set_title(f'{cpg_names[idx]}\nr = {correlations[idx]:.3f}')
    z = np.polyfit(y, X[:, idx], 1)
    ax.plot(sorted(y), np.polyval(z, sorted(y)), 'r--', lw=2)
plt.suptitle('Top Age-Correlated CpGs', fontsize=14)
plt.tight_layout()
plt.savefig('top_cpgs.png', dpi=150)
print("    Saved: top_cpgs.png", flush=True)

# Correlation histogram
fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(correlations, bins=100, edgecolor='none', alpha=0.7)
ax.axvline(0.5, color='red', linestyle='--', label='|r| = 0.5 threshold')
ax.axvline(-0.5, color='red', linestyle='--')
ax.set_xlabel('Pearson r with Age')
ax.set_ylabel('Number of CpGs')
ax.set_title('CpG-Age Correlation Distribution')
ax.legend()
plt.tight_layout()
plt.savefig('correlation_distribution.png', dpi=150)
print("    Saved: correlation_distribution.png", flush=True)

# Coefficient distribution
if n_coefs > 0:
    fig, ax = plt.subplots(figsize=(8, 5))
    nonzero_coefs = model.coef_[model.coef_ != 0]
    ax.hist(nonzero_coefs, bins=min(30, n_coefs), edgecolor='black', alpha=0.7)
    ax.axvline(0, color='red', linestyle='--', lw=2)
    ax.set_xlabel('Coefficient Value', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(f'Clock Coefficient Distribution (n={len(nonzero_coefs)} CpGs)', fontsize=14)
    plt.tight_layout()
    plt.savefig('coefficient_distribution.png', dpi=150)
    print("    Saved: coefficient_distribution.png", flush=True)

# 8. Save coefficients
print("\n[8] SAVING RESULTS", flush=True)
coef_df = pd.DataFrame({
    'CpG': cpg_filtered,
    'Coefficient': model.coef_,
    'Correlation': corr_filtered
})
coef_df = coef_df[coef_df['Coefficient'] != 0].sort_values('Coefficient', key=abs, ascending=False)
coef_df.to_csv('muscle_clock_coefficients.csv', index=False)
print(f"    Saved: muscle_clock_coefficients.csv ({len(coef_df)} CpGs)", flush=True)

print("\n    Top 15 Clock CpGs:", flush=True)
print(coef_df.head(15).to_string(index=False), flush=True)

# Summary
print("\n" + "="*60, flush=True)
print("FINAL SUMMARY", flush=True)
print("="*60, flush=True)
print(f"Samples:              {len(y)}", flush=True)
print(f"Total CpGs:           {len(cpg_names):,}", flush=True)
print(f"Pre-filtered CpGs:    {n_sig:,}", flush=True)
print(f"Final clock CpGs:     {n_coefs}", flush=True)
print(f"LOO-CV Pearson r:     {r:.4f}", flush=True)
print(f"LOO-CV R-squared:     {r2:.4f}", flush=True)
print(f"LOO-CV MAE:           {mae:.2f} years", flush=True)
print(f"LOO-CV RMSE:          {rmse:.2f} years", flush=True)
print(f"Train-CV Gap (r):     {gap:.4f}", flush=True)
print("="*60, flush=True)
