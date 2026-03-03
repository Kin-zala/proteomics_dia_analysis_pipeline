import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import os

# Ensure results folder exists
os.makedirs("results", exist_ok=True)

# --- 1. Load data ---
df = pd.read_csv("data/proteinGroups.txt", sep="\t", low_memory=False)
print(f"Raw data shape: {df.shape}")

# --- 2. Filter contaminants and reverse hits ---
df = df[
    df['Reverse'].isna() &
    df['Potential contaminant'].isna() &
    df['Only identified by site'].isna()
].copy()
print(f"After filtering contaminants: {df.shape}")

# --- 3. Select LFQ intensity columns ---
lfq_cols = [col for col in df.columns if col.startswith("LFQ intensity ")]
data = df[lfq_cols].replace(0, np.nan)
data.columns = data.columns.str.replace("LFQ intensity ", "", regex=False)
data.index = df['Gene names'].fillna(df['Protein IDs']).astype(str)

# --- 4. Log2 transform and Filter ---
log_data = np.log2(data)
# Keep proteins with at least 50% valid values across all samples
log_data = log_data.dropna(thresh=int(len(log_data.columns) / 2))
print(f"After validity filter: {log_data.shape}")

# --- 5. Downshift Imputation (Perseus-style) ---
def downshift_imputation(df, width=0.3, shift=1.8):
    df_imputed = df.copy()
    for col in df_imputed.columns:
        valid = df_imputed[col].dropna()
        if len(valid) > 1:
            mu = valid.mean()
            sigma = valid.std()
            impute_mean = mu - (shift * sigma)
            impute_std = sigma * width
            mask = df_imputed[col].isna()
            df_imputed.loc[mask, col] = np.random.normal(
                impute_mean, impute_std, size=mask.sum()
            )
    return df_imputed

log_data = downshift_imputation(log_data)

# Safety fill for any remaining NaNs
if log_data.isna().any().any():
    print("Detected remaining NaNs. Applying global minimum imputation...")
    global_min = log_data.min().min()
    log_data = log_data.fillna(global_min)

print(f"Remaining NaNs: {log_data.isna().sum().sum()}")
print(f"Data shape for analysis: {log_data.shape}")

# --- 6. Define Groups ---
samples = log_data.columns.tolist()

group_undiff = [c for c in samples if "undiff" in c]
group_RA     = [c for c in samples if c.startswith("RA")]
group_PMA    = [c for c in samples if "PMA" in c]

print(f"\nGroup sizes:")
print(f"  Undiff: {len(group_undiff)} -> {group_undiff}")
print(f"  RA:     {len(group_RA)}     -> {group_RA}")
print(f"  PMA:    {len(group_PMA)}    -> {group_PMA}")

# --- 7. PCA Analysis ---
pca = PCA(n_components=2)
pca_input = log_data.T
pca_input.columns = pca_input.columns.astype(str)
pca_result = pca.fit_transform(pca_input)

# Color by group
groups = []
for c in pca_input.index:
    if "undiff" in c:       groups.append("Undiff")
    elif c.startswith("RA"): groups.append("RA")
    elif "PMA" in c:         groups.append("PMA")
    else:                    groups.append("Other")

unique_groups = ["Undiff", "RA", "PMA"]
colors = ["#2196F3", "#4CAF50", "#F44336"]
color_map = dict(zip(unique_groups, colors))

plt.figure(figsize=(10, 7))
for group, color in zip(unique_groups, colors):
    idx = [i for i, g in enumerate(groups) if g == group]
    plt.scatter(pca_result[idx, 0], pca_result[idx, 1],
                label=group, s=120, color=color, edgecolors='k', alpha=0.85)

for i, sample in enumerate(pca_input.index):
    plt.text(pca_result[i, 0] + 0.3, pca_result[i, 1], str(sample), fontsize=7)

plt.title("PCA: SH-SY5Y Undiff vs RA vs PMA (Log2 LFQ)", fontsize=13)
plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
plt.legend(title="Condition")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("results/pca.png", dpi=150)
plt.close()
print("PCA plot saved.")

# --- 8. Differential Analysis (Welch's T-test) ---
def run_ttest(log_data, group1, group2, label):
    results_list = []
    for i in range(len(log_data)):
        g1_v = log_data.iloc[i][group1]
        g2_v = log_data.iloc[i][group2]
        stat, p = ttest_ind(g1_v, g2_v, equal_var=False, nan_policy='omit')
        results_list.append({
            "protein": log_data.index[i],
            "pvalue": p,
            "log2FC": g2_v.mean() - g1_v.mean()
        })
    results = pd.DataFrame(results_list).dropna(subset=['pvalue'])
    results["adj_pvalue"] = multipletests(results["pvalue"], method="fdr_bh")[1]
    sig_mask = (results["adj_pvalue"] < 0.05) & (abs(results["log2FC"]) > 1)
    print(f"Significant proteins ({label}): {sig_mask.sum()}")
    results.to_csv(f"results/all_proteins_{label}.csv", index=False)
    results[sig_mask].sort_values("adj_pvalue").to_csv(
        f"results/significant_proteins_{label}.csv", index=False)
    return results

results_RA  = run_ttest(log_data, group_undiff, group_RA,  "undiff_vs_RA")
results_PMA = run_ttest(log_data, group_undiff, group_PMA, "undiff_vs_PMA")

# --- 9. Volcano Plots ---
def plot_volcano(results, label, title):
    results = results.copy()
    results["color"] = "grey"
    results.loc[(results["log2FC"] > 1)  & (results["adj_pvalue"] < 0.05), "color"] = "red"
    results.loc[(results["log2FC"] < -1) & (results["adj_pvalue"] < 0.05), "color"] = "blue"

    sig_count = (results["color"] != "grey").sum()

    plt.figure(figsize=(9, 6))
    for color, label_name in zip(["grey", "red", "blue"], ["Not significant", "Upregulated", "Downregulated"]):
        subset = results[results["color"] == color]
        plt.scatter(subset["log2FC"], -np.log10(subset["adj_pvalue"]),
                    c=color, alpha=0.5, s=12, label=label_name)

    plt.axvline(x=1,  color='black', linestyle='--', linewidth=0.8)
    plt.axvline(x=-1, color='black', linestyle='--', linewidth=0.8)
    plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', linewidth=0.8)
    plt.xlabel("log2 Fold Change", fontsize=12)
    plt.ylabel("-log10 Adjusted P-value", fontsize=12)
    plt.title(f"{title}\nSignificant proteins: {sig_count}", fontsize=13)
    plt.legend(fontsize=9)
    plt.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(f"results/volcano_{label}.png", dpi=150)
    plt.close()
    print(f"Volcano plot saved: results/volcano_{label}.png")

plot_volcano(results_RA,  "undiff_vs_RA",  "Undiff vs RA-differentiated")
plot_volcano(results_PMA, "undiff_vs_PMA", "Undiff vs RA+PMA-differentiated")

# --- 10. Heatmap of top significant proteins ---
def plot_heatmap(results, log_data, group1, group2, label, top_n=50):
    sig = results[(results["adj_pvalue"] < 0.05) & (abs(results["log2FC"]) > 1)]
    sig = sig.sort_values("adj_pvalue").head(top_n)
    if len(sig) == 0:
        print(f"No significant proteins for heatmap ({label})")
        return
    cols = group1 + group2
    heatmap_data = log_data.loc[sig["protein"], cols]
    plt.figure(figsize=(12, max(6, len(sig) * 0.25)))
    sns.heatmap(heatmap_data, cmap="RdBu_r", center=heatmap_data.mean().mean(),
                yticklabels=True, xticklabels=True, linewidths=0.3)
    plt.title(f"Top {len(sig)} significant proteins: {label}", fontsize=12)
    plt.tight_layout()
    plt.savefig(f"results/heatmap_{label}.png", dpi=150)
    plt.close()
    print(f"Heatmap saved: results/heatmap_{label}.png")

plot_heatmap(results_RA,  log_data, group_undiff, group_RA,  "undiff_vs_RA")
plot_heatmap(results_PMA, log_data, group_undiff, group_PMA, "undiff_vs_PMA")

# --- 11. Export matrix for limma ---
log_data_export = log_data[group_undiff + group_RA + group_PMA].copy()
# Make index unique to avoid duplicate row.names error in R
log_data_export.index = [f"{name}__{i}" for i, name in enumerate(log_data_export.index)]
log_data_export.to_csv("results/log2_lfq_matrix_for_limma.csv")
print("\nExported matrix for limma.")
print("Python analysis complete.")