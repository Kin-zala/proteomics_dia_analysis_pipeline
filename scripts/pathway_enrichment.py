import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gseapy as gp
import os

os.makedirs("results/enrichment", exist_ok=True)

# ─── 1. Load significant proteins ───────────────────────────────────────────

def load_sig_proteins(filepath):
    df = pd.read_csv(filepath)
    # Strip any __index suffix added during export
    df["protein"] = df["protein"].str.replace(r"__\d+$", "", regex=True)
    # Some entries have multiple gene names (e.g. "GENE1;GENE2") — take first
    df["protein"] = df["protein"].str.split(";").str[0]
    return df["protein"].dropna().unique().tolist()

sig_RA  = load_sig_proteins("results/significant_proteins_undiff_vs_RA.csv")
sig_PMA = load_sig_proteins("results/significant_proteins_undiff_vs_PMA.csv")

print(f"Significant proteins - RA:  {len(sig_RA)}")
print(f"Significant proteins - PMA: {len(sig_PMA)}")

# ─── 2. Run Enrichment Analysis ──────────────────────────────────────────────

GENE_SETS = [
    "GO_Biological_Process_2023",
    "GO_Cellular_Component_2023",
    "KEGG_2021_Human",
    "Reactome_2022",
]

def run_enrichment(gene_list, label):
    print(f"\nRunning enrichment for: {label} ({len(gene_list)} proteins)")
    all_results = []

    for gene_set in GENE_SETS:
        try:
            enr = gp.enrichr(
                gene_list  = gene_list,
                gene_sets  = gene_set,
                organism   = "human",
                outdir     = None,
                cutoff     = 0.05,
            )
            df = enr.results
            df["gene_set"] = gene_set
            all_results.append(df)
            sig = df[df["Adjusted P-value"] < 0.05]
            print(f"  {gene_set}: {len(sig)} significant terms")
        except Exception as e:
            print(f"  {gene_set}: failed ({e})")

    if not all_results:
        print(f"No enrichment results for {label}")
        return None

    combined = pd.concat(all_results, ignore_index=True)
    combined = combined[combined["Adjusted P-value"] < 0.05]
    combined = combined.sort_values("Adjusted P-value")
    combined.to_csv(f"results/enrichment/enrichment_{label}.csv", index=False)
    print(f"  Total significant terms saved: {len(combined)}")
    return combined

results_RA  = run_enrichment(sig_RA,  "RA_vs_Undiff")
results_PMA = run_enrichment(sig_PMA, "PMA_vs_Undiff")

# ─── 3. Dotplot of Top Pathways ──────────────────────────────────────────────

def plot_dotplot(enr_df, label, title, top_n=20):
    if enr_df is None or len(enr_df) == 0:
        print(f"No data to plot for {label}")
        return

    df = enr_df.sort_values("Adjusted P-value").head(top_n).copy()

    # Clean up term names
    df["Term"] = df["Term"].str.replace(r"\(GO:\d+\)", "", regex=True).str.strip()
    df["Term"] = df["Term"].str[:60]

    df["-log10(FDR)"]   = -np.log10(df["Adjusted P-value"].clip(lower=1e-10))
    df["Overlap_ratio"] = df["Overlap"].apply(
        lambda x: int(x.split("/")[0]) / int(x.split("/")[1]) if "/" in str(x) else 0
    )
    df["Gene_count"] = df["Overlap"].apply(
        lambda x: int(x.split("/")[0]) if "/" in str(x) else 0
    )
    df = df.sort_values("-log10(FDR)", ascending=True)

    fig, ax = plt.subplots(figsize=(11, max(6, len(df) * 0.38)))

    scatter = ax.scatter(
        df["-log10(FDR)"],
        range(len(df)),
        c=df["Overlap_ratio"],
        s=df["Gene_count"] * 10,
        cmap="RdYlBu_r",
        alpha=0.85,
        edgecolors="grey",
        linewidths=0.4,
        vmin=0, vmax=1
    )

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["Term"], fontsize=9)
    ax.set_xlabel("-log10(Adjusted P-value)", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=12)
    ax.axvline(x=-np.log10(0.05), color="black", linestyle="--",
               linewidth=0.8, alpha=0.6)
    ax.grid(axis="x", alpha=0.25)

    cbar = plt.colorbar(scatter, ax=ax, pad=0.01)
    cbar.set_label("Gene Ratio", fontsize=9)

    # Size legend
    handles = [plt.scatter([], [], s=s*10, c="grey", alpha=0.6, label=f"{s} genes")
               for s in [5, 10, 20]]
    ax.legend(handles=handles, title="Gene count", fontsize=8, loc="lower right")

    plt.tight_layout()
    plt.savefig(f"results/enrichment/dotplot_{label}.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Dotplot saved: results/enrichment/dotplot_{label}.png")

plot_dotplot(results_RA,  "RA_vs_Undiff",  "Top Enriched Pathways: RA-differentiated vs Undiff")
plot_dotplot(results_PMA, "PMA_vs_Undiff", "Top Enriched Pathways: RA+PMA-differentiated vs Undiff")

# ─── 4. Highlight DLR-relevant Pathways ─────────────────────────────────────

DLR_KEYWORDS = [
    "cytoskeleton", "actin", "tubulin", "microtubule",
    "neurodegeneration", "neuron", "axon", "synap",
    "mitochondri", "oxidative", "stress", "autophagy",
    "ubiquitin", "proteasome", "apoptosis",
    "metabolic", "glycolysis", "ATP",
    "muscle", "motor"
]

def filter_dlr_relevant(enr_df, label):
    if enr_df is None or len(enr_df) == 0:
        return
    pattern = "|".join(DLR_KEYWORDS)
    relevant = enr_df[enr_df["Term"].str.contains(pattern, case=False, na=False)]
    relevant = relevant.sort_values("Adjusted P-value")
    relevant.to_csv(f"results/enrichment/DLR_relevant_pathways_{label}.csv", index=False)
    print(f"\nDLR-relevant pathways ({label}): {len(relevant)}")
    if len(relevant) > 0:
        print(relevant[["Term", "Adjusted P-value", "Overlap", "gene_set"]].head(15).to_string(index=False))

filter_dlr_relevant(results_RA,  "RA_vs_Undiff")
filter_dlr_relevant(results_PMA, "PMA_vs_Undiff")

print("\nPathway enrichment analysis complete.")
print("Results saved in: results/enrichment/")