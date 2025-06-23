import scanpy as sc
import pandas as pd

# 1. load data
adata = sc.read_h5ad("Amygdala.h5ad", backed="r")

# 2. choose superlcuster of your interest
target_clusters = ["Amygdala excitatory", "Hippocampal CA1-3"]  
target_idx = adata.obs[adata.obs["supercluster_term"].isin(target_clusters)].index

adata_sub = adata[target_idx].to_memory()

# 3. extract tsne adn metadata
tsne_df = pd.DataFrame({
    "CellID": adata_sub.obs_names,
    "TSNE1": adata_sub.obsm["X_tSNE"][:, 0],
    "TSNE2": adata_sub.obsm["X_tSNE"][:, 1],
    "Region": adata_sub.obs["supercluster_term"],  # Amygdala/Hippocampus
    "Cluster": adata_sub.obs["cluster_id"],       # Np. "sldsc_v154"
     
})

# 4. save
tsne_df[tsne_df["Region"] == "Amygdala excitatory"].to_csv("amygdala_tsne.csv.gz", index=False)

tsne_df[tsne_df["Region"] == "Hippocampal CA1-3"].to_csv("hipp_tsne.csv.gz", index=False)

print("save amy_tsne.csv.gz i hipp_tsne.csv.gz")