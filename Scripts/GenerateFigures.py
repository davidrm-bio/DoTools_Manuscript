# Import libraries
from pathlib import Path

import scanpy as sc
import dotools_py as do
from dotools_py.pp._importer import _qc_vln
import celltypist
import anndata as ad
import pandas as pd

import matplotlib.pyplot as plt

sc.set_figure_params(scanpy=False, vector_friendly=True)
do.settings.session_settings()

data_path = Path('/Volumes/Zmm_shared/David/Papers/dotools/Objects')
figure_path = Path('/Volumes/Zmm_shared/David/Papers/dotools/Figures')

adata = ad.read_h5ad(data_path / 'adata.h5ad')

########################################################################################################################
# Figure 1 - Workflow DOTools
########################################################################################################################

# Figure 1A - ViolinPlots of QC metrics
_qc_vln(adata, title='', path=figure_path, filename='Figure1A_ViolinPlot_QCMetrics.svg')

# Figure 1B - UMAP splitting by condition
do.pl.umap(adata, 'leiden', split_by='condition', labels='leiden', path=figure_path,
           filename='Figure1B_UMAP_SplittingCondition.svg', figsize=(6, 4), size=15)


# Figure 1B - Celltypist Probability Dotplot
model_loaded = celltypist.models.Model.load(model="Healthy_COVID19_PBMC.pkl")
adata_copy = adata.copy()
adata_copy.X = adata_copy.X.toarray()
pred = celltypist.annotate(adata_copy, model=model_loaded, majority_voting=False, over_clustering="leiden")
spred = pred.to_adata()

# Only visualise of a subset of cell type
spred = spred[spred.obs.predicted_labels.isin(["B_naive", "B_immature", "CD14_mono", "CD16_mono", "CD4.Naive",
                                               "CD4.CM", "CD8.TE", "NK_16hi", "pDC"])]

# Correct the object
pred.adata = spred
pred.decision_matrix = pred.decision_matrix[pred.decision_matrix.index.isin(list(spred.obs_names))]
pred.decision_matrix = pred.decision_matrix[[ "B_naive", "B_immature", "CD14_mono", "CD16_mono", "CD4.Naive", "CD4.CM", "CD8.TE", "NK_16hi", "pDC"]]

pred.predicted_labels = pred.predicted_labels[pred.predicted_labels.index.isin(list(spred.obs_names))]
pred.probability_matrix = pred.probability_matrix[pred.probability_matrix.index.isin(list(spred.obs_names))]
pred.probability_matrix = pred.probability_matrix[[ "B_naive", "B_immature", "CD14_mono", "CD16_mono", "CD4.Naive", "CD4.CM", "CD8.TE", "NK_16hi", "pDC"]]
pred.cell_count = 1960
pred.predicted_labels['predicted_labels'] = pd.Categorical(pred.predicted_labels['predicted_labels'].astype(str))

fig, axs = plt.subplots(1, 1, figsize=(5, 4))
axs = celltypist.dotplot(pred, use_as_prediction="predicted_labels", use_as_reference="leiden", title="",
                         show=False, ax=axs, swap_axes=True)
axs["mainplot_ax"].spines[["top", "right"]].set_visible(True)
axs["mainplot_ax"].set_xticklabels(axs["mainplot_ax"].get_xticklabels(), fontweight='bold', rotation=45, ha='right', va='top')
plt.savefig( figure_path / 'Figure1B_DotplotCelltypist.svg', bbox_inches="tight")


# Figure 1B - Barplot cell proportions
do.pl.cell_props(adata, 'annotation_recluster', 'condition', 'batch',
                 transform='arcsin', figsize=(3.5, 4), cond_order=['healthy', 'disease'], legend_fontsize=10,
                 title="Annotation", path=figure_path, filename="Figure1B_CellProportions.svg")

# Figure 1C - 3D Dotplot
do.pl.dotplot(adata, 'condition', 'NKG7', 'annotation_recluster', add_stats='x_axis',
              figsize=(3, 4), path=figure_path, filename='Figure1C_3Dotplot_AnnotationCondition.svg',
              )
do.pl.dotplot(adata, 'condition', 'NKG7', 'annotation_recluster', add_stats='x_axis',
              figsize=(3, 4), path=figure_path, filename='Figure1C_3Dotplot_AnnotationCondition.pdf',
              )
# Figure 1C - Heatmap
do.pl.heatmap(adata, 'annotation_recluster', ['CD14', 'CD3D', 'CD79A', 'CD8A', 'NKG7'],
              xticks_rotation=45, figsize=(4, 6), add_stats=True, cluster_x_axis=True, cluster_y_axis=True,
              path=figure_path, filename='Figure1C_Heatmap.svg', stats_x_size=18,
             )

# Figure 1C - Barplot NKG7
do.pl.barplot(adata, 'condition', 'NKG7', ctrl_cond='healthy', groups_cond=['disease'],
              palette={'healthy': 'sandybrown', 'disease':'royalblue'}, path=figure_path,
              filename='Figure1C_BarplotNKG7.svg')

# Figure 1C - Split bar GSEA
do.tl.rank_genes_groups(adata, groupby='condition')
table = do.get.dge_results(adata)
df = do.tl.go_analysis(table, gene_key='GeneName', pval_key='padj', log2fc_key='log2fc')
df = df[df['Adjusted P-value'] < 0.05]
do.pl.split_bar_gsea(df, 'Term',  'Combined Score', 'state', 'enriched',
                     figsize=(5, 4), path=figure_path, filename='Figure1C_SplitGSEA.svg')

# Figure 1C - Split Boxplot RPL11
do.utility.save_rds('/Volumes/Zmm_shared/David/Papers/dotools/Objects/adata.rds', 'batch', adata)

do.pl.boxplot(adata, "annotation", "RPL11", hue="condition", figsize=(6, 4),
              palette={"healthy":'sandybrown', "disease":"royalblue"}, txt_size=10, txt="",
              ctrl_cond='healthy', groups_cond=['disease'], hue_order=['healthy', 'disease'],
              path=figure_path, filename='Figure1C_Boxplot_GeneAnnotationHueCondition.svg')


# Figure 2 - Split bar GSEA readable
do.tl.rank_genes_groups(adata, groupby='condition')
table = do.get.dge_results(adata)
df = do.tl.go_analysis(table, gene_key='GeneName', pval_key='padj', log2fc_key='log2fc')
df = df[df['Adjusted P-value'] < 0.05]
df["Term"] = df["Term"].str.split(" \(GO").str[0]
do.pl.split_bar_gsea(df, 'Term',  'Combined Score', 'state', 'enriched',
                     figsize=(8, 6), path=figure_path, filename='Figure2G_SplitGSEA.svg')


# Figure 2C - ViolinPlot
do.pl.violin(adata, 'condition', 'NKG7', ctrl_cond='healthy', groups_cond=['disease'],
              palette={'healthy': 'sandybrown', 'disease':'royalblue'}, path=figure_path,
              filename='Figure2C_ViolinNKG7.svg', scatter=True, figsize=(3, 4.2))

