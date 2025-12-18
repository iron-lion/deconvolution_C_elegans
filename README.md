# Deconvolution of single-worm omics data
Inferring tissue or cell-type proportions from transcriptomic and proteomic data at single-worm resolution.

<img src="https://github.com/iron-lion/deconvolution_C_elegans/blob/main/misc/Deconvolution.png" width=50% height=50%>

- [Main Analysis Code](notebooks/deconvolution_for_pub.ipynb)
- [Colab](https://colab.research.google.com/github/iron-lion/deconvolution_C_elegans/blob/main/colab_deconvolution.ipynb) : Deconvolution with simple linear SVR

Cell type specific expression profiles of C. elegans [data](https://github.com/iron-lion/deconvolution_C_elegans/blob/main/data/top_100_markers_postsub.csv) from [Ghaddar et al.](https://doi.org/10.1126/sciadv.adg0506) (obtained from [WormSeq](https://wormseq.org/)).

## Deconvolution with single-worm RNA sequencing data
Applying deconvolution to single-worm RNA-sequencing data demonstrates that deconvolution-derived cellular and tissue proportions can serve as robust readouts of dynamic changes in cells and tissues during aging.

The public single-worm RNA sequencing from [Eder et al.](https://doi.org/10.1016/j.cell.2024.05.050)

#### minimum number of marker genes = 0 & minimum marker gene score = 0.8
<img src="https://github.com/iron-lion/deconvolution_C_elegans/blob/main/results/swRNA_N2_proportion_minnum_0_minscore_0.08.svg" width=75%>

## Deconvolution with single-worm proteomics
Deconvolution can show proteome-level dynamics as well.

The public single-worm proteomics from [Zhu et al.](https://doi.org/10.1111/acel.14055)

![Alt text](https://github.com/iron-lion/deconvolution_C_elegans/blob/main/results/Zhu_proportion_minnum_5_minscore_0.08.svg)

## Deconvolution results with different parameters (swRNAseq)
#### minimum number of marker genes = 0 & minimum marker gene score = 0
<img src="https://github.com/iron-lion/deconvolution_C_elegans/blob/main/results/swRNA_N2_proportion_minnum_0_minscore_0.0.svg" width=75%>

#### minimum number of marker genes = 0 & minimum marker gene score = 1.0
<img src="https://github.com/iron-lion/deconvolution_C_elegans/blob/main/results/swRNA_N2_proportion_minnum_0_minscore_0.1.svg" width=75%>

#### minimum number of marker genes = 15 & minimum marker gene score = 1.0
<img src="https://github.com/iron-lion/deconvolution_C_elegans/blob/main/results/swRNA_N2_proportion_minnum_15_minscore_0.1.svg" width=75%>
