# Deconvolution of single-worm omics data
This project explores deconvolution of bulk omics data (especially single-worm proteomics) to resolve cellular and tissue composition in C. elegans at single-organism resolution. By integrating transcriptomic and proteomic datasets, we demonstrate how inferred cell-type proportions can serve as proxy readouts for aging-associated dynamics, opening new avenues for high-resolution organismal (populational) study.

- [Analysis Code for the manuscript](notebooks/deconvolution_for_pub.ipynb)
- [Colab](https://colab.research.google.com/github/iron-lion/deconvolution_C_elegans/blob/main/misc/colab_deconvolution.ipynb) : Deconvolution with simple linear SVR

Cell type specific expression profiles of C. elegans [data](https://github.com/iron-lion/deconvolution_C_elegans/blob/main/data/top_100_markers_postsub.csv) from [Ghaddar et al.](https://doi.org/10.1126/sciadv.adg0506) (obtained from [WormSeq](https://wormseq.org/)).

## Deconvolution with single-worm RNA sequencing data
Applying deconvolution to single-worm RNA-sequencing data demonstrates that deconvolution-derived cellular and tissue proportions can serve as robust readouts of dynamic changes in cells and tissues during aging.

The public single-worm RNA sequencing from [Eder et al.](https://doi.org/10.1016/j.cell.2024.05.050)

[All results](results/Eder/) with different parameters

## Deconvolution with single-worm proteomics
Deconvolution can show proteome-level dynamics as well.

The public single-worm proteomics from [Zhu et al.](https://doi.org/10.1111/acel.14055)

[All results](results/Zhu/) with different parameters

## Comparison with ground truth
[Summary results](results/ground_truth/) with different parameters

[Individual tissue results](results/Tissue_proportions/) with different parameters

# Cite
Park Y, Antebi A. 2026. Deconvolution of Single-Organism Omics Resolves Cellular and Tissue Dynamics during Aging. microPublication Biology. [10.17912/micropub.biology.001989](https://doi.org/10.17912/micropub.biology.001989).
