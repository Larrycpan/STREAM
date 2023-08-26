# STREAM

STREAM (**S**ingle-cell enhancer regula**T**ory netwo**R**k inference from gene **E**xpression **A**nd chro**M**atin accessibility) is a computational 
method specifically designed to infer enhancer-driven gene regulatory networks (eGRNs) by analyzing transcriptome and chromatin
accessibility jointly profiled from the same single cells. The method aims to improve the mapping accuracy of relationships between transcription factors (TFs),
enhancers, and target genes. The algorithm employs two main strategies. First, it uses a 
[Steiner Forest Problem (SFP) model](https://www.sciencedirect.com/science/article/pii/S1570866709000628) based on a heterogeneous graph to
identify enhancer-gene relations that are highly relevant in a specific context. Second, it incorporates a hybrid biclustering approach combined with 
[submodular optimization](https://link.springer.com/article/10.1007/s10626-019-00308-7) to determine eGRNs by selecting an optimal subset 
from sets of co-regulated genes and co-accessible enhancers. 
Based on these strategies, STREAM operates through an iterative framework to find enhancer-driven regulons (eRegulons) based on both transcriptome 
and chromatin accessibility data. eRegulons specific to a cell type collectively compose an eGRN. 
It is based on several popular packages for single-cell transcriptome and chromatin accessibility analysis, particularly, 
[Seurat](https://satijalab.org/seurat/), [Signac](https://satijalab.org/signac/index.htm), [cicero](https://cole-trapnell-lab.github.io/cicero-release/docs_m3/), 
and [IRIS-FGM](https://bioconductor.org/packages/release/bioc/html/IRISFGM.html). 
Picture below showed the methodological details of STREAM for eGRN inference.

![](images/details.jpg)
