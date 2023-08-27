# Bioinformatics
## T Cell Differentiation Process Analysis with R Programming

**Immunology:** DP -> DN -> SP -> B, T Cells

### My Role
- Conducted pairwise comparisons in in-vivo data
- Drew mRNA expression correlation map
- Processed STAR-aligned data with DESEQ2 and GSEA
- Counted significant numbers for each step
- Identified Ribosome, Oxidative Phosphorylation, and Cell Cycle as common pathways for all stages
- Utilized Cytoscape for path enrichment network analysis to explore pathway interactions

### Pathway Insights
- Nodes: KEGG Gene Sets
- Edges: Common number of gene sets
- Color-coded nodes based on significant adjusted p-values:
  - Green: Immune cells (e.g., Allograft Rejection, IGA Production, Cytokine-Cytokine Receptor Interaction, Jackstat Signaling Pathway)
  - Orange: Receptor and Immune Cell Signaling
  - Blue: DNA, RNA, Protein Synthesis or Regulation related to Cell Development


![network_invivo_final](https://user-images.githubusercontent.com/62280486/128388215-b041596a-6510-4a61-ae37-f48887dc6bea.png)
