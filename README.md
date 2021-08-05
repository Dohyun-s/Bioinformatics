# Bioinformatics
T cell differentiation process Analysis with R programming

Immunology: DP-> DN->SP ->B, T cells

My role: pairwise comparison in in-vivo Data,
->draw mRNA expression correlation map, 
run the preprocessed data from STAR alignment with DESEQ2, and GSEA.

->count the significant number for each step.
In result, Ribosome, oxidative phosphorylation, cell cycle are the most common for all stages.

->Then run Cytoscape path enrichment network in order to know about pathway interaction.

Node: KEGG gene set, Edge: common number of gene set 

color: significant adjust p-value

Green: immune cell such as allograft rejection, IGA production, cytokine-cytokine receptor interaction, jackstat signaling pathway

orange: receptor signaling and Immune cell signaling

Blue: DNA, RNA, protein synthesis or regulation related to cell development

