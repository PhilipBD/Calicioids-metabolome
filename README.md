In this repository are all the .csv files required to run the R script.

If you are interested in accessing raw data, please refer to:

    - MassIVE repository (MSV000097280) for metabolomics
    - GenBank (BioProject: PRJNA1219214) for microbiota data.
    - GNPS outputs https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=2e0f4a3abaae4e6b832722725a9b69a4

**Calicioids_metabolome.R** = R script from raw files to data analyses

**bacteria_abundance.csv** = Normalized bacteria abundance values per sample used to compile richness and diversity

**bacteria_presence-absence.csv** = Presence-absence data for bacteria per sample used to compile accumulation curves

**bacteria_taxa.csv** = Taxonomic information corresponding to each bacterial ASV

**fungi_abundance.csv** = Normalized fungi abundance values per sample used to compile richness and diversity

**fungi_presence-absence.csv** = Presence-absence data for fungi per sample used to compile accumulation curves

**fungi_taxa.csv** = Taxonomic information corresponding to each fungal ASV

**metabolome_abundance.csv** = Normalized fungi abundance values per sample used to compile richness and diversity

**metabolome_phyloseq_import.csv** = This is basically the MzMine2 quantification table formatted for uploading into a phyloseq object

**metabolome_presence-absence.csv** = Presence-absence data for metabolites per sample used to compile accumulation curves

**metabolome_taxa.csv** = Summary of annotations for inclusion as "taxa" in the phyloseq object

**metabolome_taxa_backup.csv** = Same as "metabolome_taxa.csv" but with more detailed annotations

**richness_calicioids.csv** = richness of ASVs and metabolites per sample used to build Figure 3

**shannon_calicioids.csv** = Shannon diversity index for ASVs and metabolites per sample used to build Figures 6 and S2

**taxa_16S_voronoi_C.csv** = Bacteria annotations used to build the voronoi treemap for axenic cultures (Figure 2d)

**taxa_16S_voronoi_F.csv** = Bacteria annotations used to build the voronoi treemap for non-lichenized fungi (Figure 2e)

**taxa_16S_voronoi_L.csv** = Bacteria annotations used to build the voronoi treemap for lichens (Figure 2f)

**taxa_ITS_voronoi.csv** = Fungi annotations used to build all fungal voronoi treemaps (Figure 2a-c)

**Venn_diagram_16S.csv** = Data used to build Venn diagram for bacteria (Figure S3b)

**Venn_diagram_ITS.csv** = Data used to build Venn diagram for fungi (Figure S3a)

**Venn_diagram_MET.csv** = Data used to build Venn diagram for metabolites (Figure S3c)

**Lichen_substances_piechart_13mars2025.cys** = Cytoscape file to visualize molecular network presented in Figure 4






