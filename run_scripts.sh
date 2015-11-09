#!/bin/bash
#
# Run the mutation and copy number data integration
# scripts in ./tools
#
# You will need to obtain mutation and copy number
# data sets from the original sources (COSMIC, EGA
# or CCLE GISTIC calls from cBioPortal
#

perl ./tools/process_CCLEgistic_data.pl \
--ccle_cna ./resources/CCLE_cancer_gene_gistic_calls_150121.txt \
--genes ./resources/gene_name_dictionary_v03.txt \
--cell_lines ./resources/cell_line_name_dictionary.txt \
--output_genes_list ./resources/cancer_gene_classifications_v0.4.txt \
--mut_consequence_classes ./resources/mutation_consequence_classifications.txt \
--output ./data_tables/processed_CCLE_cnv_data_150225



perl ./tools/process_WTSIcnv_data.pl \
--wtsi_cnv_data ~/Dropbox/Intercell_II_analysis/15.Variant_call_data/WTSI/Variant_Catalogue_Chris_lord.snp6_exonscreen_calls.txt \
--genes ./resources/gene_name_dictionary_v03.txt \
--cell_lines ./resources/cell_line_name_dictionary.txt \
--output_genes_list ./resources/cancer_gene_classifications_v0.4.txt \
--mut_consequence_classes ./resources/mutation_consequence_classifications.txt \
--output ./data_tables/processed_WTSI_cnv_data_150225


perl ./tools/process_WTSIexprn_data.pl \
--wtsi_expression_z_data ./resources/expression_copy/expression_zscores_within_tissues_141001.txt \
--genes ./resources/gene_name_dictionary_v03.txt \
--cell_lines ./resources/cell_line_name_dictionary.txt \
--output_genes_list ./resources/cancer_gene_classifications_v0.4.txt \
--mut_consequence_classes ./resources/mutation_consequence_classifications.txt \
--output ./data_tables/processed_WTSI_exprn_data_150225


perl ./tools/process_cosmic_data.pl \
--cosmic_muts ~/Documents/113_data_sources/COSMIC/CosmicCLP_CompleteExport_v68.tsv \
--mut_freqs ./resources/Davoli_mutation_data/Mutation_Dataset.txt.mut_counts.txt \
--genes ./resources/gene_name_dictionary_v03.txt \
--cell_lines ./resources/cell_line_name_dictionary.txt \
--output_genes_list ./resources/cancer_gene_classifications_v0.4.txt \
--mut_consequence_classes ./resources/mutation_consequence_classifications.txt \
--output ./data_tables/processed_cosmic_exome_150225


perl ./tools/process_WTSIexome_data.pl \
--wtsi_exome_data ~/Dropbox/Intercell_II_analysis/15.Variant_call_data/WTSI/Variant_Catalogue_Chris_lord.caveman_pindel_calls.txt \
--mut_freqs ./resources/Davoli_mutation_data/Mutation_Dataset.txt.mut_counts.txt \
--genes ./resources/gene_name_dictionary_v03.txt \
--cell_lines ./resources/cell_line_name_dictionary.txt \
--output_genes_list ./resources/cancer_gene_classifications_v0.4.txt \
--mut_consequence_classes ./resources/mutation_consequence_classifications.txt \
--output ./data_tables/processed_WTSI_exome_data_150225


perl ./tools/process_ovarian_exome_data.pl \
--ovarian_muts ~/Dropbox/Intercell_II_analysis/15.Variant_call_data/ovarian/ovarian_exome_variants_140521.txt \
--mut_freqs ./resources/Davoli_mutation_data/Mutation_Dataset.txt.mut_counts.txt \
--genes ./resources/gene_name_dictionary_v03.txt \
--cell_lines ./resources/cell_line_name_dictionary.txt \
--output_genes_list ./resources/cancer_gene_classifications_v0.4.txt \
--mut_consequence_classes ./resources/mutation_consequence_classifications.txt \
--output ./data_tables/processed_ovarian_exome_data_150225


perl ./tools/process_bone_exome_data.pl \
--bone_muts ~/Dropbox/Intercell_II_analysis/15.Variant_call_data/wtsi_bone/Osteosarcoma_cell_line_exomes_WTSI_141008.txt \
--mut_freqs ./resources/Davoli_mutation_data/Mutation_Dataset.txt.mut_counts.txt \
--genes ./resources/gene_name_dictionary_v03.txt \
--cell_lines ./resources/cell_line_name_dictionary.txt \
--output_genes_list ./resources/cancer_gene_classifications_v0.4.txt \
--mut_consequence_classes ./resources/mutation_consequence_classifications.txt \
--output ./data_tables/processed_bone_exome_data_150225


