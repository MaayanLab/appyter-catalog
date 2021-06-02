# Constants used to convert CELLECT output to pandas df and correct pvalue
CELLECT_OUTDIR = 'out/'
METHODS = ['H-MAGMA','LDSC', 'MAGMA']
PVAL_CORRECTION = 'bonferroni'
GWAS_GROUP_DICT = {
    'Structural MRI':['volume','thickness','area'],
    'DTI Tracts':['DTI'],
    'rs-fMRI Network':['Net'],
    "Schizophrenia": ['SCZ_PGC3_2020'],
    "Alzheimer's disease": ['AD_JANSENS2019'],
    "Bipolar": ['BIP_PGC3'],
    "Depression": ['PGC_depression2019'],
    "Autism spectrum disorder": ['ASD_2019'],
} # GWAS group name : [(regex) keywords]

# GSEA Constants
GENE_SET_LIST = [
    'ARCHS4_Tissues',
    'Aging_Perturbations_from_GEO_down',
    'Aging_Perturbations_from_GEO_up',
    'Allen_Brain_Atlas_down',
    'Allen_Brain_Atlas_up',
    'ClinVar_2019',
    'Disease_Signatures_from_GEO_down_2014',
    'Disease_Signatures_from_GEO_up_2014',
    'ENCODE_Histone_Modifications_2013',
    'ENCODE_Histone_Modifications_2015',
    'Epigenomics_Roadmap_HM_ChIP-seq',
    'GO_Biological_Process_2018',
    'GO_Cellular_Component_2018',
    'GO_Molecular_Function_2018',
    'GTEx_Tissue_Sample_Gene_Expression_Profiles_down',
    'GTEx_Tissue_Sample_Gene_Expression_Profiles_up',
    'GWAS_Catalog_2019',
    'HMDB_Metabolites',
    'Human_Gene_Atlas',
    'Jensen_COMPARTMENTS',
    'Jensen_DISEASES',
    'Jensen_TISSUES',
    'KEGG_2019_Human',
    'Reactome_2016',
    'Tissue_Protein_Expression_from_ProteomicsDB',
    'WikiPathways_2019_Human',
    'miRTarBase_2017',
                ]
TOP_ANNOT = 5 # analyses TOP_ANNOT number of cell-types ranking the highest per GWAS group
TOP_FREQ = 0.01 # only takes genes within the TOP_FREQ ES value for each cell-type
OVERWRITE_GSEA_ANALYSIS = False # overwrites old gsea files if True
