import numpy as np
import pandas as pd

def get_all_genesets():
    """
    """
    # all genes
    f = "/u/home/f/f7xiesnm/project-zipursky/v1-bb/v1/data/merfish/merfish_genes.txt" 
    genes = np.loadtxt(f, dtype='str')
    
    # ABC genes
    f = '/u/home/f/f7xiesnm/project-zipursky/v1-bb/v1/results/MERFISH_gene_panel_Version1_March9.csv'
    df1 = pd.read_csv(f)
    f = '/u/home/f/f7xiesnm/project-zipursky/v1-bb/v1/results/gene_ptime_P28_L23_Mar27.tsv'
    df2 = pd.read_csv(f).sort_values('gene_ptime')
    df = pd.merge(df1, df2, left_on='gene_name_data', right_on='gene').sort_values('gene_ptime')

    genes_abco = df['gene_name_vizgen'].values
    agenes = df[df['P17on_x']=='A']['gene_name_vizgen'].values
    bgenes = df[df['P17on_x']=='B']['gene_name_vizgen'].values
    cgenes = df[df['P17on_x']=='C']['gene_name_vizgen'].values

    
    iegs = np.array([ # this list focuses on L2/3 ERGs # n=16
        'Arc',
        'Crem',        
        'Egr1',
        'Egr2',        
        'Egr4',
        'Fos',        
        'Fosb',
        'Fosl2',
        'Junb',
        'Npas4',        
        'Nr4a1',
        'Nr4a2',
        'Nr4a3',
        'Per1',
        
        'Tiparp',
        'Tnfaip6',

        # 'Cebpb',
        # 'Klf4',
        # 'Atf3',
        # 'Klf2',
        # 'Maff',
        # 'Klf10',
        # 'Jund',
        # 'Mafk',
        # 'Srf',
    ])
    
    iegs_old = np.array([ # this list focuses on both neuron and non-neuronal ERGs - excluding some non-canonical L2/3 IEGs # n=22?
        'Fos',
        'Fosl2',
        'Nr4a1',
        'Egr1',
        'Nr4a2',
        'Nr4a3',
        'Crem',
        'Egr4',
        'Fosb',
        'Egr2',
        'Per1',
        'Junb',
        'Cebpb',
        'Klf4',
        'Atf3',
        'Klf2',
        'Maff',
        'Klf10',
        'Jund',
        'Mafk',
        'Srf',
        'Arc',
    ])

    easifish_genes = np.array([
        'Slc17a7', 'Gria3', 'Rorb', 'Kcnq5',

        'Adamts2', 'Sorcs3', 'Chrm2', 
        'Rfx3', 
        'Cdh13', 'Cdh12', 
        'Epha10',
        'Trpc6',


        'Kcnip3',  'Cntn5',
        'Cntnap2', 'Gabrg3',
        'Kcnh5',  'Grm8', 
        'Ncam2',  'Baz1a',
    ])

    up_agenes = np.array([
              'Vwc2l',
              'Tmem117', 
              'Kcnk13', 
              'Airn',
              'Ttc28', 
              'Col23a1',
              'Cdh13', 
              'Glis3', 
              'Iqgap2',
              'Nckap5',
              'Cpne6',
              '1700086L19Rik',
              'Otof',
              'Pcdh19',
              'Sema6a',
              'Gm42722',
              'Gabrg3',
              'Grm8', 
              'Palm2', 
              'Syt17',
             ])
    up_agenes = np.intersect1d(up_agenes, genes_abco)
    
    genesets = {
        'allmerfish': genes,
        'a': agenes, 
        'b': bgenes,
        'c': cgenes,
        'i': iegs,
        'a_up': up_agenes,
    }
    return genesets, df


def rename_genes(genes):
    """
    """
    rename_dict = {
        'Fam19a1': 'Tafa1', # snRNA-seq to MERFISH
    }
    
    genes_renamed = genes.copy()
    for i, g in enumerate(genes):
        if g in rename_dict.keys():
            genes_renamed[i] = rename_dict[g]
            
    return genes_renamed 