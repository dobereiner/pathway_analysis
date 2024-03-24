import pandas as pd
import numpy as np

LRDB = pd.read_csv('./pathway_analysis/ligand_receptor/human_lr_pair.txt', sep='\t')

LIGANDS = LRDB['ligand_gene_symbol'].unique()
RECEPTORS = LRDB['receptor_gene_symbol'].unique()


def annotate(df: pd.DataFrame, lr_choice: str) -> pd.DataFrame:    
    match lr_choice:
        case 'receptors':
            df = df.loc[df.index.isin(RECEPTORS)]
        case 'ligands':
            df = df.loc[df.index.isin(LIGANDS)]
        case 'both':
            genes = np.concatenate([LIGANDS, RECEPTORS])
            df = df.loc[df.index.isin(genes)]
            df.index = df.index.map(lambda gene: gene + ' (receptor)' 
                                    if gene in RECEPTORS else gene + ' (ligand)')
        case _:
            raise ValueError
    
    return df
