from typing import List

import os
import numpy as np
import pandas as pd
import seaborn as sns

class PathwayDatabase:
    DATABASES = {
        'hallmark': ['h.all.v2023.2.Hs.symbols.gmt', 'hallmark.csv'],
        'canonical pathways': ['c2.cp.v2023.2.Hs.symbols.gmt', 'cp.csv'],
        'reactome': ['c2.cp.reactome.v2023.2.Hs.symbols.gmt', 'reactome.csv'],
        'transcription factor targets': ['c3.tft.v2023.2.Hs.symbols.gmt', 'tf.csv'],
        'GOBP': ['c5.go.bp.v2023.2.Hs.symbols.gmt', 'go_bp.csv'],
        'GOCC': ['c5.go.cc.v2023.2.Hs.symbols.gmt', 'go_cc.csv'],
        'GOMF': ['c5.go.mf.v2023.2.Hs.symbols.gmt', 'go_mf.csv']
    }

    def __init__(self, organ: str, database: str, pathway: str=None, genes: List[str]=None) -> None:
        self.pathways = {}
        self.organ = organ
        
        if database != 'custom':
            with open(f'./pathway_analysis/data/databases/{self.DATABASES[database][0]}') as file:
                for line in file:
                    pathway, _, *genes = line.removesuffix('\n').split('\t')
                    self.pathways[pathway] = genes
        else:
            self.pathways[pathway] = genes

        self.data = pd.read_csv(f'./pathway_analysis/data/{self.organ}/{self.DATABASES[database][1]}', index_col=0)

    def __getitem__(self, pathway: str) -> List[str]:
        if pathway in self.pathways:
            return self.pathways[pathway]
    
    def search(self, query: str) -> List[str]:
        return [pathway for pathway in self.pathways if query in pathway]

    def get_clusters(self, clusters: str, highlight: str) -> pd.DataFrame:
        files = os.listdir(f'./pathway_analysis/data/{self.organ}/diff_genes/')

        raw_diff_genes = []
        for cluster in files:
            if highlight == 'z-score':
                raw_diff_genes.append(pd.read_csv(f'./pathway_analysis/data/{self.organ}/diff_genes/{cluster}', 
                                                index_col=0).iloc[:, 0].rename(f'{cluster.split("_")[0]}'))
            else:
                raw_diff_genes.append(pd.read_csv(f'./pathway_analysis/data/{self.organ}/diff_genes/{cluster}', 
                                                index_col=0).iloc[:, 1].rename(f'{cluster.split("_")[0]}'))

        self.diff_genes = pd.concat(raw_diff_genes, axis=1, join='inner')

        if clusters != 'all':
            clusters = clusters.split()
            self.diff_genes = self.diff_genes[clusters]

        self.diff_genes.columns = self.diff_genes.columns.astype('int32')
        self.diff_genes = self.diff_genes.reindex(sorted(self.diff_genes.columns), axis=1)

        return self.diff_genes


def heatmap(data: pd.DataFrame, method: str, metric: str):
    n_rows, n_cols = data.shape
    # data = data.applymap(np.nan_to_num)    

    return sns.clustermap(data=data,
                          method=method,
                          metric=metric,
                          figsize=(n_cols * 2, n_rows * 0.5),
                          cmap='bwr',
                          center=0,
                          fmt='',
                          col_cluster=False)
