from typing import List

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
            with open(f'./pathway_analysis/databases/{self.DATABASES[database][0]}') as file:
                for line in file:
                    pathway, _, *genes = line.removesuffix('\n').split('\t')
                    self.pathways[pathway] = genes
        else:
            self.pathways[pathway] = genes

        self.data = pd.read_csv(f'./pathway_analysis/{self.organ}/{self.DATABASES[database][1]}', index_col=0)

    def __getitem__(self, pathway: str) -> List[str]:
        if pathway in self.pathways:
            return self.pathways[pathway]
    
    def search(self, query: str) -> List[str]:
        return [pathway for pathway in self.pathways if query in pathway]


def heatmap(data: pd.DataFrame, method: str, metric: str):
    n_rows, n_cols = data.shape

    return sns.clustermap(data=data,
                          method=method,
                          metric=metric,
                          figsize=(n_cols * 2, n_rows * 0.5),
                          cmap='bwr',
                          center=0,
                          fmt='',
                          col_cluster=False)
