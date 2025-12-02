import pandas as pd
import numpy as np


def load_convert_dict(filepath='../data/c_elegans.PRJNA13758.WS294.geneIDs.txt'):
    ws294 = pd.read_csv(filepath, header=None, index_col=None)
    convert_dict = dict()
    for _, row in ws294[[1, 2, 3, 4, 5]].iterrows():
        if row[4] == 'Dead':
            continue
        if row[5] != 'protein_coding_gene':
            continue

        if type(row[2]) is str:
            convert_dict[row[2]] = row[1]
        if type(row[3]) is str:
            convert_dict[row[3]] = row[1]

    return convert_dict


def zhu_ref_data(filepath='../data/Zhu_et_al.txt'):
    target_df = pd.read_csv(filepath, header=0, index_col=5, sep='\t')

    # table restructure
    genename_rows = [True if type(x) is str and 'GN=' in x else False for x in target_df.index]
    expression_columns = [True if 'LFQ' in x else False for x in target_df.columns]
    target_df = target_df.loc[genename_rows, expression_columns]
    target_df.index = [x.split('GN=')[1].split(' ')[0] for x in target_df.index]

    convert_dict = load_convert_dict()
    # NOTE: RAW signal log transformation
    target_df.index = [convert_dict[x.replace('CELE_', '')] if x in convert_dict.keys() else 'NA' for x in target_df.index]
    target_df = np.log(target_df + 1.0)
    target_df = target_df[target_df.index != 'NA']
    target_df = target_df[~target_df.index.duplicated(keep='first')]
    target_df = target_df.sort_index()

    # label extraction
    labels = [int(x.split(' ')[2].split('_')[0][1:]) for x in target_df.columns]

    target_df = target_df.T

    return target_df, labels


def load_cell_to_group(filepath='../data/abbas_cell_markers.csv'):
    markers = pd.read_csv(filepath, index_col=None, header=0)

    marker_genes_convert = {}
    last_type = ''
    for i, row in markers.iterrows():
        if type(row['Final annotation']) is str and len(row['Final annotation']) > 0:
            group = row['CeNGEN_annotation']
            last_type = row['Final annotation']
            marker_genes_convert[last_type] = group

    return marker_genes_convert
