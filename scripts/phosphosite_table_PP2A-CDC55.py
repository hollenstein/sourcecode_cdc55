import os
import sys
import pandas as pd
import itertools


# Import phospho site tables
root = 'C:/Users/hollend5/mass_spectrometry/pp2a'
table_dir = os.path.join(root, 'tables')
setup_names = ['cdc55D', 'cdc55-degron', 'SR-cdc55D', 'SR-igoD',
               'SR', 'SR-hog1as', 'rts1D']
tables = []
for setup_name in setup_names:
    path = os.path.join(table_dir, 'phospho_' + setup_name + '.tsv')
    tables.append(pd.read_csv(path, sep='\t', index_col=False))
combined_table = pd.concat(tables)


# Add M-Track results annotations
mtrack_info = {
    '+': ['REG1', 'HOG1', 'NET1', 'IGO1', 'CHS5', 'NUP60',
          'STE5', 'EPO1', 'SHE3', 'NUP2', 'SEC16', 'GET2',
          'EAP1', 'KEL1', 'NAB2', 'YAP1', 'GIS1', 'RPH1'],
    '-': ['BLM10', 'HXK2', 'BFA1', 'RIM15', 'IGO2', 'RGC1',
          'ASK10', 'OXR1', 'SMY2'],
}
combined_table['M-Track results'] = ''
for mtrack_text, proteins in mtrack_info.items():
    for protein in proteins:
        m = (combined_table['Protein Standard Name'] == protein)
        combined_table.loc[m, 'M-Track results'] = mtrack_text


# Generate a pivote table of phosphorylation sites with mean, count, std 
pivoted_table = pd.pivot_table(
    data=combined_table, columns='Setup', values='Ratio Log2 normalized',
    aggfunc={'Ratio Log2 normalized': ['mean', 'count', 'std']},
    index=[
        'Protein Standard Name', 'Protein Systematic Name',
        'Phospho sites', 'Phospho (STY)', 'S/T-P motif', 'M-Track results'
    ],
)
pivoted_table.columns = pivoted_table.columns.swaplevel(0, 1)


# Replace missing and empty values
mean_fillna = '-'
count_fillna = 0
std_fillna = '-'
for setup in setup_names:
    pivoted_table[(setup, 'mean')].fillna(mean_fillna, inplace=True)
    pivoted_table[(setup, 'count')].fillna(count_fillna, inplace=True)
    pivoted_table[(setup, 'std')].fillna(std_fillna, inplace=True)


# Reorganize column order
col_order = itertools.chain(
    *[[(s, 'mean'), (s, 'count'), (s, 'std')] for s in setup_names]
)
pivoted_table = pivoted_table.ix[:, col_order]
pivoted_table.rename(
    columns={'mean': 'SILAC Ratio', 'count': 'Num', 'std': 'Std'}, inplace=True
)
pivoted_table.reset_index(inplace=True)


# Export table
path = os.path.join(root, 'tables', 'Supplementary_table_2_CDC55.tsv')
pivoted_table.to_csv(path, sep='\t', index=False)
