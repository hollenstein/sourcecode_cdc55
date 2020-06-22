import os
import sys
import itertools
import numpy as np
import pandas as pd


def protein_phospho_islands(ev, protein, quantified=True, prob_cutoff=0.7):
    """ Define protein phosphorylation islands.

    A phosphorylation islands contains all phosphorylation sites that are
    connected by overlapping peptides. Each island is defined by its first
    (entry 'lo') and last (entry 'hi')  phosphorylation site, as well as the
    protein amino acid positions (entries 'start' and 'end') covered by the
    corresponding peptides.

    :param ev: evidence table (pandas.DataFrame)
    :param protein: protein standard name
    :param quantified: bool, if True only use quantified evidence entries
    :param prob_cutoff: 0 to 1, use only phospho peptides with a 
        'Combined Phospho Probability' greater or equal to this value
    :returns: {island_id: {'lo': int, 'hi': int, 'start': int, 'end': int}}
    """
    def _split_to_int_tuples(list_of_strings, split_by=','):
        return [tuple(map(int, s.split(split_by))) for s in list_of_strings]

    m = np.all([
        ev['Protein Standard Name'] == protein,
        ev['Phospho (STY)'] > 0,
        ev['Combined Phospho Probability'] >= prob_cutoff,
    ], axis=0)
    if quantified:
        m = m & ev['Valid Quantification']

    phospho_positions = _split_to_int_tuples(ev['Phospho positions'][m])
    start_sequences = ev['Sequence start'][m]
    end_sequences = ev['Sequence end'][m]

    phospho_islands = _iterative_phospho_islands(
        phospho_positions, start_sequences, end_sequences
    )
    return phospho_islands


def find_island_id(phospho_islands, phospho_site):
    """ Find the phospho island containing the specified phospho site. 

    :returns: the phospho island ID of the matching phospho island
    """
    matching_island_id = None
    for island_id, island in phospho_islands.items():
        if phospho_site <= island['hi'] and phospho_site >= island['lo']:
            matching_island_id = island_id
            break
    return matching_island_id


def _iterative_phospho_islands(phospho, start, end):
    """ Iteratively refine phospho islands until the minimal number of 
    unconnected islands is found. """
    old_islands = {}
    phospho_islands = _extract_phospho_islands(phospho, start, end)
    while len(old_islands) != len(phospho_islands):
        old_islands = phospho_islands
        phospho_sequence_positions = []
        for island in phospho_islands.values():
            phospho_sequence_positions.append([
                (island['lo'], island['hi']), island['start'], island['end'],
            ])
        phospho, start, end = zip(*phospho_sequence_positions)
        phospho_islands = _extract_phospho_islands(phospho, start, end)
    return phospho_islands


def _extract_phospho_islands(phospho, start, end):
    """ phospho = [(3, 4),]; start = [1, ]; end = [5, ] """
    last_lo_site = 0
    last_hi_site = 0
    last_seq_start = 0
    last_seq_end = 0
    island_id = 0
    islands = {}
    phospho_sequence_positions = sorted(set(zip(phospho, start, end)))
    for phos_sites, seq_start, seq_end, in phospho_sequence_positions:
        lo_site = min(phos_sites)
        hi_site = max(phos_sites)
        if any([(seq_start <= last_hi_site and seq_end >= last_lo_site),
                (last_seq_start <= hi_site and last_seq_end >= lo_site)]):
            if lo_site < last_lo_site:
                last_lo_site = lo_site
                islands[island_id]['lo'] = last_lo_site
            if hi_site > last_hi_site:
                last_hi_site = hi_site
                islands[island_id]['hi'] = last_hi_site
            if seq_start < last_seq_start:
                last_seq_start = seq_start
                islands[island_id]['start'] = last_seq_start
            if seq_end > last_seq_end:
                last_seq_end = seq_end
                islands[island_id]['end'] = last_seq_end
        else:
            # Generate new Island
            island_id += 1
            last_lo_site = lo_site
            last_hi_site = hi_site
            last_seq_start = seq_start
            last_seq_end = seq_end
            islands[island_id] = {
                'lo': last_lo_site, 'hi': last_hi_site,
                'start': last_seq_start, 'end': last_seq_end
            }
    return islands


###########################################################
# IMPORT EVIDENCE TABLES AND GENERATE PHOSPHO ISLANDS #
###########################################################

# Load evidence tables to generate the phospho islands
root = '../'
table_dir = os.path.join(root, 'tables')
setup_names = [
    'cdc55D', 'cdc55-degron', 'SR', 'SR-hog1as',
    'rts1D', 'SR-cdc55D', 'SR-igoD', 'htb-gis1-rph1'
]

tables = []
for setup_name in setup_names:
    path = os.path.join(table_dir, 'evidence_' + setup_name + '.tsv')
    tables.append(pd.read_csv(path, sep='\t', index_col=False))
evidence_table = pd.concat(tables)


# Define phospho islands and calculate the unmodified counter peptide ratio
setups = ['15min', '30min', '60min']
protein_islands = {}
for protein in ['GIS1', 'RPH1']:
    phospho_islands = protein_phospho_islands(
        evidence_table, protein, quantified=True, prob_cutoff=0.7
    )

    # Extract counter peptides
    m = np.all([
        evidence_table['Protein Standard Name'] == protein,
        evidence_table['Valid Quantification'],
        evidence_table['Phospho (STY)'] == 0,
        np.any([(evidence_table['Setup'] == s) for s in setups], axis=0)
    ], axis=0)

    for island_id, phospho_island in phospho_islands.items():
        for setup in setups:
            phospho_islands[island_id][setup] = []
        for exp, group in evidence_table[m].groupby('Experiment'):
            group_mask = np.all([
                group['Sequence start'] <= phospho_island['hi'],
                group['Sequence end'] >= phospho_island['lo']
            ], axis=0)
            group_mean = group['Ratio Log2 normalized'][group_mask].mean()
            setup = group['Setup'].unique()[0]
            phospho_islands[island_id][setup].append(group_mean)
    protein_islands[protein] = phospho_islands


###############################################
# GENERATE PHOSPHORYLATION AND COUNTER TABLES #
###############################################

# Import phospho site tables
setup_names = ['cdc55D', 'cdc55-degron', 'SR', 'SR-igoD', 'htb-gis1-rph1']
tables = []
for setup_name in setup_names:
    path = os.path.join(root, 'tables', 'phospho_' + setup_name + '.tsv')
    tables.append(pd.read_csv(path, sep='\t', index_col=False))
phospho_table = pd.concat(tables)

protein_mask = np.any([
    (phospho_table['Protein Standard Name'] == 'GIS1'),
    (phospho_table['Protein Standard Name'] == 'RPH1'),
], axis=0)
phospho_table = phospho_table[protein_mask]


# Add phosphorylation site conservation annotation
def lookup_site(entry):
    aa, site = entry.split(':')
    return ''.join([site, '(', aa, ')'])
def lookup_pos(entry):
    aa, site = entry.split(':')
    return int(site)


gis1_rph1_conservation = [
    ('S:70', 'S:72'), ('T:152', 'T:175'), ('S:167', 'S:190'),
    ('S:212', 'S:243'), ('S:226', 'S:257'), ('S:262', 'S:293'),
    ('T:286', 'T:317'), ('S:304', 'S:335'), ('S:340', 'S:370'),
    ('S:398', 'S:405'), ('S:399', 'S:406'), ('S:400', 'S:407'),
    ('S:403', 'S:410'), ('T:404', 'T:411'), ('S:421', 'S:426'),
    ('S:424', 'S:429'), ('S:425', 'S:430'), ('S:429', 'S:434'),
    ('S:435', 'S:440'), ('T:441', 'T:446'), ('S:569', 'S:490'),
    ('S:576', 'S:497'), ('S:580', 'S:501'), ('S:590', 'S:511'),
    ('S:664', 'S:555'), ('S:674', 'S:568'), ('T:686', 'T:580'),
    ('S:690', 'S:584'), ('S:694', 'S:588'), ('S:696', 'S:590'),
    ('S:719', 'S:612'), ('T:720', 'T:613'), ('S:725', 'S:618'),
    ('S:772', 'S:664'), ('T:789', 'T:669'), ('S:812', 'S:693'),
    ('S:822', 'S:703'), ('S:838', 'S:719'), ('S:839', 'S:720'),
    ('T:844', 'T:725'), ('S:849', 'S:730'), ('S:852', 'S:733'),
    ('S:858', 'S:739')
]
conservation_lookup = {
    'GIS1': {lookup_pos(g): lookup_site(r) for g, r in gis1_rph1_conservation},
    'RPH1': {lookup_pos(r): lookup_site(g) for g, r in gis1_rph1_conservation}
}

conservation_annotation = []
for prot, sites in zip(
        phospho_table['Protein Standard Name'],
        phospho_table['Phospho positions']
    ):
    is_conserved = False
    annotation = []
    for site in map(int, sites.split(',')):
        if site in conservation_lookup[prot]:
            annotation.append(conservation_lookup[prot][site])
            is_conserved = True
        else:
            annotation.append('-')
    if is_conserved:
        conservation_annotation.append(' / '.join(annotation))
    else:
        conservation_annotation.append('')
phospho_table['Phospho site conservation'] = conservation_annotation


# Create a counter group table
setups = ['15min', '30min', '60min']
counter_group_table = []
counter_group_headers = [
    'Protein phospho sites', 'Phospho island ID',
    'Ratio Log2 normalized', 'Setup'
]
for protein_site, group in phospho_table.groupby('Protein phospho sites'):
    protein_name = group['Protein Standard Name'].get_values()[0]
    phospho_positions = group['Phospho positions'].get_values()[0]
    first_site = int(phospho_positions.split(',')[0])

    island_id = find_island_id(protein_islands[protein_name], first_site)
    phospho_island = protein_islands[protein_name][island_id]
    for setup in setups:
        setup_name = ' '.join([setup, '(CG)'])
        for ratio in phospho_island[setup]:
            counter_group_table.append(
                [protein_site, str(island_id), ratio, setup_name]
            )
counter_group_table = pd.DataFrame(
    counter_group_table, columns=counter_group_headers
)


# Define missing label columns for the phospho and counter group tables 
phospho_labels = phospho_table.groupby('Protein phospho sites').first()
phospho_labels.drop(
    ['Setup', 'Experiment', 'Ratio Log2 normalized', 'Phospho positions'],
    axis=1, inplace=True
)
counter_labels = counter_group_table.groupby('Protein phospho sites').first()
counter_labels.drop(['Setup', 'Ratio Log2 normalized'], axis=1, inplace=True)


# Add missing labels to the phospho_table and the counter_group_table
phospho_table = phospho_table.join(
    counter_labels, on='Protein phospho sites', how='left'
)
counter_group_table = counter_group_table.join(
    phospho_labels, on='Protein phospho sites', how='left'
)


###########################################
# GENERATE THE PHOSPHORYLATION SITE TABLE #
###########################################

combined_table = pd.concat([phospho_table, counter_group_table])
final_table = pd.pivot_table(
    data=combined_table,
    columns='Setup', values='Ratio Log2 normalized',
    aggfunc={'Ratio Log2 normalized': ['mean', 'count', 'std']},
    index=[
        'Protein Standard Name', 'Protein Systematic Name',
        'Phospho sites', 'Phospho island ID', 'Phospho (STY)',
        'Phospho site conservation'
    ],
)
final_table.columns = final_table.columns.swaplevel(0, 1)

setup_order = [
    '15min', '30min', '60min', '15min (CG)', '30min (CG)', '60min (CG)',
    'cdc55D', 'cdc55-degron', 'SR-igoD', 'SR'
]


# Replace missing and empty values
mean_fillna = '-'
count_fillna = 0
std_fillna = '-'
for setup in setup_order:
    final_table[(setup, 'mean')].fillna(mean_fillna, inplace=True)
    final_table[(setup, 'count')].fillna(count_fillna, inplace=True)
    final_table[(setup, 'std')].fillna(std_fillna, inplace=True)


# Reorganize column order
col_order = itertools.chain(
    *[[(s, 'mean'), (s, 'count'), (s, 'std')] for s in setup_order]
)
final_table = final_table.ix[:, col_order]
final_table.rename(
    columns={'mean': 'SILAC Ratio', 'count': 'Num', 'std': 'Std'}, inplace=True
)
final_table.reset_index(inplace=True)


# Export table
path = os.path.join(root, 'tables', 'Supplementary_table_3_CDC55.tsv')
final_table.to_csv(path, sep='\t', index=False)


"""
# Some test cases for defining phospho islands #

phospho_sequence_positions = [
    (1, 5, (3,)),
    (4, 10, (5, 10)),
    (11, 15, (15,)),
    (14, 20, (20,)),
]
start, end, phospho = zip(*phospho_sequence_positions)
phospho_islands = _iterative_phospho_islands(phospho, start, end)
# phospho_islands = {
#  1: {'hi': 10, 'lo': 3},
#  2: {'hi': 15, 'lo': 20},
# }

phospho_sequence_positions = [
    (1, 5, (3,)),
    (4, 10, (10,)),
    (5, 10, (5,)),
]
start, end, phospho = zip(*phospho_sequence_positions)
phospho_islands = _iterative_phospho_islands(phospho, start, end)
# phospho_islands = {
#  1: {'hi': 10, 'lo': 3},
# }

start = [1, 5, 15, 17, 20, 22]
end = [7, 12, 19, 19, 24, 26]
phospho = [(3,), (10,), (16, 18), (18,), (21, 23), (25,)]
phospho_islands = _iterative_phospho_islands(phospho, start, end)
# phospho_islands = {
#  1: {'hi': 3, 'lo': 3},
#  2: {'hi': 10, 'lo': 10},
#  3: {'hi': 18, 'lo': 16},
#  4: {'hi': 25, 'lo': 21}
# }

phospho_sequence_positions = [
    (1, 6, (3,)),
    (3, 8, (8,)),
    (3, 10, (5,)),

    (9, 13, (11, 12)),
    (9, 15, (15,)),

    (16, 19, (17,)),
    (16, 19, (19,)),
    (16, 19, (18,)),

    (20, 23, (23,)),
    (23, 26, (26,)),
    (27, 29, (27,)),
    (26, 28, (28,)),

    (30, 30, (30,))
]
start, end, phospho = zip(*phospho_sequence_positions)
phospho_islands = _iterative_phospho_islands(phospho, start, end)
# islands = {
#  1: {'hi': 8, 'lo': 3},
#  2: {'hi': 15, 'lo': 11},
#  3: {'hi': 19, 'lo': 17},
#  4: {'hi': 28, 'lo': 23}
#  5: {'hi': 30, 'lo': 30}
# }
"""
