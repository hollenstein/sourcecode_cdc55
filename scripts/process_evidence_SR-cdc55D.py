import os
import max_quant as MQ


# Process the MaxQuant results of the stress response cdc55 delta experiment
if __name__ == '__main__':
    root = '../'
    table_dir = os.path.join(root, 'tables')
    fasta_path = os.path.join(root, 'yeast_cont_20140324.fasta')
    mq_txt_folder = 'MQ_txt_SR-cdc55D'
    evidence_path = os.path.join(
        root, 'MaxQuant', mq_txt_folder, 'evidence.txt'
    )
    experiments_toreverse = ['05_cdc55_nacl', '06_cdc55_nacl']
    setup = 'SR-cdc55D'

    # Import Evidence
    evidence, normalization_data = MQ.process_evidence(evidence_path, fasta_path)
    evidence['Setup'] = setup
    MQ.reverse_ratios(evidence, experiments_toreverse)

    # Generate output files
    evidence_outpath = os.path.join(table_dir, 'evidence_' + setup + '.tsv')
    phosphosite_outpath = os.path.join(table_dir, 'phospho_' + setup + '.tsv')
    protein_outpath = os.path.join(table_dir, 'protein_' + setup + '.tsv')

    MQ.write_evidence_table(evidence, evidence_outpath)
    MQ.write_phosphosite_table(evidence, phosphosite_outpath, prob_cutoff=0.7)
    MQ.write_protein_table(evidence, protein_outpath)
