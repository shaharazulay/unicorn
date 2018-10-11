import unittest
import os
import pandas as pd

import magellan

curr_path = os.path.abspath(os.path.dirname(__file__))
_sample_data_dir = os.path.join(curr_path, '_sample_data')


class _TestDropDuplication(unittest.TestCase):
    def test_basic(self):
        
        deldet = magellan.DeletionDetector(
            os.path.join(_sample_data_dir, 'client_gene_reads_with_dup.pkl'),
            os.path.join(_sample_data_dir, 'sample_abundances.pkl'),
            os.path.join(_sample_data_dir, 'sample_total_reads.pkl'),
            os.path.join(_sample_data_dir, 'gene_lengths_with_dup.pkl'),
        )
        origin_df = pd.read_pickle(
            os.path.join(_sample_data_dir, 'client_gene_reads_with_dup.pkl'))
        
        self.assertEqual(
            deldet.client_gene_reads.shape[0],
            origin_df.shape[0])
        
        deldet.drop_duplicate_uniref_in_the_same_bacteria()

        self.assertLess(
            deldet.client_gene_reads.shape[0],
            origin_df.shape[0])
            
        current_read_index = deldet.client_gene_reads.index.values.tolist()
        self.assertEqual(len(current_read_index), len(set(current_read_index)))

        current_gene_length_index = deldet.gene_length.index.values.tolist()
        self.assertEqual(len(current_gene_length_index), len(set(current_gene_length_index)))


class _TestShortGeneRemoval(unittest.TestCase):
    def test_basic(self):
        deldet = magellan.DeletionDetector(
            os.path.join(_sample_data_dir, 'client_gene_reads_with_dup.pkl'),
            os.path.join(_sample_data_dir, 'sample_abundances.pkl'),
            os.path.join(_sample_data_dir, 'sample_total_reads.pkl'),
            os.path.join(_sample_data_dir, 'gene_lengths_with_dup.pkl'),
        )
                
        deldet.drop_duplicate_uniref_in_the_same_bacteria()

        dedup_df = deldet.client_gene_reads.copy()

        min_len = 20
        deldet.remove_short_genes_and_get_lengths_vector(min_len)

        self.assertLess(
            deldet.client_gene_reads.shape[0],
            dedup_df.shape[0])

        relevant_genes = deldet.relevant_gene_lengths

        self.assertEqual(len(relevant_genes), len(deldet.client_gene_reads))
        self.assertLess(len(relevant_genes), len(deldet.gene_length))
    
        self.assertGreaterEqual(relevant_genes.length.max(), min_len / 1000.)

        
if __name__ == '__main__':
    unittest.main()
