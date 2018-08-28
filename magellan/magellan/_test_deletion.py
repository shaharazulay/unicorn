import unittest
import os
import pandas as pd

_sample_data_dir = './sample_data/'


class _TestDropDuplication(unittest.TestCase):
    def test_basic(self):
        import _deletion
        deldet = _deletion.DeletionDetector(
            os.path.join(_sample_data_dir, 'client_gene_reads_with_dup.pkl'),
            os.path.join(_sample_data_dir, 'sample_abundances.pkl'),
            os.path.join(_sample_data_dir, 'sample_total_reads.pkl'),
            os.path.join(_sample_data_dir, 'gene_lengths_with_dup.pkl'),
        )
        origin_df = pd.read_pickle(os.path.join(_sample_data_dir, 'client_gene_reads_with_dup.pkl'))
        
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


if __name__ == '__main__':
    unittest.main()
