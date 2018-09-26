import unittest
import os
import pandas as pd
import glob

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

        
class _TestBowtie2Index(unittest.TestCase):
    def test_basic(self):
        magellan.bowtie2_build_index(
            os.path.join(_sample_data_dir, 'sample_bac_genome_short.fa'),
            os.path.join(_sample_data_dir, 'test_index'))

        index_files = glob.glob(os.path.join(_sample_data_dir, 'test_index.*'))
        for i in index_files:
            os.remove(i)
            
        magellan.bowtie2_build_index(
            os.path.join(_sample_data_dir, 'sample_bac_genome_short.fa'),
            os.path.join(_sample_data_dir, 'test_index_offrate'),
            offrate=2
        )

        offrate_index_files = glob.glob(os.path.join(_sample_data_dir, 'test_index_offrate.*'))
        for i in offrate_index_files:
            os.remove(i)

    def test_missing_fasta(self):
        with self.assertRaises(RuntimeError):
            magellan.bowtie2_build_index(
                os.path.join(_sample_data_dir, 'madeup.fa'),
                os.path.join(_sample_data_dir, 'test_index'))


class _TestBowtie2Mapping(unittest.TestCase):
    def test_basic(self):
        magellan.bowtie2_build_index(
            os.path.join(_sample_data_dir, 'sample_bac_genome_short.fa'),
            os.path.join(_sample_data_dir, 'test_index'))

        magellan.bowtie2_map(
            os.path.join(_sample_data_dir, 'sample_reads_for_mapper.fastq'),
            os.path.join(_sample_data_dir, 'test_index'),
            os.path.join(_sample_data_dir, 'sample_reads_for_mapper.sam'),
            os.path.join(_sample_data_dir, 'sample_reads_for_mapper.log'))
        
        index_files = glob.glob(os.path.join(_sample_data_dir, 'test_index'))
        
        for i in index_files:
            os.remove(i)

        os.remove(os.path.join(_sample_data_dir, 'sample_reads_for_mapper.sam'))
        os.remove(os.path.join(_sample_data_dir, 'sample_reads_for_mapper.log'))
        
    def test_missing_inputs(self):
        magellan.bowtie2_build_index(
            os.path.join(_sample_data_dir, 'sample_bac_genome_short.fa'),
            os.path.join(_sample_data_dir, 'test_index'))

        with self.assertRaises(RuntimeError):
            magellan.bowtie2_map(
                os.path.join(_sample_data_dir, 'madeup.fastq'),
                os.path.join(_sample_data_dir, 'test_index'),
                os.path.join(_sample_data_dir, 'sample_reads_for_mapper.sam'),
                os.path.join(_sample_data_dir, 'sample_reads_for_mapper.log'))

        with self.assertRaises(RuntimeError):
            magellan.bowtie2_map(
                os.path.join(_sample_data_dir, 'sample_reads_for_mapper.fastq'),
                os.path.join(_sample_data_dir, 'madeup_index'),
                os.path.join(_sample_data_dir, 'sample_reads_for_mapper.sam'),
                os.path.join(_sample_data_dir, 'sample_reads_for_mapper.log'))
            
        index_files = glob.glob(os.path.join(_sample_data_dir, 'test_index.*'))
        for i in index_files:
            os.remove(i)

        s = os.path.join(_sample_data_dir, 'sample_reads_for_mapper.sam')
        if os.path.isfile(s):
            os.remove(s)
        l = os.path.join(_sample_data_dir, 'sample_reads_for_mapper.sam')
        if os.path.isfile(l):
            os.remove(l)

            
class _TestAlignmentParser(unittest.TestCase):
    @classmethod
    def setUpClass(cls):

        magellan.bowtie2_build_index(
            os.path.join(_sample_data_dir, 'sample_bac_genome_short.fa'),
            os.path.join(_sample_data_dir, 'test_index'))

        magellan.bowtie2_map(
            os.path.join(_sample_data_dir, 'sample_reads_for_mapper.fastq'),
            os.path.join(_sample_data_dir, 'test_index'),
            os.path.join(_sample_data_dir, 'sample_reads_for_mapper.sam'),
            os.path.join(_sample_data_dir, 'sample_reads_for_mapper.log'))

    @classmethod
    def tearDownClass(cls):

        index_files = glob.glob(os.path.join(_sample_data_dir, 'test_index'))
        
        for i in index_files:
            os.remove(i)

        os.remove(os.path.join(_sample_data_dir, 'sample_reads_for_mapper.sam'))
        os.remove(os.path.join(_sample_data_dir, 'sample_reads_for_mapper_si.bam'))
        os.remove(os.path.join(_sample_data_dir, 'sample_reads_for_mapper_si.bam.bai'))
        os.remove(os.path.join(_sample_data_dir, 'sample_reads_for_mapper.log'))

    def test_basic(self):
        ap = magellan.AlignmentParser(os.path.join(_sample_data_dir, 'sample_reads_for_mapper.sam'))

        self.assertTrue(ap._alignment_file.check_index)

        ap.write_reads(os.path.join(_sample_data_dir, 'sample_reads_for_mapper_copy.sam'), write_sam=True)
        os.remove(os.path.join(_sample_data_dir, 'sample_reads_for_mapper_copy.sam'))

    def test_count(self):
        ap = magellan.AlignmentParser(os.path.join(_sample_data_dir, 'sample_reads_for_mapper.sam'))
        read_df = ap.count_all_references()
        read_df.set_index('name', inplace=True)
        read_dict = read_df['count'].to_dict()
        sam_df = pd.read_csv(os.path.join(_sample_data_dir, 'sample_reads_for_mapper.sam'), sep='\t', skiprows=5, header=None)
        sam_dict = sam_df.groupby(2)[0].count().to_dict()
        self.assertDictEqual(read_dict, sam_dict)

    def test_filter(self):
        ap = magellan.AlignmentParser(os.path.join(_sample_data_dir, 'sample_reads_for_mapper.sam'))
        
        it = ap.filter(read_filter_callback=lambda read: magellan.is_high_map_quality(read, 4))
        filtered_reads = []
        for r in it:
            filtered_reads.append(r.query_name)

        sam_df = pd.read_csv(os.path.join(_sample_data_dir, 'sample_reads_for_mapper.sam'), sep='\t', skiprows=5, header=None)
        filtered_sam_reads = sam_df[sam_df[4] > 4][0].values.tolist()
        self.assertListEqual(sorted(filtered_reads), sorted(filtered_sam_reads))
        

if __name__ == '__main__':
    unittest.main()
