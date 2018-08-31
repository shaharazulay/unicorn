import pandas as pd
import numpy as np


class DeletionDetector:
    def __init__(
            self,
            client_gene_reads_path,
            species_abund_path,
            client_total_reads_path,
            gene_length_path):
        self.client_gene_reads_path = client_gene_reads_path
        self.species_abund_path = species_abund_path
        self.client_total_reads_path = client_total_reads_path
        self.gene_length_path = gene_length_path

        self.client_gene_reads = pd.read_pickle(self.client_gene_reads_path).astype('float32')
        self.species_abund = pd.read_pickle(self.species_abund_path).astype('float32')
        self.client_total_reads = pd.read_pickle(self.client_total_reads_path).astype('float32')
        # Tmp tali - why is this the only one with different structure (species/uniref not in index)
        self.gene_length = pd.read_pickle(self.gene_length_path)
        
    def set_table_of_k(self):
        self.df_k = pd.DataFrame(
            self.species_abund.as_matrix() * np.atleast_2d(self.relevant_gene_lengths+0.075).T * self.arr_exp_reads_per_kb,
            columns=self.species_abund.columns,
            index=self.species_abund.index).astype('float32')
        
    def remove_all_table_entries_by_correction_matrix(self, corr_mat, df):
        df_mat = df.as_matrix()
        df_mat[np.isnan(corr_mat)] = np.nan
        df = pd.DataFrame(df_mat, columns=df.columns, index=df.index)
        return df
        
    def set_deletions_table_and_clean_low_abunds(
            self,
            min_abund,
            max_0_zscore,
            max_deletion_zscore):
        
        species_abund_mat = self.species_abund.as_matrix()
        df_expected_means = self.df_k*self.client_total_reads*1
        df_expected_means = df_expected_means.astype('float32')
        df_expected_std = np.sqrt(self.df_k*self.client_total_reads*1*(1-self.df_k*1))
        df_expected_std = df_expected_std.astype('float32')
        
        z_1_mat = (self.client_gene_reads.as_matrix()
                   - df_expected_means.as_matrix()) / df_expected_std.as_matrix()
        z_of_0_for_1 = (0 - df_expected_means.as_matrix()) / df_expected_std.as_matrix()
        
        z_1_mat[z_of_0_for_1 > max_0_zscore] = np.nan
        z_1_mat[species_abund_mat < min_abund] = np.nan
        
        z_1_mat_deletions = z_1_mat.copy()
        z_1_mat_deletions[z_1_mat_deletions >= max_deletion_zscore] = np.nan
        z_1_mat[z_1_mat < max_deletion_zscore] = np.nan
        
        del_or_not_mat = z_1_mat.copy()
        del_or_not_mat[~np.isnan(z_1_mat_deletions)] = 0
        del_or_not_mat[~np.isnan(z_1_mat)] = 1
        
        self.df_del_or_not = pd.DataFrame(
            del_or_not_mat,
            columns=self.client_gene_reads.columns,
            index=self.client_gene_reads.index)
            
    def sort_genes_by_deletions(
            self,
            output_path,
            min_gene_length=700.0,  # Tmp tali - needed? will be nan anyway for genes in which we aren't certain and doesn't need to be a param to set expected
            max_num_reads_for_search_exp_reads=0.0002,
            num_bins_for_search_exp_reads=50,
            start_bin_search_exp_reads=10,
            min_abund=0.004,
            max_0_zscore=-2,
            max_deletion_zscore=-2,
            min_num_of_client_for_gene=67):  # Tmp tali - set by pct or just output number of genes and let user filter?

        # Tmp tali - new flow to implement:
        # 1. specify bacteria (mandatory?)
        # 2. compute expected reads based on all bacteria
        # 3. compute deletions only on specified bacteria
        # 4. sort?
        
        #Tmp tali - copied this also to init. Currently, we run members over, so need to keep it still here. best if we don't
        self.client_gene_reads = pd.read_pickle(self.client_gene_reads_path).astype('float32')
        self.species_abund = pd.read_pickle(self.species_abund_path).astype('float32')
        self.client_total_reads = pd.read_pickle(self.client_total_reads_path).astype('float32')
        # Tmp tali - why is this the only one with different structure (species/uniref not in index)
        self.gene_length = pd.read_pickle(self.gene_length_path)

        self.drop_duplicate_uniref_in_the_same_bacteria()
        self.remove_short_genes_and_get_lengths_vector(min_gene_length)

        self.inflate_matrices()
        
        self.set_expected_num_reads_per_client(
            max_num_reads_for_search_exp_reads,
            num_bins_for_search_exp_reads,
            start_bin_search_exp_reads)
        
        self.set_nan_deletions_as_zero()
        self.set_table_of_k()

        # Tmp Tali - many cases in which members are modified. Make sure that if function is run with different bacteria/param chnages to members still hold
        self.client_gene_reads = pd.DataFrame(
            self.client_gene_reads.as_matrix()*np.atleast_2d(self.relevant_gene_lengths).T,
            columns=self.client_gene_reads.columns,
            index=self.client_gene_reads.index)\
            .astype('float32')
        
        self.set_deletions_table_and_clean_low_abunds(
            min_abund,
            max_0_zscore,
            max_deletion_zscore)
        
        self.df_del_or_not = self.df_del_or_not[
            self.df_del_or_not.count(1) > min_num_of_client_for_gene]
        
        self.sort_deletion_table_minority()
        self.df_del_or_not.to_pickle(output_path)
    
    def set_nan_deletions_as_zero(self):
        self.client_gene_reads[(self.client_gene_reads.isnull() & ~self.species_abund.isnull())] = 0.00001
    
    def set_expected_num_reads_per_client(
            self,
            max_num_reads_for_search,
            num_bins_for_search_exp_reads,
            start_bin_search_exp_reads):
        
        lst_exp_reads = []
        df_exp_reads = (self.client_gene_reads / (self.species_abund * self.client_total_reads)).astype('float32')

        # Tmp Tali - apply function to all columns instead of iterating? also, perhaps this can be simplified
        for i in df_exp_reads.columns:
            x = np.histogram(
                df_exp_reads[[i]].dropna()[df_exp_reads[[i]].dropna()[i] < max_num_reads_for_search][i],
                num_bins_for_search_exp_reads)
            lst_exp_reads.append(
                x[1][x[0][start_bin_search_exp_reads:x[0].shape[0]].argmax() + start_bin_search_exp_reads])
            
        self.arr_exp_reads_per_kb = np.array(lst_exp_reads)

    def check_same_genes_and_clients(self):
        res = True
        res = res and all(self.client_gene_reads.columns == self.species_abund.columns)
        res = res and all(self.species_abund.columns == self.client_total_reads.columns)
        res = res and all(self.client_gene_reads.index == self.species_abund.index)
        res = res and all(self.species_abund.index == self.client_total_reads.index)

        return res

    def drop_duplicate_uniref_in_the_same_bacteria(self):
        self.client_gene_reads = self.client_gene_reads\
            .reset_index()\
            .drop_duplicates(['species', 'uniref90'], keep=False)\
            .set_index(['species', 'uniref90'])
                
        self.gene_length = self.gene_length.reset_index()\
            .drop_duplicates(['species', 'uniref90'], keep=False)
        
    def remove_short_genes_and_get_lengths_vector(self, min_length):
        # Tmp tali - does this need to be a member?
        self.relevant_gene_lengths = self.gene_length[self.gene_length['length'] >= min_length]
        self.relevant_gene_lengths['length'] = self.relevant_gene_lengths['length']/1000.
        
        self.client_gene_reads = self.client_gene_reads.reset_index()
        self.client_gene_reads = pd.merge(
            self.client_gene_reads,
            self.relevant_gene_lengths[['species', 'uniref90']],
            on=['species', 'uniref90'],
            how='inner')
        
        self.client_gene_reads = self.client_gene_reads.set_index(['species', 'uniref90'])
        self.client_gene_reads = self.client_gene_reads.astype('float32')

        # Tmp tali - this is turned into an array. Safer if kept as DF
        self.relevant_gene_lengths.set_index(['species', 'uniref90'], inplace=True)
        self.relevant_gene_lengths = self.relevant_gene_lengths.loc[self.client_gene_reads.index.values.tolist()]
    
    def inflate(self):
        self.species_abund = self.client_gene_reads.reset_index()[['species', 'uniref90']].merge(
            self.species_abund.reset_index(), on=['species'], how='left')
        self.species_abund.set_index(['species', 'uniref90'], inplace=True)

        self.client_total_reads = self.client_gene_reads.reset_index()[['species', 'uniref90']].merge(
            self.client_total_reads.reset_index(), on=['species'], how='left')
        self.client_total_reads.set_index(['species', 'uniref90'], inplace=True)

        if not self.check_same_genes_and_clients():
            raise RuntimeError('the input DataFrames do not have the same columns/index in the same order')

    def sort_deletion_table_minority(self):
        total_count = self.df_del_or_not.count(1)
        ones_count = self.df_del_or_not.sum(1)
        zero_count = total_count-ones_count
        self.df_del_or_not['zero_count'] = zero_count
        self.df_del_or_not['ones_count'] = ones_count
        self.df_del_or_not['minority'] = self.df_del_or_not[['zero_count', 'ones_count']].min(1)
        self.df_del_or_not = self.df_del_or_not.sort_values('minority', ascending=False)
    
