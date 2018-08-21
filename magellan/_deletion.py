import pandas as pd
import numpy as np


class Cnv:
    def __init__(
            self,
            df_genes_spec_path,
            df_abund_path,
            df_nr_path,
            df_length_path):
        self.df_genes_spec_path = df_genes_spec_path
        self.df_abund_path = df_abund_path
        self.df_nr_path = df_nr_path
        self.df_length_path = df_length_path
        
    def set_table_of_k(self):
        self.df_k = pd.DataFrame(
            self.df_abund.as_matrix() * np.atleast_2d(self.gene_lengths+0.075).T * self.arr_exp_reads_per_kb,
            columns=self.df_abund.columns,
            index=self.df_abund.index).astype('float32')
        
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
        
        df_abund_mat = self.df_abund.as_matrix()
        df_expected_means = self.df_k*self.df_nr*1
        df_expected_means = df_expected_means.astype('float32')
        df_expected_std = np.sqrt(self.df_k*self.df_nr*1*(1-self.df_k*1))
        df_expected_std = df_expected_std.astype('float32')
        
        z_1_mat = (self.df_genes_spec.as_matrix()
                   - df_expected_means.as_matrix()) / df_expected_std.as_matrix()
        z_of_0_for_1 = (0 - df_expected_means.as_matrix()) / df_expected_std.as_matrix()
        
        z_1_mat[z_of_0_for_1 > max_0_zscore] = np.nan
        z_1_mat[df_abund_mat < min_abund] = np.nan
        
        z_1_mat_deletions = z_1_mat.copy()
        z_1_mat_deletions[z_1_mat_deletions >= max_deletion_zscore] = np.nan
        z_1_mat[z_1_mat < max_deletion_zscore] = np.nan
        
        del_or_not_mat = z_1_mat.copy()
        del_or_not_mat[~np.isnan(z_1_mat_deletions)] = 0
        del_or_not_mat[~np.isnan(z_1_mat)] = 1
        
        self.df_del_or_not = pd.DataFrame(
            del_or_not_mat,
            columns=self.df_genes_spec.columns,
            index=self.df_genes_spec.index)
        
    def clean_deletions_and_low_abunds(
            self,
            min_abund,
            max_0_zscore,
            max_deletion_zscore):
        
        df_abund_mat = self.df_abund.as_matrix()
        df_expected_means = self.df_k*self.df_nr*1
        df_expected_means = df_expected_means.astype('float32')
        df_expected_std = np.sqrt(self.df_k*self.df_nr*1*(1-self.df_k*1))
        df_expected_std = df_expected_std.astype('float32')
        
        z_1_mat=(self.df_genes_spec.as_matrix()
                 - df_expected_means.as_matrix())/df_expected_std.as_matrix()
        z_of_0_for_1 = (0 - df_expected_means.as_matrix()) / df_expected_std.as_matrix()
        
        z_1_mat[z_of_0_for_1 > max_0_zscore] = np.nan
        z_1_mat[df_abund_mat < min_abund] = np.nan
        z_1_mat[z_1_mat < max_deletion_zscore] = np.nan
        
        self.df_genes_spec = self.remove_all_table_entries_by_correction_matrix(
            z_1_mat,
            self.df_genes_spec)
        
        self.df_abund = self.remove_all_table_entries_by_correction_matrix(
            z_1_mat,
            self.df_abund)
        
        self.df_nr = self.remove_all_table_entries_by_correction_matrix(
            z_1_mat,
            self.df_nr)
        
    def remove_genes_with_same_uniref_in_same_client(self):
        df_uniq = self.df_genes_spec.reset_index()\
            .set_index('specie')\
            .groupby('uniref90', as_index=False)\
            .count()
        
        df_map = self.df_genes_spec.reset_index()[['specie', 'uniref90']]
        
        df_uniq = pd.merge(df_uniq, df_map, on='uniref90')
        df_uniq = df_uniq.set_index(['specie', 'uniref90'])
        df_uniq = df_uniq.astype('float32')
        df_uniq_mat = df_uniq.as_matrix()
        df_uniq_mat[df_uniq_mat != 1] = np.nan
        
        self.df_genes_spec = self.remove_all_table_entries_by_correction_matrix(
            df_uniq_mat,
            self.df_genes_spec)
        
        self.df_abund = self.remove_all_table_entries_by_correction_matrix(
            df_uniq_mat,
            self.df_abund)
        
        self.df_nr = self.remove_all_table_entries_by_correction_matrix(
            df_uniq_mat,
            self.df_nr)
        
    def sort_genes_by_cnv(
            self,
            output_path,
            min_gene_length=700.0,
            max_num_reads_for_search_exp_reads=0.0002,
            num_bins_for_search_exp_reads=50,
            start_bin_search_exp_reads=10,
            min_abund=0.004,
            max_0_zscore=-2,
            max_deletion_zscore=-2,
            min_num_of_client_for_gene=67,
            max_possible_copy_number=60,
            max_frac_for_sort=0.95,
            min_frac_for_sort=0.5,
            interval_frac_for_sort=0.05,
            num_of_bins_for_sort=50,
            outliers_perc_for_sort=2.5):
        
        self.df_genes_spec = pd.read_pickle(self.df_genes_spec_path).astype('float32')
        self.df_abund = pd.read_pickle(self.df_abund_path).astype('float32')
        self.df_nr = pd.read_pickle(self.df_nr_path).astype('float32')
        self.df_length = pd.read_pickle(self.df_length_path)

        if not self.check_same_genes_and_clients():
            print 'the DataFrames do not have the same columns/index in the same order'
            return None
        
        self.drop_duplicate_uniref_in_the_same_bacteria()
        self.remove_short_genes_and_get_lengths_vector(min_gene_length)
        
        self.set_expected_num_reads_per_client(
            max_num_reads_for_search_exp_reads,
            num_bins_for_search_exp_reads,
            start_bin_search_exp_reads)
        
        self.set_table_of_k()
        
        self.df_genes_spec = pd.DataFrame(
            self.df_genes_spec.as_matrix() * np.atleast_2d(self.gene_lengths).T,
            columns=self.df_genes_spec.columns,
            index=self.df_genes_spec.index)\
            .astype('float32')
        
        self.clean_deletions_and_low_abunds(
            min_abund,
            max_0_zscore,
            max_deletion_zscore)
        
        self.remove_genes_with_same_uniref_in_same_client()
        self.set_table_of_k()
        self.remove_genes_with_low_clients_num(min_num_of_client_for_gene)
        
        self.set_cn_table(
            min_num_of_client_for_gene,
            max_possible_copy_number)
        
        self.sort_cn_table_by_second()
        self.df_cn.to_pickle(output_path)
    
    def sort_genes_by_deletions(
            self,
            output_path,
            min_gene_length=700.0,
            max_num_reads_for_search_exp_reads=0.0002,
            num_bins_for_search_exp_reads=50,
            start_bin_search_exp_reads=10,
            min_abund=0.004,
            max_0_zscore=-2,
            max_deletion_zscore=-2,
            min_num_of_client_for_gene=67):

        self.df_genes_spec = pd.read_pickle(self.df_genes_spec_path).astype('float32')
        self.df_abund = pd.read_pickle(self.df_abund_path).astype('float32')
        self.df_nr = pd.read_pickle(self.df_nr_path).astype('float32')
        self.df_length = pd.read_pickle(self.df_length_path)
        
        if not self.check_same_genes_and_clients():
            print 'the DataFrames do not have the same columns/index in the same order'
            return None
        
        self.drop_duplicate_uniref_in_the_same_bacteria()
        self.remove_short_genes_and_get_lengths_vector(min_gene_length)

        self.set_expected_num_reads_per_client(
            max_num_reads_for_search_exp_reads,
            num_bins_for_search_exp_reads,
            start_bin_search_exp_reads)
        
        self.set_nan_deletions_as_zero()
        self.set_table_of_k()
        
        self.df_genes_spec = pd.DataFrame(
            self.df_genes_spec.as_matrix()*np.atleast_2d(self.gene_lengths).T,
            columns=self.df_genes_spec.columns,
            index=self.df_genes_spec.index)\
            .astype('float32')
        
        self.set_deletions_table_and_clean_low_abunds(
            min_abund,
            max_0_zscore,
            max_deletion_zscore)
        
        self.df_del_or_not = self.df_del_or_not[
            self.df_del_or_not.count(1) > min_num_of_client_for_gene]
        
        self.sort_deletion_table_minority()
        self.df_del_or_not.to_pickle(output_path)
        
    def sort_deletion_table_by_enthropy(self):
        total_count = self.df_del_or_not.count(1)
        ones_count = self.df_del_or_not.sum(1)
        zero_count = total_count-ones_count

        zero_ratio = zero_count / total_count
        ones_ratio = ones_count / total_count

        self.df_del_or_not['zero_ratio'] = zero_ratio
        self.df_del_or_not['ones_ratio'] = ones_ratio
        
        self.df_del_or_not['ratio'] = self.df_del_or_not[['zero_ratio', 'ones_ratio']].min(1)
        self.df_del_or_not = self.df_del_or_not.sort_values('ratio', ascending=False)
    
    def set_nan_deletions_as_zero(self):

        gs_mat = self.df_genes_spec.as_matrix()
        abund_mat = self.df_abund.as_matrix()

        gs_mat[np.isnan(gs_mat) & ~np.isnan(abund_mat)] = 0.00001
        
        self.df_genes_spec = pd.DataFrame(
            gs_mat,
            columns=self.df_genes_spec.columns,
            index=self.df_genes_spec.index)
    
    def remove_genes_with_low_clients_num(self, min_num_of_client_for_gene):
        
        self.df_genes_spec = self.df_genes_spec[
            self.df_genes_spec.count(1) > min_num_of_client_for_gene]
        
        self.df_abund = self.df_abund[
            self.df_abund.count(1) > min_num_of_client_for_gene]
        
        self.df_nr = self.df_nr[self.df_nr.count(1) > min_num_of_client_for_gene]
        self.df_k = self.df_k[self.df_k.count(1) > min_num_of_client_for_gene]
        
    def set_cn_table(self, min_num_of_client_for_gene, max_possible_copy_number):
        cn_mat = self.df_genes_spec.as_matrix().copy()
        cn_mat[~np.isnan(self.df_genes_spec.as_matrix())] = 1

        df_expected_means = self.df_k * self.df_nr * 1
        df_expected_means = df_expected_means.astype('float32')
        df_expected_std = np.sqrt(self.df_k * self.df_nr * 1 * (1 - self.df_k * 1))
        df_expected_std = df_expected_std.astype('float32')
        best_z_mat = abs(
            (self.df_genes_spec.as_matrix() - df_expected_means.as_matrix())
            / df_expected_std.as_matrix())

        for i in range(1, max_possible_copy_number):
            
            df_expected_means = self.df_k*self.df_nr*i
            df_expected_means = df_expected_means.astype('float32')
            df_expected_std = np.sqrt(self.df_k*self.df_nr*i*(1-self.df_k*i))
            df_expected_std = df_expected_std.astype('float32')

            curr_z_mat = abs(
                (self.df_genes_spec.as_matrix() - df_expected_means.as_matrix())
                / df_expected_std.as_matrix())
            
            cn_mat[best_z_mat > curr_z_mat] = i
            best_z_mat[best_z_mat > curr_z_mat] = curr_z_mat[best_z_mat > curr_z_mat]
            
        self.df_cn = pd.DataFrame(
            cn_mat,
            columns=self.df_genes_spec.columns,
            index=self.df_genes_spec.index)

    def set_expected_num_reads_per_client(
            self,
            max_num_reads_for_search,
            num_bins_for_search_exp_reads,
            start_bin_search_exp_reads):
        
        lst_exp_reads = []
        df_exp_reads = (self.df_genes_spec / (self.df_abund * self.df_nr)).astype('float32')
        
        for i in df_exp_reads.columns:
            x = np.histogram(
                df_exp_reads[[i]].dropna()[df_exp_reads[[i]].dropna()[i] < max_num_reads_for_search][i],
                num_bins_for_search_exp_reads)
            lst_exp_reads.append(
                x[1][x[0][start_bin_search_exp_reads:x[0].shape[0]].argmax() + start_bin_search_exp_reads])
            
        self.arr_exp_reads_per_kb = np.array(lst_exp_reads)

    def check_same_genes_and_clients(self):
        res = True
        res = res and all(self.df_genes_spec.columns == self.df_abund.columns)
        res = res and all(self.df_abund.columns == self.df_nr.columns)
        res = res and all(self.df_genes_spec.index == self.df_abund.index)
        res = res and all(self.df_abund.index == self.df_nr.index)

        return res

    def drop_duplicate_uniref_in_the_same_bacteria(self):
        self.df_genes_spec = self.df_genes_spec\
            .reset_index()\
            .drop_duplicates(['specie', 'uniref90'], keep=False)\
            .set_index(['specie', 'uniref90'])
        
        self.df_abund = self.df_abund\
            .reset_index()\
            .drop_duplicates(['specie', 'uniref90'], keep=False)\
            .set_index(['specie', 'uniref90'])
        
        self.df_nr = self.df_nr\
            .reset_index()\
            .drop_duplicates(['specie', 'uniref90'], keep=False)\
            .set_index(['specie', 'uniref90'])
        
        self.df_length = self.df_length\
            .drop_duplicates(['specie', 'uniref90'], keep=False)
        
    def remove_short_genes_and_get_lengths_vector(self, min_length):
        self.df_genes_spec, self.gene_lengths = self.remove_short_genes_for_df(
            min_length,
            self.df_genes_spec)
        
        self.df_abund, self.gene_lengths = self.remove_short_genes_for_df(
            min_length,
            self.df_abund)
        
        self.df_nr, self.gene_lengths = self.remove_short_genes_for_df(
            min_length,
            self.df_nr)
        
    def remove_short_genes_for_df(self, min_length, df):
        
        df = df.reset_index()
        df = pd.merge(df, self.df_length, on=['specie', 'uniref90'])
        df = df[df.length > min_length]
        
        gene_lengths = np.array(df.length)
        df = df.drop(['length'], 1)
        df = df.set_index(['specie', 'uniref90'])
        gene_lengths = gene_lengths/1000.0
        
        return df.astype('float32'), gene_lengths
    
    def smallestSubWithSum(self, arr, n, x):
        min_len = n + 1

        for start in range(0, n):
            curr_sum = arr[start]
            if (curr_sum > x):
                return 1
            for end in range(start + 1, n):
                curr_sum += arr[end]
                if curr_sum > x and (end - start + 1) < min_len:
                    min_len = (end - start + 1)

        return min_len
    
    def sort_cn_table_by_smalest_window(
            self,
            max_frac,
            min_frac,
            interval_frac,
            num_of_bins,
            outliers_perc):
        
        df_cn_mat = self.df_cn.as_matrix().copy()
        M = np.nanmean(df_cn_mat, 1)
        df_cn_mat = df_cn_mat / np.atleast_2d(M).T
        
        lst_sorting_vals = []
        for i in np.arange(max_frac, min_frac - interval_frac, interval_frac*(-1)):

            lst_res = []
            for row in range(df_cn_mat.shape[0]):
                r = df_cn_mat[row]
                r= r[~np.isnan(r)]
                q025 = np.nanpercentile(r, outliers_perc)
                q0975 = np.nanpercentile(r, 100 - outliers_perc)
                r = r[r >= q025]
                r = r[r <= q0975]
                arr = np.histogram(r, num_of_bins)[0]
                lst_res.append(self.smallestSubWithSum(arr, len(arr), len(r) * i))
                
            self.df_cn['data'+str(i)] = lst_res
            lst_sorting_vals.append('data'+str(i))
        self.df_cn = self.df_cn.sort_values(lst_sorting_vals, ascending=False)
    
        def sort_deletion_table_minority(self):

            total_count = self.df_del_or_not.count(1)
            ones_count = self.df_del_or_not.sum(1)
            zero_count = total_count - ones_count
            
            self.df_del_or_not['zero_count'] = zero_count
            self.df_del_or_not['ones_count'] = ones_count
            
            self.df_del_or_not['minority'] = self.df_del_or_not[['zero_count', 'ones_count']].min(1)
            self.df_del_or_not = self.df_del_or_not.sort_values('minority', ascending=False)

        def sort_cn_table_by_second(self):
            self.df_cn['sec'] = [self.get_second(list(self.df_cn.iloc[i]))
                                 for i in range(self.df_cn.shape[0])]
            
            self.df_cn = self.df_cn.sort_values('sec', ascending=False)

        def get_second(self, lst):
            
            d = {x: lst.count(x) for x in lst}
            res_lst = sorted(d, key=d.get, reverse=True)
            
            if len(res_lst) == 1:
                return 0
            else:
                return d[res_lst[1]]
