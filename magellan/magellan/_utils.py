from _deletion import Cnv
import argparse


def _add_client_gene_reads_path_to_parser(parser):
    parser.add_argument(
        '--client_gene_reads_path',
        dest='client_gene_reads_path',
        required='True',
        help='the path of the reads dataframe')

    
def _add_abundance_path_to_parser(parser):
    parser.add_argument(
        '--df_abund_path',
        dest='df_abund_path',
        required='True',
        help='the path of the species abandance dataframe')

    
def _add_reads_path_to_parser(parser):
    parser.add_argument(
        '--df_nr_path',
        dest='df_nr_path',
        required='True',
        help='the path of the total_reads dataframe')

    
def _add_gene_length_path_to_parser(parser):
    parser.add_argument(
        '--df_length_path',
        dest='df_length_path',
        required='True',
        help='the path of the gene_length dataframe')


def _add_output_path_to_parser(parser):
    parser.add_argument(
        '--output_path',
        dest='output_path',
        required='True',
        help='the path of the output dataframe(pickle)')
    

def _add_deletion_flag_to_parser(parser):
    parser.add_argument(
        '--deletions',
        dest='deletions',
        action='store_true',
        default=False,
        help='if stated, the genes will be sorted by deletions')


def _add_min_gene_length_to_parser(parser):
    parser.add_argument(
        '--min_gene_length',
        dest='min_gene_length',
        type=int,
        default=700,
        help='genes shorter than this value will be ignored')


def _add_max_num_reads_to_parser(parser):
    parser.add_argument(
        '--max_num_reads_for_search_exp_reads',
        dest='max_num_reads_for_search_exp_reads',
        type=float,
        default=0.0002,
        help='max num of reads to look at when searching expected reads per client')


def _add_num_bins_to_parser(parser):
    parser.add_argument(
        '--num_bins_for_search_exp_reads',
        dest='num_bins_for_search_exp_reads',
        type=int,
        default=50,
        help='num of bins to use when searching for expected reads per client')

    
def _add_start_bin_to_parser(parser):
    parser.add_argument(
        '--start_bin_search_exp_reads',
        dest='start_bin_search_exp_reads',
        type=int,
        default=10,
        help='bin to start from when searching for expected reads per client')

    
def _add_abundance_lower_bound_to_parser(parser):
    parser.add_argument(
        '--min_abund',
        dest='min_abund',
        type=float,
        default=0.004,
        help='species with lower abundance will be ignored')


def _add_max_zero_z_score_to_parser(parser):
    parser.add_argument(
        '--max_0_zscore',
        dest='max_0_zscore',
        type=int,
        default=-2,
        help='if 0 gets a z score higher than this value, the gene will be ignored')

    
def _add_max_deletion_z_score_bound_to_parser(parser):
    parser.add_argument(
        '--max_deletion_zscore',
        dest='max_deletion_zscore',
        type=int,
        default=-2,
        help='z scores below this value will be considered as deletions')

    
def _add_min_clients_per_gene_parser(parser):
    parser.add_argument(
        '--min_num_of_client_for_gene',
        dest='min_num_of_client_for_gene',
        type=int,
        default=67,
        help='gene that appear in less then this number of clients will be ignored')

    
def _add_max_expeted_copy_number_to_parser(parser):
    parser.add_argument(
        '--max_possible_copy_number',
        dest='max_possible_copy_number',
        type=int,
        default=60,
        help='max copy number considered (above this value it will be considered this value)')

    
def _add_max_sort_frac_to_parser(parser):
    parser.add_argument(
        '--max_frac_for_sort',
        dest='max_frac_for_sort',
        type=float,
        default=0.95,
        help='the highst fraction of data for smallest window sort')

    
def _add_min_sort_frac_to_parser(parser):
    parser.add_argument(
        '--min_frac_for_sort',
        dest='min_frac_for_sort',
        type=float,
        default=0.5,
        help='the lowest fraction of data for smallest window sort')

    
def _add_interval_sort_frac_to_parser(parser):
    parser.add_argument(
        '--interval_frac_for_sort',
        dest='interval_frac_for_sort',
        type=float,
        default=0.05,
        help='the interval between fractions of data for smallest window sort')

    
def _add_outlier_sort_perc_to_parser(parser):
    parser.add_argument(
        '--outliers_perc_for_sort',
        dest='outliers_perc_for_sort',
        type=float,
        default=2.5,
        help='outlier percent to remove from each size for window sort')

    
def _add_num_sort_bins_to_parser(parser):
    parser.add_argument(
        '--num_of_bins_for_sort',
        dest='num_of_bins_for_sort',
        type=int,
        default=50,
        help='number of bins for smallest window sort')


def _main():
    parser = argparse.ArgumentParser(description='main_parser')
    
    _add_abundance_lower_bound_to_parser(parser)
    _add_abundance_path_to_parser(parser)
    _add_deletion_flag_to_parser(parser)
    _add_gene_length_path_to_parser(parser)
    _add_client_gene_reads_path_to_parser(parser)
    _add_interval_sort_frac_to_parser(parser)
    _add_max_deletion_z_score_bound_to_parser(parser)
    _add_max_expeted_copy_number_to_parser(parser)
    _add_max_num_reads_to_parser(parser)
    _add_max_zero_z_score_to_parser(parser)
    _add_min_clients_per_gene_parser(parser)
    _add_min_gene_length_to_parser(parser)
    _add_min_sort_frac_to_parser(parser)
    _add_num_bins_to_parser(parser)
    _add_num_sort_bins_to_parser(parser)
    _add_outlier_sort_perc_to_parser(parser)
    _add_output_path_to_parser(parser)
    _add_reads_path_to_parser(parser)
    _add_start_bin_to_parser(parser)

    args = parser.parse_args()

    cnv = Cnv(
        args.client_gene_reads_path,
        args.df_abund_path,
        args.df_nr_path,
        args.df_length_path)

    if args.deletions:
        cnv.sort_genes_by_deletions(
            args.output_path,
            args.min_gene_length,
            args.max_num_reads_for_search_exp_reads,
            args.num_bins_for_search_exp_reads,
            args.start_bin_search_exp_reads,
            args.min_abund,
            args.max_0_zscore,
            args.max_deletion_zscore,
            args.min_num_of_client_for_gene)
        
    else:
        cnv.sort_genes_by_cnv(
            args.output_path,
            args.min_gene_length,
            args.max_num_reads_for_search_exp_reads,
            args.num_bins_for_search_exp_reads,
            args.start_bin_search_exp_reads,
            args.min_abund,
            args.max_0_zscore,
            args.max_deletion_zscore,
            args.min_num_of_client_for_gene,
            args.max_possible_copy_number,
            args.max_frac_for_sort,
            args.min_frac_for_sort,
            args.interval_frac_for_sort,
            args.num_of_bins_for_sort,
            args.outliers_perc_for_sort)

