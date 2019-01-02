import os

import numpy as np

import gemreads_daytwo_adapt


_this_dir = os.path.split(os.path.abspath(__file__))[0]


# how to deal with alchemy version of gemsim?
def generate_se_non_circular_metagenome_reads(
        reference_dir,
        abundance_file,
        output_prefix,
        error_model,
        read_num=1000,
        read_len=75,
        quals=33):
    '''
    Generate single end reads from reference files based on abundance file.
    Arguments:
       reference_dir - directory containing reference genomes (fasta),
       abundance_file - file with genome abundances. Format: <ref_file_name_incl_ext>\tab<abundance_scale_0_to_1>
       output_prefix - output filename prefix,
       read_num - number of reads to generate. Default: 1000,
       read_len - length of reads to generate, default 75,
       error_model - read error model. *_single.gzip or *_paired.gzip. default: 'models/ill100v5_s.gzip' relative to install path
       quals - quality score offset. Default: 33

    Example:
    >>> import os
    >>> import magellan
    >>>
    >>> _this_dir = os.path.split(os.path.abspath(__file__))[0]
    >>> 
    >>> magellan.generate_se_non_circular_metagenome_reads(
    ...     os.path.join(_this_dir, '../test/_sample_data/community_genomes/'),
    ...     os.path.join(_this_dir, '../test/_sample_data/abund/abund.txt'),
    ...     os.path.join(_this_dir, '../test/_sample_data/gemsim_results'),
    ...     None)
    >>>
    >>>
    '''
    rel_ref_dir = os.path.relpath(reference_dir)

    if not error_model:
        error_model = os.path.join(
            _this_dir,
            'models',
            'ill100v5_s.gzip')
    if not os.path.exists(error_model):
        raise RuntimeError('model file does not exist at %s, please provide path to model file' % (error_model))

    gemreads_daytwo_adapt.run_simulation(
        None,
        rel_ref_dir,
        abundance_file,
        read_num,
        read_len,
        None,
        '',
        error_model,
        quals,
        output_prefix,
        None,
        None,
        False,
        True,
        False
        )
    
    if not os.path.exists(output_prefix + '_single.fastq'):
        raise RuntimeError('GemSIM run failed for ' + output_prefix)

    
def create_communities(
        simulation_names,
        community_pool,
        community_genomes_dir,
        abundance_files_base_dir,
        fastqs_base_dir,
        community_num_members=50,
        read_num=5000000,
        read_len=75):
    '''
    create mock communities from a specified member pool using log normal abundances
    Arguments:
       simulation_names: list of simulation_names, will be used in output file names
       community_pool: list of genome fasta file names under community_genomes_dir that contain the genomes of potential members in the mock community
       community_genomes_dir: base directory, under which all genomes of potential mock community members are located
       abundance_files_base_dir: directory under which abundance files for the mock community will be written, in the form of <genome_name>\t<abundance>
       fastqs_base_dir: directory under which the generated fastq files will be written
       community_num_members: number of members in each community (default: 50)
       read_num: number of reads simulated from each community (default: 5M)
       read_len: length of each simulated read (default: 75bp)

    Example:
    >>> import os
    >>> import magellan
    >>> import glob
    >>>
    >>> _this_dir = os.path.split(os.path.abspath(__file__))[0]
    >>>
    >>> magellan.create_communities(
    ... [0, 1, 2],
    ... ['GCF_000185525.fna', 'GCF_000185545.fna', 'GCF_001314995.fna'],
    ... os.path.join(_this_dir, '../test/_sample_data/community_genomes/'),
    ... os.path.join(_this_dir, '../test/_sample_data/abund/'),
    ... os.path.join(_this_dir, '../test/_sample_data/fastqs/'),
    ... community_num_members=2,
    ... read_num=1000)
    >>>
    >>>
    '''
    for sn in simulation_names:
        # create random abundances for community
        abund = np.random.lognormal(size=community_num_members)
        abund = abund / abund.sum()
        random_members = np.random.choice(community_pool, community_num_members, replace=False)
        abund_file = os.path.join(abundance_files_base_dir, 'abundance_' + str(sn) + '.txt')
        with open(abund_file, 'w') as f:
            for m, a in zip(random_members, abund):
                f.write(m + '\t' + str(a) + '\n')

        # create read files with adapted GemSIM
        out_prefix = os.path.join(fastqs_base_dir, 'fastq_' + str(sn))
        generate_se_non_circular_metagenome_reads(
            community_genomes_dir,
            abund_file,
            out_prefix,
            None,
            read_num=read_num,
            read_len=read_len)

