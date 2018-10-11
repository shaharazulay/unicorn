
import subprocess
import pysam  # tested with version 0.13
import os
import pandas as pd
import glob


def _check_bowtie_install():
    '''
    Check if Bowtie2 is installed. Assumes linux platform
    '''
    try:
        bowtie2 = subprocess.Popen(['which', 'bowtie2'], stdout=subprocess.PIPE).communicate()[0].decode()
        bowtie2_version = subprocess.Popen(['bowtie2', '--version'], stdout=subprocess.PIPE).communicate()[0]
        bowtie2_version = bowtie2_version.decode().split()[2]
        print('Bowtie2 version ' + str(bowtie2_version) + ';  path: ' + str(bowtie2).strip())
    except Exception as err:
        print 'Bowtie2 installation error: ', str(err)
        raise
    

def bowtie2_build_index(fasta_file, index, offrate=None):
    '''
    Create bowtie2 index
    Arguments:
    fasta_file - containing all sequences to be indexed in fasta format
    index - base path for index
    offrate - bowtie2's offrate. Controls lookup speed at the expense of memory. Lower values lead to faster lookups.
              Default - None. Use bowtie2's default for this parameter
    '''
    cmd = 'bowtie2-build %s %s' % (fasta_file, index)
    if offrate:
        cmd = cmd + ' -o ' + str(offrate)
    print 'running: ', cmd
    _check_bowtie_install()
    
    _, err = subprocess.Popen(cmd,
                              stderr=subprocess.PIPE,
                              shell=True).communicate()

    if len(glob.glob(index + '.*')) == 0:
        raise RuntimeError('bowtie2 index files %s for %s are missing' % (index, fasta_file))
    
    
def bowtie2_map(
        fastq_file,
        index,
        out_sam,
        out_log,
        **kwargs):
    '''
    SE mapping using bowtie2.
    Arguments:
       fastq_file - file containing SE reads
       index - bowtie index to map to
       out_sam - path to sam file output
       out_log - path to log file
       kwargs - additional arguments provided to the bowtie mapper according to bowtie options. Optional. example {'L': 32}
    '''
    short_args = ['q', 'f', 'r', 'c', 's', 'u', '3', '5', 'N', 'L', 'i', 'k', 'a', 'D', 'R', 'I', 'X', 't', 'p']
    
    cmd = 'bowtie2 -U %s -x %s -S %s' % (fastq_file, index, out_sam)
    for k, v in kwargs.items():
        if k in short_args:
            cmd += ' -' + k + ' ' + str(v)
        else:
            cmd += ' --' + k + ' ' + str(v)
    cmd += ' 2> %s' % out_log
    print cmd
    _check_bowtie_install()

    _, err = subprocess.Popen(cmd,
                              stderr=subprocess.PIPE,
                              shell=True).communicate()
    if err:
        raise RuntimeError('bowtie mapping failed for fastq %s' % (fastq_file))

    if not os.path.isfile(out_log) or not os.path.isfile(out_sam):
        raise RuntimeError('output files missing for bowtie2 run on %s' % (fastq_file))
 

class AlignmentParser():
    '''
    A class for parsing, filtering and aggregating reads from sam and bam files.

    Example:
    >>> import os
    >>> import pandas as pd
    >>> import magellan
    >>>
    >>> _sample_data_dir = 'magellan/test/_sample_data'
    >>>
    >>> magellan.bowtie2_build_index(
    ...     os.path.join(_sample_data_dir, 'sample_bac_genome_short.fa'),
    ...     os.path.join(_sample_data_dir, 'test_index'))
    running:  bowtie2-build magellan/test/_sample_data/sample_bac_genome_short.fa magellan/test/_sample_data/test_index
    ...
    >>> magellan.bowtie2_map(
    ...     os.path.join(_sample_data_dir, 'sample_reads_for_mapper.fastq'),
    ...     os.path.join(_sample_data_dir, 'test_index'),
    ...     os.path.join(_sample_data_dir, 'sample_reads_for_mapper.sam'),
    ...     os.path.join(_sample_data_dir, 'sample_reads_for_mapper.log'))
    bowtie2 -U magellan/test/_sample_data/sample_reads_for_mapper.fastq -x magellan/test/_sample_data/test_index -S magellan/test/_sample_data/sample_reads_for_mapper.sam 2> magellan/test/_sample_data/sample_reads_for_mapper.log
    Bowtie2 version 2.3.4.3;  path: /Users/taliraveh/Downloads/bowtie2-2.3.4.3-macos-x86_64/bowtie2
    >>>
    >>> sam_reads = pd.read_csv(os.path.join(_sample_data_dir, 'sample_reads_for_mapper.sam'), sep='\t', skiprows=5, header=None)
    >>> sam_reads[[0, 2, 3, 4, 5]]
    0                          2    3   4    5
    0              read1_NZ_GL622401.1              NZ_GL622401.1    1  38  80M
    1              read2_NZ_GL622401.1              NZ_GL622401.1   81   1  80M
    2              read3_NZ_GL622401.1              NZ_GL622401.1  161   1  80M
    3  read4_NZ_GL622401.1_fake_break1  NZ_GL622401.1_fake_break1    1  42  80M
    4  read5_NZ_GL622401.1_fake_break1  NZ_GL622401.1_fake_break1   81  42  80M
    5  read6_NZ_GL622401.1_fake_break2  NZ_GL622401.1_fake_break2    1  42  80M
    >>>
    >>>
    >>>
    >>> ap = magellan.AlignmentParser(os.path.join(_sample_data_dir, 'sample_reads_for_mapper.sam'))
    >>>
    >>> ap.count_by_region('NZ_GL622401.1', 1, 100)
    2
    >>>
    >>> read_df = ap.count_all_references()
    >>> read_df
    count  length                       name
    0      3     400              NZ_GL622401.1
    1      2     640  NZ_GL622401.1_fake_break1
    2      1  329244  NZ_GL622401.1_fake_break2
    >>>
    >>> read_df = ap.count_all_references(read_filter_callback=lambda read: magellan.is_high_map_quality(read, 4))
    >>> read_df
    count  length                       name
    0      1     400              NZ_GL622401.1
    1      2     640  NZ_GL622401.1_fake_break1
    2      1  329244  NZ_GL622401.1_fake_break2
    >>>
    >>>
    >>> it = ap.filter(read_filter_callback=lambda read: magellan.is_high_map_quality(read, 39))
    >>> for r in it:
    ...     print r.query_name
    ... 
    read4_NZ_GL622401.1_fake_break1
    read5_NZ_GL622401.1_fake_break1
    read6_NZ_GL622401.1_fake_break2
    >>>
    >>> ap.write_reads(
    ...     os.path.join(_sample_data_dir, 'sample_reads_for_mapper_filtered_mapq4.sam'),
    ...     write_sam=True,
    ...     read_filter_callback=lambda read: magellan.is_high_map_quality(read, 4))
    >>> sam_reads_filtered = pd.read_csv(
    ...     os.path.join(_sample_data_dir, 'sample_reads_for_mapper_filtered_mapq4.sam'),
    ...     sep='\t',
    ...     skiprows=5,
    ...     header=None)
    >>> sam_reads_filtered[[0, 2, 4, 5]]
    0                          2   4    5
    0              read1_NZ_GL622401.1              NZ_GL622401.1  38  80M
    1  read4_NZ_GL622401.1_fake_break1  NZ_GL622401.1_fake_break1  42  80M
    2  read5_NZ_GL622401.1_fake_break1  NZ_GL622401.1_fake_break1  42  80M
    3  read6_NZ_GL622401.1_fake_break2  NZ_GL622401.1_fake_break2  42  80M
    
    '''
    @staticmethod
    def _sort_and_index(in_file, out_file):
        pysam.sort('-o', out_file, in_file)
        pysam.index(out_file)

    def __init__(self, in_file):
        '''
        in_file could be bam or sam format. an extension of si is added to the sorted, indexed bam file that is generated
        '''
        self._base_path, fname = os.path.split(in_file)
        f_base, f_ext = os.path.splitext(fname)
        bam_name = f_base + '_si.bam'
        self._bam_file = os.path.join(self._base_path, bam_name)
        self._sort_and_index(in_file, self._bam_file)

        self._alignment_file = pysam.AlignmentFile(self._bam_file, 'rb')

    def count_by_region(self,
                        ref_name,
                        start,
                        stop,
                        read_filter_callback='nofilter'):
        '''
        Count reads mapped to a specified region (specified by reference name, start and end positions.
        A callback function for filtering reads can be provided with read_filter_callback (default 'nofilter' for no filtering).
        Reads returning True, will be retained
        Returns - integer representing total number of reads mapped to region
        '''
        return self._alignment_file.count(contig=ref_name, start=start, stop=stop, read_callback=read_filter_callback)
    
    def count_all_references(self,
                             read_filter_callback='nofilter'):
        '''
        Counts reads for all references appearing in the alignemnt file across their entire length. Internally calls count_by_region.
        A callback function for filtering reads can be provided with read_filter_callback (default 'nofilter' for no filtering).
        Reads returning True, will be retained
        Returns a pandas dataframe with the following fields: name, length, count
        '''
        counts = []
        for r, l in zip(self._alignment_file.references, self._alignment_file.lengths):
            count_entry = {
                'name': r,
                'length': l,
                'count': self.count_by_region(r, 0, l, read_filter_callback)}
            counts.append(count_entry)
        return pd.DataFrame(counts)

    def filter(self,
               contig=None,
               start=None,
               stop=None,
               read_filter_callback='nofilter'):
        '''
        Creates an iterator over all reads either in a region or in the entire alignemnt file that meet a specific criteria specified by read_filter_callback
        (default 'nofilter' for no filtering).
        '''
        reads_iterator = self._alignment_file.fetch(contig, start, stop)
        for r in reads_iterator:
            if read_filter_callback != 'nofilter' and read_filter_callback(r):
                yield r

    def write_reads(self,
                    out_file,
                    write_sam=False,
                    contig=None,
                    start=None,
                    stop=None,
                    read_filter_callback='nofilter'):
        '''
        write reads to file, optionally using read filtering specified by read_filter_callback (default 'nofilter' for no filtering).
        Can also be used to convert bam to sam format
        '''
        if write_sam:
            filtered_reads = pysam.AlignmentFile(out_file, 'w', template=self._alignment_file)
        else:
            filtered_reads = pysam.AlignmentFile(out_file, 'wb', template=self._alignment_file)
        reads_iterator = self._alignment_file.fetch(contig, start, stop)
        for r in reads_iterator:
            if read_filter_callback == 'nofilter' or read_filter_callback(r):
                filtered_reads.write(r)
        filtered_reads.close()
        

def is_high_map_quality(read, min_qual_threshold=5):
    '''
    Example filter to be used as read_filter_callback when calling AlignmentParser functions.
    Returns True is MAPQ for read is >= min_qual_threshold (default 5), False otherwise
    '''
    return read.mapping_quality >= min_qual_threshold


def has_low_mismatches(read, max_num_mismatches=5):
    '''
    Example filter to be used as read_filter_callback when calling AlignmentParser functions.
    Returns True is number of mismatches for read is <= max_num_mismatches (default 5), False otherwise
    '''
    if not read.has_tag('NM'):  # unaligned
        return False
    return read.get_tag('NM') <= max_num_mismatches
