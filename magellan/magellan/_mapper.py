import os
import subprocess

# Do we want to install bowtie? do we want to adapt to windows?
# use logger?


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
    cmd = 'bowtie2-build %s %s' % (fasta_file, index)
    if offrate:
        cmd = cmd + '-o ' + offrate
    print 'running: ', cmd
    _check_bowtie_install()
    
    '''
    _, err = subprocess.Popen(cmd,
                              stderr=subprocess.PIPE,
                              shell=True).communicate()
    if err:
        raise RuntimeError('bowtie index building failed for fasta %s' % (fasta_file))
    '''
    

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

    
