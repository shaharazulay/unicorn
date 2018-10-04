import os
import sys
from setuptools import setup, Command
from setuptools.command.install import install
import zipfile
import subprocess

try:
    from urllib.request import urlretrieve
except ImportError:
    from urllib import urlretrieve

_python = 'python%d' % sys.version_info.major


def find_install(software_name):
    '''
    Check if software is installed.
    '''
    try:
        software_path = subprocess.Popen(['which', software_name], stdout=subprocess.PIPE).communicate()[0].decode()
        if len(software_path) == 0:
            return None
        if not os.access(software_path, os.X_OK):
            return None
        print 'Found software at: ' + software_path
        return software_path
    except Exception:
        return None
    

def find_exe_in_path(software_name):
    """
    Check that an executable exists in $PATH
    """
    
    paths = os.environ["PATH"].split(os.pathsep)
    for path in paths:
        check_exe_full_path = os.path.join(path, software_name)
        if os.path.exists(check_exe_full_path):
            if os.access(check_exe_full_path, os.X_OK):
                return path
    return None


def download_unpack_zip(url, dest_folder, filename, software_name):
    """
    Download the url to the file and decompress into the folder
    """
    
    # Check for write permission to the target folder
    if not os.access(dest_folder, os.W_OK):
        print 'WARNING: Destination directory is not writeable: ' + dest_folder

    print 'downloading ' + url
    urlretrieve(url, os.path.join(dest_folder, filename))
    
    try:
        zipfile_handle = zipfile.ZipFile(os.path.join(dest_folder, filename))
        zipfile_handle.extractall(path=dest_folder)
        zipfile_handle.close()
        return True
    except EnvironmentError:
        print 'WARNING: Unable to extract ' + software_name
        return False
        

def install_bowtie2(base_folder=None, mac_os=False):
    """
    Download and install the bowtie2 software if not already installed. Based on HumanN code at:
    https://bitbucket.org/biobakery/humann2/src/76e92a1250f0265dff390ba903ba22b8a74b6c8e/setup.py?at=default&fileviewer=file-view-default
    """

    bowtie2_exe = 'bowtie2'
    bowtie2_index_build = 'bowtie2-build'
    
    # Check if bowtie2 is already installed
    bowtie2_path = find_install(bowtie2_exe)  # find_exe_in_path(bowtie2_exe)

    if not bowtie2_path:
        bowtie2_url = 'https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.3/bowtie2-2.3.4.3-linux-x86_64.zip/download'
        if mac_os:
            bowtie2_url = 'https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.3/bowtie2-2.3.4.3-macos-x86_64.zip/download'

        if not base_folder:
            base_folder = os.getcwd()
        bowtie2_install_folder = os.path.join(base_folder, 'bowtie2-2.3.4.3')
        bowtie2_download_file = os.path.split(os.path.split(bowtie2_url)[0])[1]
        print 'DEBUG: Did not find BOWTIE2, This is found url - ' + bowtie2_url + '. This is download dir - ' + bowtie2_install_folder

        if not os.path.exists(bowtie2_install_folder):
            os.mkdir(bowtie2_install_folder)
        
        # install the bowtie2 software
        print('Bowtie2 missing, installing')
        install_error = False
        exe_dir = os.path.join(bowtie2_install_folder, os.path.splitext(bowtie2_download_file)[0])
        print 'This is the expected exe dir ' + exe_dir

        download_ok = download_unpack_zip(bowtie2_url, bowtie2_install_folder, bowtie2_download_file, bowtie2_exe)
        
        if not download_ok:
            install_error = True
        else:
            os.environ["PATH"] += os.pathsep + exe_dir

            files = []
            try:
                files = os.listdir(exe_dir)
                if bowtie2_exe not in files:
                    raise EnvironmentError('Bowtie2 exe not found')
                os.chmod(os.path.join(exe_dir, bowtie2_exe), 0o755)
                os.chmod(os.path.join(exe_dir, bowtie2_index_build), 0o755)
                
            except EnvironmentError:
                print("WARNING: Bowtie2 files not found.")
                install_error = True
                
        if install_error:
            print 'WARNING: Unable to install bowtie2.'
        else:
            print 'Installed bowtie2 at ' + exe_dir
        
    else:
        print 'Found bowtie2 at ' + bowtie2_path


class _InstallCommand(install):
    def run(self):
        # find out the platform
        mac_os = False
        if sys.platform in ['darwin', 'os2', 'os2emx']:
            mac_os = True
        install_bowtie2(mac_os=mac_os)
        install.run(self)

        
class _TestCommand(Command):
    user_options = [
        ]

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        install_bowtie2()
        run_str = "%s -m unittest discover test *test.py" % _python
        os.system(run_str)


setup(
    name='magellan',
    version='0.0.1',
    author='Shahar Azulay, Tali Raveh, Ariel Hanemann, Yossi Cohen',
    author_email='shahar4@gmail.com',
    url='https://github.com/shaharazulay/unicorn/magallen',
    packages=[
        'magellan'
    ],
    license='bsd',
    description='Detection of gene-level deletions, to uncover new bacterial strains',
    long_description=open('docs/README.rst').read(),
    install_requires=[],
    zip_safe=False,
    package_data={},
    include_package_data=True,
    cmdclass={
        'test': _TestCommand,
        'install': _InstallCommand
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Information Analysis',
    ],
)
