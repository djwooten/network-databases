import setuptools

import codecs
import os.path

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

with open("README.md", "r") as fh:
    long_description = fh.read()

def get_packages():
    packages=setuptools.find_packages(where='src')
    return packages

setuptools.setup(
    name="networkdb", # Replace with your own username
    version=get_version("src/networkdb/__init__.py"),
    author="David J. Wooten",
    author_email="dwooten@psu.edu",
    description="Python package interacting with network databases and gene names",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/djwooten/networkdb",

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
    ],
    python_requires='>=3.5',
    keywords='systems biology network gene regulatory',
    packages=get_packages(),
    package_dir={'': 'src'},
    install_requires=[
        "networkx >= 2.0",
        "pandas >= 1.0.0"
    ],
    # package_data VS data_files VS ???
)
