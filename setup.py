from setuptools import setup, find_packages

setup(
    name='metaDMM',
    version='1.0.0',
    description='A tool for generating synthetic DNA reads from complete genomes',
    author='Swati Tak Nielsen',
    author_email='swtak@proton.me',
    url='https://github.com/aswatib/metaDMM',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'pyfaidx',
        'biopython',
        'click',
	'pandas',
    ],
    entry_points={
        'console_scripts': [
            'metaDMM = metaDMM.cli:cli',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
