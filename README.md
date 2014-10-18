License:
--------
ibdadmx.py
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

Program available for download at: 
----------------------------------
https://github.com/itamareskin/IBDAdmixed

About:
------
ibadmx.py is a python program for detecting segments of identity-by-descent (IBD) on datasets of admixed individuals. The program was developed by Itamar (Eskin) Afek at Tel Aviv University. Full description of the algorithm is included in the paper "IBDAdmixed: Fine-Scale Detection of Identity-by-Descent in Admixed Populations".

Dependencies:
-------------

Python 2.6
Cython

Installation:
-------------

After downloading the source code, navigate to the lib subdirectory and run the following in a unix command line:
python setup.py build_ext --inplace

Usage:
------

For detailed usage instruction, run
python ibadmx.py --help

Required Input:
---------------

input         prefix of input files (.dat,.map in IBDAdmixed/LAMP format)
out           output prefix
hapmodelfile  .dag beagle model file name (one for each ancestry)

Options:
--------

  -h, --help            show this help message and exit
  -k K, --num-anc K     set number of ancestries
  -a ALPHAS [ALPHAS ...], --set-alphas ALPHAS [ALPHAS ...]
                        set the alphas
  -g GENERATIONS, --generations GENERATIONS
                        set the number of generations of admixture
  -n NUM_SNPS, --max-snps NUM_SNPS
                        maximal number of snps to be used
  -p NUM_CPUS, --num-cpus NUM_CPUS
                        number of cpus to be used
  -e EPSILON, --epsilon EPSILON
                        epsilon for error
  -m MIN_SCORE, --min-score MIN_SCORE
                        minimal score to report as IBD
  -w WIN_SIZE, --win-size WIN_SIZE
                        window size (in number of SNPs)
  -o OFFSET, --offset OFFSET
                        offset for windows (in number of SNPs)
  --pair PAIR PAIR      single pair to process
  --pairs-file PAIRS_FILE
                        file containing pairs of individuals to process
  --germline-file GERMLINEFILE
                        germline results file
  --set-ibd-trans IBD_TRANS [IBD_TRANS ...]
                        set ibd to no IBD probabilities
  --debug               print debugging information
  --phased              use phased mode
  --naive-model         use naive model
  --scramble            scramble phase of genotypes
  --condor              send jobs to htcondor
  --keep-temp           keep temoporary files
  --rerun               rerun failed jobs
  --recover             recover results from all finished jobs

Testing:
--------

