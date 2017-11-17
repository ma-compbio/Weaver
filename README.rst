===================================================================================
Weaver -- Allele-Specific Quantification of Structural Variations in Cancer Genomes
===================================================================================
Version 0.21

Overview
--------
Weaver is a tool which, given a whole genome sequencing sample from a tumour
sample (in BAM format) and reference FASTA file as input, returns
allele-specific copy numbers of regions and phased breakpoints of structural
variants. The principal model in the framework is a Markov Random Field, 
which takes as auxiliary inputs germline SNP data from the 1000 Genomes
database, GC content and mappability, and has the phase and copy number of
genomic loci as hidden states.

Installation
------------

Installing Weaver requires the following dependencies.

    1. `CMake <https://cmake.org>`_
    2. `Bamtools <https://github.com/pezmaster31/bamtools>`_ libraries are needed, included in Weaver_SV/lib and Weaver_SV/inc
    3. `Parallel::ForkManager <http://search.cpan.org/~szabgab/Parallel-ForkManager-1.06/lib/Parallel/ForkManager.pm>`_ perl package is needed
    4. `Bedtools <https://github.com/arq5x/bedtools>`_
    5. `Samtools <http://samtools.sourceforge.net/>`_
    6. `BOOST C++ library <http://www.boost.org/>`_
    7. `BWA <http://bio-bwa.sourceforge.net/>`_
    8. `Bowtie <http://bowtie-bio.sourceforge.net/index.shtml>`_
    9. `libz required //-lz flag`

We recommended that the system on which Weaver is installed have more than
40 GB of memory, since Weaver incorporates variant calling. For small inputs,
depending on the number of process threads, the memory used can be much lower.

In order to install Weaver, we need to run the following commands.
::

    ``export LD_LIBRARY_PATH=<PREFIX>/Weaver/Weaver_SV/lib/:$LD_LIBRARY_PATH``

Define variables ``$BOOST`` and ``${BOOST_OPT}`` as the locations of your Boost
install and the linker file for the program options library in boost. The
latter is traditionally located at:: 

``$BOOST/bin.v2/libs/program_options/build/gcc-4.8/release/link-static/threading-multi/libboost_program_options.a``

Then run the following command.
::

    ``./INSTALL.sh $BOOST ${BOOST_OPT}``

The Weaver executable will be located in ``bin/`` within the installation directory.

Auxiliary data
++++++++++++++

Weaver requires input data that is available for download.
::

    ``wget http://bioen-compbio.bioen.illinois.edu/weaver/Weaver_data.tar.gz``

The data must be stored in the folder ``data/`` in the installation directory.



Example data
++++++++++++
::

    ``wget http://bioen-compbio.bioen.illinois.edu/weaver/Weaver_example.tar.gz``


Running Weaver
--------------

We will assume the following input variables, which we need for all the
scripts.

    1. ``${BAM}``: Sorted and indexed BAM file. Index must be present in same
       directory.
    2. ``${REFDIR}``: The address and prefix to the reference FASTA file.
       Exclude the extension (``.fa`` or ``.fasta``). The location must also
       contain the Bowtie (``*.ebwt``) and BWA indices, with the same prefix.
       For the sake of exposition, we will assume that the reference FASTA
       file has extension ``.fa``.
    3. ``${GAPALPHA}``: The file ``GAP_20140416`` in the ``data/`` directory.
    4. ``${GAP}``: The file ``GAP_20140416_num`` in the ``data/`` directory.
       This is the same as ``$GAPALPHA``, but does not have the ``chr`` prefix
       preceding chromosome names.
    5. ``${1000G}``: The 1000 Genomes Phase 1 haplotypes data. This can be 
       obtained from the 1000 Genomes project website.
    6. ``${SEX}``: The sex of the subject from whom the sample is obtained.
    7. ``${MAP}``: The file ``wgEncodeCrgMapabilityAlign100mer_number.bd`` in
         the ``data/`` directory.
    8. ``${THREADS}``: User-specified number of threads the program will run
       on.

We will use ``${BIN}`` as the default directory under which the Weaver
executables are stored. This will be the full path to the directory ``bin/``
in the installation.

In order to run Weaver, we must first call the different single nucleotide and
structural variants in the sample. To do this, we will use the perl script  
``bin/Weaver_pipeline.pl``.

Using Weaver_pipeline.pl
++++++++++++++++++++++++
::

    Usage:
            $BIN/Weaver_pipeline.pl ALL <mode> \ 
                                -p/--thread     number of cores
                                -f/--fa         [MANDATORY] bowtie and bwa reference dir/name
                                -g/--gap        [MANDATORY] Gap file 
                                -b/--bam        bam file
                                -o/--output     output dir
                                -k/--onekg      1000 Gemomes Project data dir
                                -s/--sex        Female (F) or Male (M). Y chromosome will not be used if the bam is from female tissue.
                                -h/--help

``<mode>`` takes one of the following arguments: ``SV, SNP, WIG``.

Calling SVs
+++++++++++
::

      perl $BIN/Weaver_pipeline.pl ALL SV \
                               -f ${REFDIR} \
                               -g ${GAPALPHA} \
                               -b ${BAM} \
                               -k ${T1000} \
                               -s ${SEX} \
                               -p ${THREADS}

Calling SNVs
++++++++++++
::

      perl $BIN/Weaver_pipeline.pl ALL SNP \
                               -f ${REFDIR} \
                               -g ${GAPALPHA} \
                               -b ${BAM} \
                               -k ${T1000} \
                               -s ${SEX} \
                               -p ${THREADS}

Creating WIG file
+++++++++++++++++
::

      perl $BIN/Weaver_pipeline.pl ALL WIG \
                               -f ${REFDIR} \
                               -g ${GAPALPHA} \
                               -b ${BAM} \
                               -k ${T1000} \
                               -s ${SEX}

Finding the haplotype level coverage
++++++++++++++++++++++++++++++++++++

The core Weaver program needs haplotype level coverage for the cancer and
normal genomes as input. We can estimate this using the following command
from the same directory that ``Weaver_pipline.pl`` was executed. Assume
that the variable ``${NEWGAP}`` is equal to ``$GAPALPHA`` if the reference 
FASTA and BAM file have chromosome names with ``chr`` prefixed, and equal to
``$GAP`` otherwise.
::

    $BIN/Weaver PLOIDY -f ${REFDIR}.fa \
                       -S ${BAM}.Weaver.GOOD \
                       -s SNP_dens \
                       -g ${NEWGAP} \
                       -w ${BAM}.wig \
                       -z ${TILESIZE} \
                       -r 1 \
                       -m $MAP \
                       -p $THREADS

* Inputs:

    * -f reference file (fasta), should match the reference used in original bam file. Especially for most TCGA datasets, the alignment was performed on //www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta, which does not have "chr" prefix  [MANDATORY]
    * -S SV file, with format consistent with Weaver_SV. [MANDATORY]
    * -s SNP file, with ref and alt mappings [MANDATORY]
    * -w wig file from bam, storing the coverage information [MANDATORY]
    * -z tile size argument: positive integer size to partition genome. [default 5000]
    * -r 1, if first time running (generating temp files); 0 if want to use existing temp files. [default 1]
    * -m mappability file, download from http://bioen-compbio.bioen.illinois.edu/weaver/Weaver_data.tar.gz [MANDATORY]
    * -p number of cores [default 1]

* Output:
    * TARGET: File containing haplotype level coverage of different regions

The tile size argument must be set by trial and error, so that the ``TARGETi`` file is properly populated. The argument varies from sample to sample, but the usual range is from 500 to 5000, with 1000 being common for TCGA samples.
Once this is obtained, we use the following command to obtain the haplotype level coverage.
:: 

    $BIN/solo_ploidy TARGET 2

The ``2`` here indicates a diploid normal genome. This will write the estimated
haplotype level normal and tumour coverage to ``STDOUT``.

Run Weaver core program
+++++++++++++++++++++++

Finally, in order to obtain the main result, we run the following script. Here,
we assume that ``${TUMOUR_COV}`` and ``${NORMAL_COV}`` are the tumour and
normal haplotype level coverage obtained in the previous step respectively.
::

    $BIN/Weaver LITE -f ${REFDIR}.fa \
                     -S ${BAM}.Weaver.GOOD \
                     -s SNP_dens \
                     -g ${NEWGAP} \
                     -w ${BAM}.wig \
                     -z ${TILESIZE} \
                     -r 1 \
                     -m $MAP \
                     -t ${TUMOUR_COV} \
                     -n ${NORMAL_COV} \
                     -p $THREADS



File format declaritions
------------------------

Wiggle file
+++++++++++

Wiggle file need to be declared with fixedStep, step 1 and span 1
fixedStep chrom=chr1 start=9994 step=1 span=1
if a chromosome has multiple declaration lines, they need to be sorted based on position:
fixedStep chrom=chr1 start=9994 step=1 span=1
X
X
X
fixedStep chrom=chr1 start=100 step=1 span=1
X
X
X
Is not allowed



Bam file
+++++++++

Must be sorted and indexed.

SNP file:

NGS SNP link file


1KGP SNP link


SV
++++++


Genome region file:

GAP regions in assembly are annotated.


Output:
=======

REGION_CN_PHASE
+++++++++++++++
Storing phased allele specific copy number of genome

CHR	BEGIN	END	ALLELE_1_CN	ALLELE_2_CN




SV_CN_PHASE
+++++++++++

Structural variation copy number and phasing, catagory

CHR_1	POS_1	ORI_1	ALLELE_	CHR_2   POS_2   ORI_2   ALLELE_	CN	germline/somatic_post_aneuploidy/somatic_pre_aneuploidy


CONTACT
=======

`Yang Li <leofountain@gmail.com>`_
Jian Ma's Computational Genomics Lab at Carnegie Mellon
The code was developed by Yang Li when the Ma lab was at the University of Illinois at Urbana-Champaign

https://github.com/ma-compbio/Weaver

