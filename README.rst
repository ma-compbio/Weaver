Weaver

Allele specific base-pair resolution quantification of Strcutrual variations in cancer genome

Version 0.20


INSTALL
=========


`Bamtools <https://github.com/pezmaster31/bamtools>`_ libraries are needed, included in Weaver_SV/lib and Weaver_SV/inc

`Parallel::ForkManager <http://search.cpan.org/~szabgab/Parallel-ForkManager-1.06/lib/Parallel/ForkManager.pm>`_ perl package is needed

`Bedtools <https://github.com/arq5x/bedtools>`_

`Samtools <http://samtools.sourceforge.net/>`_

`BOOST C++ library <http://www.boost.org/>`_

`BWA <http://bio-bwa.sourceforge.net/>`_

`Bowtie <http://bowtie-bio.sourceforge.net/index.shtml>`_

``export LD_LIBRARY_PATH=<PREFIX>/Weaver/Weaver_SV/lib/:$LD_LIBRARY_PATH``

``libz required //-lz flag``


1	Modify the required BOOST directory in src/Makefile

2	``./INSTALL.sh``



DATA
=========


``wget http://bioen-compbio.bioen.illinois.edu/weaver/Weaver_data.tar.gz``





EXAMPLE DATA
=========

``wget http://bioen-compbio.bioen.illinois.edu/weaver/Weaver_example.tar.gz``



RUN
=========

``Weaver PLOIDY -f SIMU.fa -S FINAL_SV -s SNP -g REGION -w X.bam.wig -r 0 -m map100mer.bd -p 64
solo_ploidy TARGET 2``
``Weaver LITE -f SIMU.fa -S FINAL_SV -s SNP -g REGION -w X.bam.wig -r 0 -m map100mer.bd -p 64 -t 20 -n 0``


Weaver_SV.pl
----------------------------
SV finding

* Input: BAM file from BWA
* Output: VCF file for SV


Weaver_pipeline.pl
----------------------------
Master program to generate SV together with other inputs needed for Weaver

* Input: 1000 Genomes Project Phase 1 haplotypes




Weaver
----------------------------
Core MRF program

* Input: SV
* Outputs:

	1.	Purity and haploid-level sequencing coverage
	2.	Allele specific copy number of genomic regions
	3.	Allele specific copy number of structural variations
	4.	Relative timing of structural variations
	5.	Cancer scaffolds
	6.	Phasing of germline SNPs in CNV regions




Weaver_lite
----------------------------
Core MRF program, with SNP phasing disabled to speed up

* Inputs:

	1.	SV
	2.	reference
	3.	Mappability (available for hg19)
	4.	Region (available for hg19)
	5.	wig (from bam)




Weaver PLOIDY
-------------------------

``Weaver PLOIDY -f  -S  -s ../SNP_dens -g GAP_20140416_num -w  -r 1 -m  -p 16``



* Inputs:

    * -f reference file (fasta), should match the reference used in original bam file. Especially for most TCGA datasets, the alignment was performed on //www.broadinstitute.org/ftp/pub/seq/references/Homo_sapiens_assembly19.fasta, which does not have "chr" prefix  [MANDATORY]
    * -S SV file, with format consistent with Weaver_SV. [MANDATORY]
    * -s SNP file, with ref and alt mappings [MANDATORY]
    * -w wig file from bam, storing the coverage information [MANDATORY]
    * -r 1, if first time running (generating temp files); 0 if want to use existing temp files. [default 1]
    * -m mappability file, download from http://bioen-compbio.bioen.illinois.edu/weaver/Weaver_data.tar.gz [MANDATORY]
    * -p number of cores [default 1]



FILE FORMAT DECLARITIONS
---------------------------

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

