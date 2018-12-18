FROM ubuntu:16.04

WORKDIR /opt

# Installing softwares 
RUN apt-get update && apt-get install -y \
		wget \
		make \
		g++ \
		build-essential \
		autoconf \
		automake \
		zlib1g-dev \
		libgsl0-dev \
		python \
		perl \
        cmake \
		git \
		curl \
		libncurses5-dev \
		unzip


#RUN wget https://cmake.org/files/v3.9/cmake-3.9.0-rc5.tar.gz
#RUN tar -xf cmake-3.9.0-rc5.tar.gz
#WORKDIR /opt/cmake-3.9.0-rc5
#RUN ./bootstrap.sh && make && make install
#ENV PATH="/opt/cmake-3.9.0-rc5:$PATH"
 

RUN wget https://github.com/pezmaster31/bamtools/archive/v2.4.1.tar.gz
RUN tar -zxvf v2.4.1.tar.gz
WORKDIR /opt/bamtools-2.4.1
RUN mkdir build 
WORKDIR /opt/bamtools-2.4.1/build
RUN cmake .. 
RUN make
ENV LD_LIBRARY_PATH="/opt/bamtools-2.4.1/lib:${LD_LIBRARY_PATH}"

# Installing perl Parallel:ForkManager
WORKDIR /opt
RUN wget http://search.cpan.org/CPAN/authors/id/Y/YA/YANICK/Parallel-ForkManager-1.19.tar.gz
RUN tar -zxvf Parallel-ForkManager-1.19.tar.gz
WORKDIR /opt/Parallel-ForkManager-1.19
RUN perl Makefile.PL
RUN make && make install

# Installing bedtools 
WORKDIR /opt
RUN wget https://github.com/arq5x/bedtools2/archive/v2.26.0.tar.gz
RUN tar -zxvf v2.26.0.tar.gz
WORKDIR /opt/bedtools2-2.26.0
RUN make 

#Installing samtools 
WORKDIR /opt
RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
RUN tar -xjf samtools-1.3.1.tar.bz2
WORKDIR /opt/samtools-1.3.1
RUN make && make install

# Installing boost
WORKDIR /opt
RUN wget https://sourceforge.net/projects/boost/files/boost/1.63.0/boost_1_63_0.tar.gz
RUN tar -zxvf boost_1_63_0.tar.gz
WORKDIR /opt/boost_1_63_0
RUN ./bootstrap.sh --with-libraries=program_options,chrono,filesystem,system 
RUN ./b2 install

# Downloading Weaver
WORKDIR /opt
RUN echo "reclone" & git clone https://github.com/ma-compbio/Weaver.git
WORKDIR /opt/Weaver/src
#RUN sed -i.bak "s/BOOST = \/usr0\/home\/ashokr\/local\/boost_1_61_0\//BOOST = \/opt\/boost_1_63_0\//" Makefile 
#RUN sed -i.bak "s/BOOST_OPT = \/usr0\/home\/ashokr\/local\/boost_1_61_0\/bin\.v2\/libs\/libboost_program_options.a/BOOST_OPT = \/opt\/boost_1_63_0\/bin.v2\/libs\/program_options\/build\/*\/release\/link-static\/threading-multi\/libboost_program_options.a/" Makefile

# Commenting off partition.cpp parallelization
RUN sed -i.bak "s/#pragma /\/\/#pragma/" partition.cpp
#RUN sed -i.bak "s/\/\/ofstream otemp_1(\"EACH_REGION_1\")/ofstream otemp_1(\"EACH_REGION_1\")/g" LBP.cpp

# Small changes to make Weaver upstream pipeline run
WORKDIR /opt/Weaver/bin
#RUN sed -i.bak "s/\$TorN || \"T\"/\$TorN = \$TorN || \"T\";/g" Weaver_pipeline.pl

# Downloading necessary data
WORKDIR /opt/Weaver/data
RUN wget http://genome.compbio.cs.cmu.edu/~ashokr/data/Weaver_data.tar.gz
RUN tar -xvzf Weaver_data.tar.gz
RUN mv Weaver_data/* .
RUN rm -r Weaver_data

## Renaming old binaries
#WORKDIR /opt/Weaver/external_bin
#RUN mv /opt/Weaver/external_bin/bwa /opt/Weaver/external_bin/bwa_old
#RUN mv /opt/Weaver/external_bin/bowtie /opt/Weaver/external_bin/bowtie_old
#
## Installing bwa 
#WORKDIR /opt
#RUN 

http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.15.tar.bz2
#RUN tar -xjf bwa-0.7.15.tar.bz2
#WORKDIR /opt/bwa-0.7.15
#RUN make
## Linking from Weaver/external_bin
#RUN ln -sf /opt/bwa-0.7.15/bwa /opt/Weaver/external_bin/bwa
#
## Downloading bowtie 
#WORKDIR /opt
#RUN curl -L "http://freefr.dl.sourceforge.net/project/bowtie-bio/bowtie2/2.2.6/bowtie2-2.2.6-linux-x86_64.zip" > /opt/bowtie2-2.2.6-linux-x86_64.zip
#RUN unzip bowtie2-2.2.6-linux-x86_64.zip
#WORKDIR /opt/bowtie2-2.2.6
## Linking from Weaver/external_bin
#RUN ln -sf /opt/bowtie2-2.2.6/bowtie /opt/Weaver/external_bin/bowtie
# CHMODS
WORKDIR /opt/Weaver/external_bin
RUN chmod 777 *

# Changing environment variable
WORKDIR /opt/
ENV PATH="/opt/Weaver/bin/:$PATH"
ENV LIBRARY_PATH="/opt/boost_1_63_0/stage/lib:${LIBRARY_PATH}"
ENV INCLUDE="/opt/boost_1_63_0:$INCLUDE"
ENV LD_LIBRARY_PATH="/opt/Weaver/Weaver_SV/lib/:${LD_LIBRARY_PATH}"
#ENV PATH="/opt/bwa-0.7.15/:/opt/bowtie2-2.2.6/:$PATH"


# Installing Weaver
WORKDIR /opt/Weaver
RUN RUNARG=0 
RUN bash INSTALL.sh /opt/boost_1_63_0 /opt/boost_1_63_0/bin.v2/libs/program_options/build/*/release/link-static/threading-multi/libboost_program_options.a



WORKDIR /opt
RUN rm *.tar.gz *.bz2 #*.zip
#RUN PATH=/opt/Weaver/bin:${PATH}
#RUN PATH=/opt/Weaver/Weaver_SV/bin:${PATH}

COPY Dockerfile /opt
MAINTAINER ashokr <ashokr@cs.cmu.edu>
