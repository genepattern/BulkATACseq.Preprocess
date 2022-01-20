# copyright 2017-2021 Regents of the University of California and the Broad Institute. All rights reserved.
#version 1.0
FROM broadinstitute/picard:2.25.4

MAINTAINER Edwin Juarez <ejuarez@ucsd.edu>
ENV LANG=C LC_ALL=C

USER root

# Making alias for running picard
RUN echo 'alias picard="java -jar /usr/picard/picard.jar"' >> ~/.bashrc

# Installing Samtoolsa
## prereqs
RUN apt-get update  # Ensure the package list is up to date 2021-05-17 at 14:03:07
RUN apt-get -y install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev
## Download latest version of samtools 
RUN apt-get install wget
RUN wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2
RUN tar -vxjf samtools-1.12.tar.bz2 &&\
    cd samtools-1.12 &&\
    ./configure &&\
    make &&\
    make install

# Installing sambamba
RUN mkdir /usr/sambamba
RUN cd /usr/sambamba/
RUN wget https://github.com/biod/sambamba/releases/download/v0.8.0/sambamba-0.8.0-linux-amd64-static.gz
RUN gzip -d sambamba-0.8.0-linux-amd64-static.gz
RUN chmod u+x sambamba-0.8.0-linux-amd64-static
RUN mv sambamba-0.8.0-linux-amd64-static /usr/sambamba/
RUN echo 'alias sambamba="/usr/sambamba/sambamba-0.8.0-linux-amd64-static"' >> ~/.bashrc

# Assign Multimappers Python File used for BulkATAC Seq Preprocessing
# ADD assign_multimappers.py /
RUN mkdir /genepattern
# run python 3 launch docker locally, ask what version of python, otherwise add python3