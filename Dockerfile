#==========================================================
# PALMID BASE CONTAINER ===================================
#==========================================================
# Docker Base: amazon linux2
FROM amazonlinux:2 AS serratus-base

## Build/test container for palmid
# sudo yum install -y docker git
# sudo service docker start
# git clone https://github.com/ababaian/palmid.git; cd palmid
#
# sudo docker build -t palmid:local
#

## Push to dockerhub
# sudo docker login
# 
# sudo docker build --no-cache \
#  -t serratusbio/palmid \
#  -t serratusbio/palmid:0.0.3 \
#  -t serratusbio/palmid:latest \
#  -t palmid:latest .
#
# sudo docker push serratusbio/palmid

## Push to ecr (palmid-lambda)
#
# aws ecr-public get-login-password --region us-east-1 \
#  | sudo docker login --username AWS \
#                      --password-stdin public.ecr.aws/q4q7t4w2
#
# docker build -t palmid .
#
# docker tag palmid
#            palmid:0.0.0 \
#            palmid:latest \
#            public.ecr.aws/q4q7t4w2/palmid:latest
# 
# docker push public.ecr.aws/q4q7t4w2/palmid:latest
#

## Dev testing to enter enter
# sudo docker run --rm --entrypoint /bin/bash -it palmid:latest

#==========================================================
# Container Meta-data =====================================
#==========================================================
# Set working directory
# RUN adduser palmid
ENV BASEDIR=/home/palmid
WORKDIR $BASEDIR

# Container Build Information
ARG PROJECT='palmid'
ARG TYPE='base'
ARG VERSION='0.0.3'

# Software Versions (pass to shell)
ENV PALMIDVERSION=$VERSION

ENV SEQKITVERSION='2.0.0'
ENV DIAMONDVERSION='2.0.6-dev'
ENV MUSCLEVERSION='3.8.31'
ENV PALMSCANVERSION='1.0'
ENV PALMDBVERSION='2021-03-14'
ENV R='4'

# Additional Metadata
LABEL author="ababaian"
LABEL container.base.image="amazonlinux:2"
LABEL project.name=${PROJECT}
LABEL project.website="https://github.com/ababaian/palmid"
LABEL container.type=${TYPE}
LABEL container.version=${VERSION}
LABEL container.description="palmid-base image"
LABEL software.license="GPLv3"
LABEL tags="palmscan, diamond, muscle, R, palmid"

#==========================================================
# Dependencies ============================================
#==========================================================
# Update Core
# RUN yum -y update
RUN yum -y install tar wget gzip which sudo shadow-utils \
           util-linux byacc git

# For development
RUN yum -y install vim htop less

# htslib/samtools
RUN yum -y install gcc make \
    unzip bzip2 bzip2-devel xz-devel zlib-devel \
    curl-devel openssl-devel \
    ncurses-devel

# Python3
RUN yum -y install python3 python3-devel &&\
  alias python=python3 &&\
  curl -O https://bootstrap.pypa.io/get-pip.py &&\
  python3 get-pip.py &&\
  rm get-pip.py

# AWS S3
RUN pip install boto3 awscli &&\
  yum -y install jq

# R package dependencies
# PostgreSQL
# Leaflet
# RMarkdown
RUN \
  echo "PostgreSQL" &&\
  yum -y install libxml2-devel postgresql-devel &&\
echo "Leaflet" &&\
  yum install -y gcc-c++.x86_64 cpp.x86_64 sqlite-devel.x86_64 libtiff.x86_64 cmake3.x86_64 &&\
  wget https://download.osgeo.org/geos/geos-3.9.1.tar.bz2 &&\
  tar -xvf geos-3.9.1.tar.bz2 &&\
  cd geos-3.9.1 &&\
  ./configure --libdir=/usr/lib64 &&\
  sudo make &&\
  sudo make install &&\
  wget https://download.osgeo.org/proj/proj-6.1.1.tar.gz &&\
  tar -xvf proj-6.1.1.tar.gz &&\
  cd proj-6.1.1 &&\
  ./configure --libdir=/usr/lib64 &&\
  sudo make &&\
  sudo make install &&\
  cd .. && rm -rf proj-* &&\
  wget https://github.com/OSGeo/gdal/releases/download/v3.2.1/gdal-3.2.1.tar.gz &&\
  tar -xvf gdal-3.2.1.tar.gz &&\
  cd gdal-3.2.1 &&\
  ./configure --libdir=/usr/lib64 --with-proj=/usr/local --with-geos=/usr/local &&\
  sudo make &&\
  sudo make install &&\
  cd .. && rm -rf gdal-* &&\
echo "RMarkdown" &&\
  wget https://github.com/jgm/pandoc/releases/download/2.14.2/pandoc-2.14.2-linux-amd64.tar.gz &&\
  tar xvzf pandoc-2.14.2-linux-amd64.tar.gz --strip-components 1 -C /usr/local &&\
  rm -rf pandoc-2.14.2*

#==========================================================
# Install Software ========================================
#==========================================================

# SeqKit ========================================
RUN wget https://github.com/shenwei356/seqkit/releases/download/v${SEQKITVERSION}/seqkit_linux_amd64.tar.gz &&\
  tar -xvf seqkit* && mv seqkit /usr/local/bin/ &&\
  rm seqkit_linux*

# SAMTOOLS ====================================== (opt?)
## Download
# RUN wget -O /samtools-${SAMTOOLSVERSION}.tar.bz2 \
#   https://github.com/samtools/samtools/releases/download/${SAMTOOLSVERSION}/samtools-${SAMTOOLSVERSION}.tar.bz2 &&\
#   tar xvjf /samtools-${SAMTOOLSVERSION}.tar.bz2 &&\
#   rm /samtools-${SAMTOOLSVERSION}.tar.bz2 &&\
#   cd samtools-${SAMTOOLSVERSION} && make && make install &&\
#   cd .. && rm -rf samtools-${SAMTOOLSVERSION}

# MUSCLE ========================================
RUN wget http://drive5.com/muscle/downloads"$MUSCLEVERSION"/muscle"$MUSCLEVERSION"_i86linux64.tar.gz &&\
  tar -xvf muscle* &&\
  rm muscle*.tar.gz &&\
  mv muscle* /usr/local/bin/muscle

# DIAMOND =======================================
# RUN wget --quiet https://github.com/bbuchfink/diamond/releases/download/v"$DIAMONDVERSION"/diamond-linux64.tar.gz &&\
#   tar -xvf diamond-linux64.tar.gz &&\
#   rm    diamond-linux64.tar.gz &&\
#   mv    diamond /usr/local/bin/

# Use serratus-built dev version
RUN wget --quiet https://serratus-public.s3.amazonaws.com/bin/diamond &&\
    chmod 755 diamond &&\
    mv    diamond /usr/local/bin/

# PALMSCAN ======================================
RUN wget -O /usr/local/bin/palmscan \
  https://github.com/ababaian/palmscan/releases/download/v${PALMSCANVERSION}/palmscan-v${PALMSCANVERSION} &&\
  chmod 755 /usr/local/bin/palmscan

# PALMDB ========================================
# clone repo + make sOTU-database
RUN git clone https://github.com/rcedgar/palmdb.git &&\
  gzip -dr palmdb/* &&\
  cp "palmdb/"$PALMDBVERSION"/otu_centroids.fa" palmdb/palmdb.fa &&\
  diamond makedb --in palmdb/palmdb.fa -d palmdb/palmdb
  # db hash: 0c43dc6647b7ba99b4035bc1b1abf746

# R 4.0 =========================================
# Install R
# Note: 1 GB install
RUN amazon-linux-extras install R4

# R Packages ====================================
RUN \
  R -e 'install.packages( c("devtools"), repos = "http://cran.us.r-project.org")' &&\
  R -e 'library("devtools"); install_github("ababaian/palmid")'

#==========================================================
# palmid Initialize =======================================
#==========================================================
# scripts + test data
COPY palmid.Rmd scripts/* ./
COPY data/* inst/extdata/* img/* data/
RUN chmod 755 palmid.sh &&\
    chmod 755 fev2tsv.py

#==========================================================
# Resource Files ==========================================
#==========================================================
# Sequence resources / databases for analysis
# RUN cd /home/palmid/ &&\
#   git clone https://github.com/rcedgar/palmdb.git &&\
#   gzip -dr palmdb/*

#==========================================================
# CMD =====================================================
#==========================================================
CMD ["/home/palmid/palmid.sh"]
