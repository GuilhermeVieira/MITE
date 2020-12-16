FROM python:3.8.3-buster

WORKDIR /apps

# Installing SuperHirn

RUN apt-get update && \
    apt-get -y install git && \
    git clone https://github.com/GuilhermeVieira/SuperHirn.git && \
    cd SuperHirn/SuperHirnv03/make/ && \
    make

# Installing OpenMPI

RUN mkdir openmpi && \
    cd openmpi && \
    wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.3.tar.gz && \
    tar -xzvf openmpi-4.0.3.tar.gz && \
    cd openmpi-4.0.3/ && \
    ./configure --prefix=/apps/openmpi && \
    make all && \
    make install

ENV PATH "$PATH:/apps/openmpi/bin"

ENV LD_LIBRARY_PATH "$LD_LIBRARY_PATH:/apps/openmpi/lib"

# Setting environment variables that enables running ompi as root

ENV OMPI_ALLOW_RUN_AS_ROOT 1

ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM 1

# Installing MrBayes

RUN mkdir mrbayes && \
    cd mrbayes && \
    wget https://github.com/NBISweden/MrBayes/releases/download/v3.2.7/mrbayes-3.2.7.tar.gz && \
    tar -xzvf mrbayes-3.2.7.tar.gz && \
    cd mrbayes-3.2.7/ && \
    ./configure --prefix=/apps/mrbayes --with-mpi && \
    make && \
    make install

ENV PATH "$PATH:/apps/mrbayes/bin"

# Installing R

RUN echo "deb [trusted=yes] https://eddelbuettel.github.io/drr35/ ./" > /etc/apt/sources.list.d/debian-r-3.5.list 

ENV R_BASE_VERSION 3.5.2

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    littler \
    r-cran-littler \
    r-base=${R_BASE_VERSION}-* \
    r-base-dev=${R_BASE_VERSION}-* \
    r-recommended=${R_BASE_VERSION}-* \
    && echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site \
    && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r \
    && ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
    && ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
    && ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
    && ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
    && install.r docopt \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*

# Installing R packages

RUN R -e "install.packages('ape', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('pvclust', dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Installing Python packages

COPY requirements.txt ./

RUN pip install --upgrade pip && \
    pip install -r requirements.txt

WORKDIR /user/src/mite

ENTRYPOINT ["bash"]