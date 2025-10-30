FROM quay.io/pangeo/pangeo-notebook

# -- troubleshoot: 
#    build image with `docker build -t scepter-dev:latest .` from the same dir as the dockerfile
#    if running into issues, you can remove the error lines and get into the container 
#    to see what's in there with `docker run -it --rm scepter-dev:latest /bin/bash`. exit with `exit`.

USER root

# Install system dependencies
RUN apt-get update \
    && apt-get -y install \
        make \
        gfortran \
        curl \
        wget \
        libopenblas-dev \
    && export LD_LIBRARY_PATH=/usr/lib/openblas-base/

# install s5cmd (latest release from GitHub)
RUN wget -q https://github.com/peak/s5cmd/releases/latest/download/s5cmd_Linux-64bit.tar.gz \
    && tar -xzf s5cmd_Linux-64bit.tar.gz -C /usr/local/bin s5cmd \
    && rm s5cmd_Linux-64bit.tar.gz \
    && chmod +x /usr/local/bin/s5cmd

# Clone SCEPTER 
# (first switch user, required for quay.io/pangeo/pangeo-notebook)
USER jovyan    
WORKDIR /home/jovyan
RUN git clone https://github.com/carbonplan/SCEPTER.git 

# Update the Makefile with the proper SCEPTER data path
WORKDIR /home/jovyan/SCEPTER
RUN sed -i "s|-I/home/tykukla/SCEPTER/data|-I${PWD}/data|" makefile

# compile
RUN make clean
RUN make all
RUN make rateA
RUN make rateB
RUN make rateC
RUN make rateD