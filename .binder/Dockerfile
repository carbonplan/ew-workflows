FROM quay.io/igwn/builder:bookworm    
# FROM quay.io/pangeo/pangeo-notebook

# -- troubleshoot: 
#    build image with `docker build -t scepter-dev:latest .` from the same dir as the dockerfile
#    if running into issues, you can remove the error lines and get into the container 
#    to see what's in there with `docker run -it --rm scepter-dev:latest /bin/bash`. exit with `exit`.

USER root

# Install system dependencies
RUN apt-get update \
    && apt-get -y install \
        python3 \ 
        python3-pip \
        python3-venv \
        make \
        gfortran \
        curl \
        wget \
        libopenblas-dev \
    && export LD_LIBRARY_PATH=/usr/lib/openblas-base/

# Install Jupyter and typical notebook deps
RUN apt-get update && apt-get install -y \
    &&  python3-jupyterlab \
        python3-ipykernel \
        python3-numpy \
        python3-pandas \
        python3-dask distributed

# install s5cmd (latest release from GitHub)
RUN S5CMD_VERSION=2.3.0 \
    && wget -q https://github.com/peak/s5cmd/releases/download/v${S5CMD_VERSION}/s5cmd_${S5CMD_VERSION}_Linux-64bit.tar.gz \
    && tar -xzf s5cmd_${S5CMD_VERSION}_Linux-64bit.tar.gz -C /usr/local/bin s5cmd \
    && rm s5cmd_${S5CMD_VERSION}_Linux-64bit.tar.gz \
    && chmod +x /usr/local/bin/s5cmd

# # Clone SCEPTER 
# # (first switch user, required for quay.io/pangeo/pangeo-notebook)
# USER jovyan    
# WORKDIR /home/jovyan
# RUN git clone https://github.com/carbonplan/SCEPTER.git 

# # Update the Makefile with the proper SCEPTER data path
# WORKDIR /home/jovyan/SCEPTER
# RUN sed -i "s|-I/home/tykukla/SCEPTER/data|-I${PWD}/data|" makefile

# # Use root (default)
# USER root
# WORKDIR /root

# # Clone SCEPTER
# RUN git clone https://github.com/carbonplan/SCEPTER.git

# # Update Makefile
# WORKDIR /root/SCEPTER
# RUN sed -i "s|-I/home/tykukla/SCEPTER/data|-I${PWD}/data|" makefile

# # compile
# RUN make clean
# RUN make all
# RUN make rateA
# RUN make rateB
# RUN make rateC
# RUN make rateD
# RUN make richards