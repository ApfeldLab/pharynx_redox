FROM python:3.7-slim-buster

# Install some necessary dependencies
RUN apt-get -q update && \
    apt-get install -q -y --no-install-recommends \
        unzip \
        wget \
        curl

# Install the MATLAB runtime
RUN mkdir /mcr-install && \
    mkdir /opt/mcr && \
    cd /mcr-install && \
    wget -q http://ssd.mathworks.com/supportfiles/downloads/R2019a/Release/4/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2019a_Update_4_glnxa64.zip && \
    unzip -q MATLAB_Runtime_R2019a_Update_4_glnxa64.zip && \
    rm -f MATLAB_Runtime_R2019a_Update_4_glnxa64.zip && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && \
    rm -rf mcr-install

# Configure MATLAB environment variables
ENV LD_LIBRARY_PATH /opt/mcr/v96/runtime/glnxa64:/opt/mcr/v96/bin/glnxa64:/opt/mcr/v96/sys/os/glnxa64:/opt/mcr/v96/extern/bin/glnxa64