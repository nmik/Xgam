FROM ubuntu:18.04
MAINTAINER simone.ammazzalorso@unito.it

# Install needed python libraries
RUN apt-get update && yes|apt-get upgrade
RUN apt-get update && apt-get install -y curl git time
WORKDIR tmp
RUN curl -O https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
RUN mkdir /archive/ && mkdir /archive/home/ && mkdir /archive/home/sammazza/ && mkdir /archive/home/sammazza/fermi_data/
RUN mkdir /run_xgam/

#Install Anaconda2
RUN bash Miniconda2-latest-Linux-x86_64.sh -b -p /run_xgam/anaconda2
ENV PATH /run_xgam/anaconda2/bin:$PATH
SHELL ["/bin/bash", "-c"]

#Install healpy and fermitools
RUN conda update -n base -c defaults conda
RUN conda config --add channels conda-forge
RUN conda create -n fermi -c conda-forge/label/cf201901 -c fermi fermitools
RUN conda install -y --name fermi healpy

#Clone Xgam
WORKDIR /run_xgam
RUN git clone https://github.com/nmik/Xgam.git
COPY config/config_dataselection.py /run_xgam/Xgam/config/config_dataselection.py
#RUN less /run_xgam/Xgam/config/config_dataselection.py

#Creating bashrc file
RUN echo "echo 'Setting Xgam environment...'" > /run_xgam/.bashrc \
#RUN echo "export PATH=/run_xgam/anaconda2/bin:$PATH" >> /run_xgam/.bashrc \
   && echo "export PATH=/run_xgam/anaconda2/envs/fermi/bin:$PATH" >> /run_xgam/.bashrc \
   && echo "export PYTHONPATH=:/run_xgam/:${PYTHONPATH}" >> /run_xgam/.bashrc \
   && echo "export PATH=/run_xgam/Xgam/bin:${PATH}" >> /run_xgam/.bashrc \
   && echo "export P8_DATA=/archive/home/sammazza/fermi_data" >> /run_xgam/.bashrc \
   && echo "export X_OUT=/archive/home/sammazza/fermi_data" >> /run_xgam/.bashrc \
   && echo "export X_OUT_FIG=/archive/home/sammazza/fermi_data" >> /run_xgam/.bashrc \
   && echo "source activate fermi" >> /run_xgam/.bashrc \
   && echo "echo 'Done.'" >> /run_xgam/.bashrc \
RUN echo "bashrc file:" && less /run_xgam/.bashrc

WORKDIR /archive/home/sammazza/fermi_data
## Define entrypoint and default values for args
CMD ["/bin/bash","-c","source /run_xgam/.bashrc && /archive/home/sammazza/fermi_data/bash_script.sh"]
