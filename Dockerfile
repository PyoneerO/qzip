FROM ubuntu:16.04

LABEL maintainer="stefan dot neuhaus at awi dot de"
LABEL description="https://github.com/PyoneerO/qzip/blob/master/README.MD"

#install needed packages from apt repos
RUN apt-get clean
RUN apt-get update
#RUN apt-get -y install gedit
#RUN apt-get -y install firefox
RUN apt-get -y install parallel
RUN apt-get -y install bc
RUN apt-get -y install python-dev
RUN apt-get -y install wget
RUN apt-get -y install openjdk-8-jdk
RUN apt-get -y install gcc
RUN apt-get -y install gawk
RUN apt-get -y install zip
RUN apt-get -y install nano
RUN apt-get -y install git
RUN apt-get autoremove 
RUN apt-get clean

#install python package index tool pip
RUN wget https://bootstrap.pypa.io/get-pip.py
RUN python get-pip.py
RUN rm get-pip.py

#install python packages via pip
RUN pip install cutadapt
RUN pip install numpy
RUN pip install biom-format
RUN pip install h5py

#get needed binaries and java archives
#vsearch
RUN wget -P /opt/ https://github.com/torognes/vsearch/releases/download/v2.8.2/vsearch-2.8.2-linux-x86_64.tar.gz
RUN tar xzf /opt/vsearch-2.8.2-linux-x86_64.tar.gz -C /opt/
RUN ln -s /opt/vsearch-2.8.2-linux-x86_64/bin/vsearch /opt/vsearch
#swarm
RUN wget -P /opt/ https://github.com/torognes/swarm/releases/download/v2.2.2/swarm-2.2.2-linux-x86_64
RUN chmod u+x /opt/swarm-2.2.2-linux-x86_64
RUN ln -s /opt/swarm-2.2.2-linux-x86_64 /opt/swarm
#mothur
RUN wget -P /opt/ https://github.com/mothur/mothur/releases/download/v1.40.5/Mothur.linux_64.noReadLine.zip
RUN unzip /opt/Mothur.linux_64.noReadLine.zip -d /opt/
RUN mv /opt/mothur/ /opt/mothur-1.40.5-linux-x86_64_noReadLine/
RUN ln -s /opt/mothur-1.40.5-linux-x86_64_noReadLine/mothur /opt/mothur
#trimmomatic
RUN wget -P /opt/ http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
RUN unzip /opt/Trimmomatic-0.38.zip -d /opt/
RUN chmod u+rx /opt/Trimmomatic-0.38/trimmomatic-0.38.jar
RUN ln -s /opt/Trimmomatic-0.38/trimmomatic-0.38.jar /opt/trimmomatic
#update PATH
RUN export PATH=/opt/:$PATH
#clean
RUN rm -r /opt/__MACOSX/ 
RUN rm /opt/*zip /opt/*tar.gz
