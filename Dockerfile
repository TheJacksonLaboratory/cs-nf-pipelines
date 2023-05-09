FROM ubuntu:kinetic
 
RUN apt-get update \
    && apt-get install -y eatmydata \
    && eatmydata apt-get install -y build-essential wget bzip2 \
      ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 libz-dev \
      git unzip python3.10 pip \
    && apt-get clean
 
# make sure we can just issue "python"
RUN ln -s /usr/bin/python3.10 /usr/bin/python
 
# bowtie
RUN wget -q -O bowtie.zip https://github.com/BenLangmead/bowtie/releases/download/v1.3.1/bowtie-1.3.1-linux-x86_64.zip; \
      unzip bowtie.zip -d /opt/; \
      ln -s /opt/bowtie-1.3.1-linux-x86_64 /opt/bowtie; \
      rm bowtie.zip
ENV PATH $PATH:/opt/bowtie
 
# bowtie2
RUN wget -q -O bowtie2.zip https://github.com/BenLangmead/bowtie2/releases/download/v2.5.1/bowtie2-2.5.1-linux-x86_64.zip; \
      unzip bowtie2.zip -d /opt/; \
      ln -s /opt/bowtie2-2.5.1-linux-x86_64 /opt/bowtie2; \
      rm bowtie2.zip
ENV PATH $PATH:/opt/bowtie2
 
RUN cd / && git clone https://github.com/MikeWLloyd/emase-zero.git && cd emase-zero/src && make
RUN cd / && git clone https://github.com/churchill-lab/alntools --branch feature/py3 && cd alntools && pip install .
RUN cd / && git clone https://github.com/churchill-lab/emase.git --branch py3 && cd emase && pip install .
RUN cd / && git clone https://github.com/churchill-lab/gbrs.git --branch feature/py3 && cd gbrs && pip install .
