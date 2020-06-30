FROM ubuntu:18.04
MAINTAINER Akshay Balsubramani abalsubr@stanford.edu 

#install miniconda 
ENV PATH /opt/conda/bin:$PATH
RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ENV PATH /opt/conda/bin:$PATH

RUN mkdir /var/www
WORKDIR /var/www

#install gcsfuse and copy data to the 
RUN curl -sSL https://sdk.cloud.google.com | bash
ENV PATH $PATH:/root/google-cloud-sdk/bin



#install coessentiality browser
RUN git clone https://github.com/kundajelab/coessentiality-browser.git
WORKDIR /var/www/coessentiality-browser
RUN gsutil -m cp -r gs://coessentiality-browser/data /var/www/coessentiality-browser/
RUN pip install -r requirements.txt
ENTRYPOINT [ "python" ]
CMD [ "app.py" ]

