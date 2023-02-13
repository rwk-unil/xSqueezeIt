FROM ubuntu:20.04

# This is a Dockerfile for running xSqueezeIt
MAINTAINER Rick Wertenbroek <rick.wertenbroek@unil.ch>

# Install required software and clean as not to make the layer dirty
RUN apt-get update && apt-get -y upgrade && apt-get install -y \
	apt-utils curl gnupg gcc g++ make autoconf && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN apt-get update && apt-get -y upgrade && apt-get install -y \
	git zlib1g-dev libbz2-dev liblzma-dev && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Clone Source Code Repository and build
RUN mkdir -p /usr/src/ && \
    cd /usr/src/ && \
    git clone https://github.com/rwk-unil/xSqueezeIt.git && \
    cd /usr/src/xSqueezeIt && \
    git submodule update --init --recursive htslib && \
    cd htslib && \
    autoheader && \
    autoconf && \
    ./configure && \
    make && \
    make install && \
    ldconfig && \
    cd .. && \
    git clone https://github.com/facebook/zstd.git && \
    cd zstd && \
    make && \
    cd .. && \
    make && \
    chmod +x xsqueezeit && \
    cp xsqueezeit /usr/local/bin/ && \
    cd /usr/src/ && \
    rm -r xSqueezeIt

# Work in this temporary directory
WORKDIR /tmp/work

CMD echo "Run with the following command : docker run <tag> xsqueezeit [args]"