FROM ubuntu:14.04
COPY Docker/coverage /bin/coverage
COPY Docker/v2.3.0.tar.gz /
RUN mkdir /covtobed
COPY *.cpp *.h /covtobed/
RUN apt-get update && apt-get install -y --no-install-recommends software-properties-common=0.92.37.8 \
     && apt-get clean \
     && rm -rf /var/lib/apt/lists/*

RUN add-apt-repository universe && apt-get update
RUN apt-get install -y --no-install-recommends build-essential=11.6ubuntu6 cmake=2.8.12.2-0ubuntu3 wget=1.15-1ubuntu1  zlib1g-dev=1:1.2.8.dfsg-1ubuntu1.1 \
     && apt-get clean \
     && rm -rf /var/lib/apt/lists/*
RUN tar xvfz "v2.3.0.tar.gz"
RUN mkdir /bamtools-2.3.0/build 
WORKDIR /bamtools-2.3.0/build
RUN cmake .. 
RUN make
RUN cp /bamtools-2.3.0/lib/libbamtools* /usr/lib/x86_64-linux-gnu/

WORKDIR /covtobed 
RUN c++ -std=c++11 ./*.cpp  -I/bamtools-2.3.0/include/  /usr/lib/x86_64-linux-gnu/libbamtools.a   \
	 -o /bin/covtobed -lz;
