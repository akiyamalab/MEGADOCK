FROM ubuntu:16.04

MAINTAINER Akiyama Laboratory, Tokyo Institute of Technology <megadock@bi.cs.titech.ac.jp>

RUN apt-get update
RUN apt-get install -y --no-install-recommends \
                        wget \
                        make \
                        g++

# FFTW
ENV FFTW_TARGET fftw-3.3.5

RUN wget -P /tmp http://www.fftw.org/${FFTW_TARGET}.tar.gz  &&  \
    tar xzf /tmp/${FFTW_TARGET}.tar.gz -C /tmp  &&  \
    cd /tmp/${FFTW_TARGET}  &&  \
    ./configure --enable-float --enable-sse2 --prefix=/usr/local/${FFTW_TARGET}  &&  \
    make  &&  \
    make check  &&  \
    make install


# MEGADOCK
ENV MEGADOCK_WORK_DIR /opt/MEGADOCK

ADD . ${MEGADOCK_WORK_DIR}
COPY ./Dockerfiles/cpu/Makefile ${MEGADOCK_WORK_DIR}/Makefile

RUN cd ${MEGADOCK_WORK_DIR} && make

RUN ln -s ${MEGADOCK_WORK_DIR}/megadock /usr/local/bin/megadock \
    && ln -s ${MEGADOCK_WORK_DIR}/decoygen /usr/local/bin/decoygen  \
    && ln -s ${MEGADOCK_WORK_DIR}/block /usr/local/bin/block  \
    && ln -s ${MEGADOCK_WORK_DIR}/ppiscore /usr/local/bin/ppiscore

WORKDIR ${MEGADOCK_WORK_DIR}/data

# CONFIGURATION
ENV OMP_NUM_THREADS 4
RUN echo "export OMP_NUM_THREADS=${OMP_NUM_THREADS}" >> /etc/profile