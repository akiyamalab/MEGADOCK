FROM nvidia/cuda:11.1.1-devel-ubuntu20.04

MAINTAINER Akiyama Laboratory, Tokyo Institute of Technology <megadock@bi.cs.titech.ac.jp>

RUN apt-get update -y
RUN apt-get install -y --no-install-recommends \
                wget \
                make \
                g++ \
                cuda-samples-11-1 \
    && rm -rf /var/lib/apt/lists/*

# FFTW
ENV FFTW_TARGET fftw-3.3.8

RUN wget -P /tmp http://www.fftw.org/${FFTW_TARGET}.tar.gz && \
    tar xzf /tmp/${FFTW_TARGET}.tar.gz -C /tmp && \
    rm -rf /tmp/${FFTW_TARGET}.tar.gz && \
    cd /tmp/${FFTW_TARGET} && \
    ./configure --enable-float --enable-sse2 --prefix=/usr/local/${FFTW_TARGET} && \
    make -j$(nproc) && \
    make install


# MEGADOCK
ENV MEGADOCK_WORK_DIR /opt/MEGADOCK

ADD . ${MEGADOCK_WORK_DIR}
COPY Dockerfiles/gpu/Makefile ${MEGADOCK_WORK_DIR}/Makefile

RUN cd ${MEGADOCK_WORK_DIR} && \
    make -j$(nproc)

WORKDIR ${MEGADOCK_WORK_DIR}

RUN ln -s ${MEGADOCK_WORK_DIR}/megadock-gpu /usr/local/bin/megadock-gpu && \
    ln -s ${MEGADOCK_WORK_DIR}/decoygen /usr/local/bin/decoygen && \
    ln -s ${MEGADOCK_WORK_DIR}/block /usr/local/bin/block && \
    ln -s ${MEGADOCK_WORK_DIR}/ppiscore /usr/local/bin/ppiscore

