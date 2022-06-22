FROM ubuntu:20.04

MAINTAINER Akiyama Laboratory, Tokyo Institute of Technology <megadock@bi.cs.titech.ac.jp>

RUN apt-get update -y
RUN apt-get install -y --no-install-recommends \
                wget \
                make \
                g++ \
                openmpi-bin \
                libopenmpi-dev \
                ssh \
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
COPY Dockerfiles/mpi/Makefile ${MEGADOCK_WORK_DIR}/Makefile

RUN cd ${MEGADOCK_WORK_DIR} && \
    make -j$(nproc)

RUN ln -s ${MEGADOCK_WORK_DIR}/megadock-dp /usr/local/bin/megadock-dp && \
    ln -s ${MEGADOCK_WORK_DIR}/decoygen /usr/local/bin/decoygen && \
    ln -s ${MEGADOCK_WORK_DIR}/block /usr/local/bin/block && \
    ln -s ${MEGADOCK_WORK_DIR}/ppiscore /usr/local/bin/ppiscore

WORKDIR ${MEGADOCK_WORK_DIR}

# CONFIGURATION
ENV OMP_NUM_THREADS 8
RUN echo "export OMP_NUM_THREADS=${OMP_NUM_THREADS}" >> /etc/profile
