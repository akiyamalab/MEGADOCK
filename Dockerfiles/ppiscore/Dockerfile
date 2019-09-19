FROM python:3.6-alpine

RUN apk --no-cache --update-cache add g++ && \
    pip install numpy

ADD Dockerfiles/ppiscore/ppiscore.py /ppiscore.py

RUN chmod +x /ppiscore.py && ln -s /ppiscore.py /usr/local/bin/ppiscore