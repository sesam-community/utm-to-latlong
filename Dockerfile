FROM ubuntu:16.04

MAINTAINER Tom Bech "tom.bech@sesam.io"

EXPOSE 5001/tcp

RUN \
apt-get update && \
DEBIAN_FRONTEND=noninteractive apt-get install -y locales && \
apt-get clean all && \
apt-get -y autoremove --purge && \
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN localedef -i en_US -f UTF-8 en_US.UTF-8
ENV LANGUAGE en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8
RUN locale-gen en_US.UTF-8
RUN DEBIAN_FRONTEND=noninteractive dpkg-reconfigure locales
ENV PYTHON_EGG_CACHE /tmp
ENV PYTHONIOENCODING UTF-8

RUN \
apt-get update && \
apt-get install -y \
  apt-utils \
  curl \
  python3 \
  python3-dev \
  build-essential

RUN \
  apt-get install -y \
  libffi-dev \
  libyajl-dev

RUN curl -sSL https://bootstrap.pypa.io/get-pip.py | python3
COPY ./service /service
WORKDIR /service
RUN pip install -r requirements.txt
CMD ["/usr/bin/python3", "transform-service.py"]
