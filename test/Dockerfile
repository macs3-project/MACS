ARG BASE_IMAGE=ubuntu:bionic
FROM ${BASE_IMAGE}

ENV TEST_USER macs

RUN uname -a
RUN apt-get update -qq && \
  apt-get install -yq \
  python3 \
  python3-pip \
  python3-setuptools \
  python3-wheel

RUN useradd -m "${TEST_USER}"
WORKDIR /build
COPY . .
RUN chown -R "${TEST_USER}:${TEST_USER}" /build

USER "${TEST_USER}"
RUN env | grep ^PATH
ENV PATH "/home/${TEST_USER}/.local/bin:${PATH}"

RUN python3 --version
RUN pip3 --version

RUN pip3 install --user --upgrade -r requirements.txt
RUN pip3 install --user .

RUN pytest --version
