# We will use Ubuntu for our image
FROM ubuntu:latest

ENV DEBIAN_FRONTEND "noninteractive"
ENV TZ "America/New_York"

RUN set -x \
  && apt-get update \
  && apt-get install -y\
    tzdata \
    python3 \
    python3-pip \
    git \
    r-base \
    make \
    gcc \
    build-essential \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev

WORKDIR /usr/src/app

COPY setup.R ./
RUN R -e "source('setup.R')"

COPY requirements.txt ./
RUN pip3 install --no-cache-dir -r requirements.txt

COPY . .
COPY run.sh ./
RUN set +x run.sh

ENV PREFIX "/"

EXPOSE 80
CMD [ "sh", "./run.sh" ]

