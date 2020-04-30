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
    gcc

WORKDIR /usr/src/app

COPY setup.R ./
RUN R -e "source('setup.R')"

COPY requirements.txt ./
RUN pip3 install --no-cache-dir -r requirements.txt

COPY . .
COPY ./biojupies_bulk2/ /usr/local/lib/python3.8/dist-packages/jupyter_template/profiles/biojupies_bulk/
COPY run.sh ./
RUN set +x run.sh

EXPOSE 81
CMD [ "sh", "./run.sh" ]

