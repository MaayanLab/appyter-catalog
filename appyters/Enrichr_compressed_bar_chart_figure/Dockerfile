FROM ubuntu

ENV DEBIAN_FRONTEND "noninteractive"
ENV TZ "America/New_York"

RUN set -x \
    && echo "Preparing system..." \
    && apt-get -y update \
    && apt-get -y install git python3-pip python3-dev \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install --upgrade pip

RUN set -x \
  && echo "Installing jupyter kernel..." \
  && pip3 install ipykernel \
  && python3 -m ipykernel install

ADD deps.txt /app/deps.txt
RUN set -x \
  && echo "Installing system dependencies from deps.txt..." \
  && apt-get -y update \
  && apt-get -y install $(grep -v '^#' /app/deps.txt) \
  && rm -rf /var/lib/apt/lists/* \
  && rm /app/deps.txt

ADD requirements.txt /app/requirements.txt
RUN set -x \
  && echo "Installing python dependencies from requirements.txt..." \
  && pip3 install -Ivr /app/requirements.txt \
  && rm /app/requirements.txt

ARG appyter_version=git+git://github.com/Maayanlab/appyter.git
RUN set -x \
  && echo "Installing appyter..." \
  && pip3 install -Iv ${appyter_version}

WORKDIR /app
EXPOSE 80
VOLUME /app/data

ENV PREFIX "/"
ENV HOST "0.0.0.0"
ENV PORT "80"
ENV DEBUG "false"

COPY . /app

CMD [ "appyter", "--profile=biojupies", "Enrichr_compressed_bar_chart_figure.ipynb" ]
