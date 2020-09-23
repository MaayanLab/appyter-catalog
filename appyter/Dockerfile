FROM ubuntu

ENV DEBIAN_FRONTEND "noninteractive"
ENV TZ "America/New_York"

RUN set -x \
  && echo "Preparing system..." \
  && apt-get -y update \
  && apt-get -y install \
    curl \
    fuse \
    git \
    nginx \
    python3-dev \
    python3-pip \
    rclone \
  && curl https://get.docker.com | sh \
  && rm -rf /var/lib/apt/lists/* \
  && pip3 install --upgrade pip \
  && chmod og+rwx -R /var/lib/nginx /var/log/nginx

RUN set -x \
  && echo "Preparing user..." \
  && useradd -ms /bin/bash -d /app app \
  && adduser app docker

RUN set -x \
  && echo "Installing jupyter kernel..." \
  && pip3 install ipykernel \
  && python3 -m ipykernel install

ARG appyter_version=git+git://github.com/Maayanlab/appyter.git
RUN set -x \
  && echo "Installing appyter..." \
  && pip3 install --upgrade ${appyter_version}

WORKDIR /app

CMD [ "appyter" ]
