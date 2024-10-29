FROM maayanlab/appyter-base:py

RUN --mount=type=cache,target=/root/.cache/pip set -x \
  && echo "Installing jupyter kernel..." \
  && pip install ipython_genutils ipykernel \
  && python -m ipykernel install

ARG appyter_version=appyter[production]@git+https://github.com/Maayanlab/appyter
RUN --mount=type=cache,target=/root/.cache/pip set -x \
  && echo "Installing appyter..." \
  && pip install ${appyter_version}

RUN set -x \
  && echo "Preparing system..." \
  && apt-get -y update \
  && apt-get -y install rclone \
  && curl https://get.docker.com | sh \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* \
  && chmod og+rwx -R /var/lib/nginx /var/log/nginx

RUN set -x \
  && echo "Preparing user..." \
  && useradd -ms /bin/bash -d /app app \
  && adduser app docker \
  && groupadd fuse \
  && adduser app fuse

WORKDIR /app

CMD [ "appyter" ]
