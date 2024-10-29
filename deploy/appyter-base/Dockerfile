FROM ubuntu AS base

ENV DEBIAN_FRONTEND="noninteractive"
ENV TZ="America/New_York"

RUN set -x \
  && apt-get update -y \
  && apt-get install -y \
    build-essential \
    curl \
    libbz2-dev \
    libffi-dev \
    liblzma-dev \
    libreadline-dev \
    libsqlite3-dev \
    libssl-dev \
    zlib1g-dev \
    fuse \
    git \
    nginx \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*

ENV CARGO_HOME="/opt/rust/cargo"
ENV RUSTUP_HOME="/opt/rust/rustup"
RUN set -x \
  && curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | bash -s -- -y
ENV PATH="/opt/rust/cargo/bin:$PATH"

RUN set -x \
  && mkdir -p /opt/rig \
  && curl -Ls https://github.com/r-lib/rig/releases/download/v0.7.0/rig-linux-0.7.0.tar.gz | tar xz -C /opt/rig
ENV PATH="/opt/rig/bin:$PATH"

RUN set -x \
  && mkdir -p /opt/pyenv \
  && curl -Ls https://github.com/pyenv/pyenv/archive/refs/tags/v2.4.17.tar.gz | tar xz -C /opt/pyenv --strip-components=1

ENV PYENV_ROOT=/opt/pyenv
ENV PATH="/opt/pyenv/shims:/opt/pyenv/bin:$PATH"
RUN set -x \
  && echo 'eval "$(pyenv init -)"' >> /root/.bashrc \
  && echo 'eval "$(pyenv init -)"' >> /root/.profile

FROM base AS py
RUN set -x \
  && echo "Installing python..." \
  && pyenv install 3.8.20 \
  && pyenv global 3.8.20 \
  && pyenv rehash \
  && pip install -U pip

FROM py as rpy
RUN set -x \
  && echo "Installing R..." \
  && rig add 4.1.3 \
  && rig default 4.1.3 \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/*
