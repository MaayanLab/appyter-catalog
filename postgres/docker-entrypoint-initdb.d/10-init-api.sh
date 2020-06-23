#!/bin/bash

set -e

psql --username "${POSTGRES_USER}" <<EOSQL

create role ${POSTGRES_ANON_USER} nologin;
grant ${POSTGRES_ANON_USER} to ${POSTGRES_USER};

create schema ${POSTGRES_SCHEMA};
grant all privileges on schema ${POSTGRES_SCHEMA} to ${POSTGRES_USER};
grant usage on schema ${POSTGRES_SCHEMA} to ${POSTGRES_ANON_USER};

create table ${POSTGRES_SCHEMA}.hits (
  url varchar primary key,
  hits bigint
);
grant all privileges on ${POSTGRES_SCHEMA}.hits to ${POSTGRES_USER};
grant select on ${POSTGRES_SCHEMA}.hits to ${POSTGRES_ANON_USER};

EOSQL
