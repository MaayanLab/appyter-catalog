FROM postgres
ADD ./migrations/01_init/up.sql /docker-entrypoint-initdb.d/01_init.sql
ADD ./migrations/02_pagehits/up.sql /docker-entrypoint-initdb.d/02_pagehits.sql