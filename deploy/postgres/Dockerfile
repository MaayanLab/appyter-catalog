FROM postgres
ADD ./migrations/01_init/up.sql /docker-entrypoint-initdb.d/01_init.sql
ADD ./migrations/02_pagehits/up.sql /docker-entrypoint-initdb.d/02_pagehits.sql
ADD ./migrations/03_report_error/up.sql /docker-entrypoint-initdb.d/03_report_error.sql
ADD ./migrations/04_storage/up.sql /docker-entrypoint-initdb.d/04_storage.sql
ADD ./migrations/05_accounts/up.sql /docker-entrypoint-initdb.d/05_accounts.sql
ADD ./migrations/06_keycloak/up.sql /docker-entrypoint-initdb.d/06_keycloak.sql