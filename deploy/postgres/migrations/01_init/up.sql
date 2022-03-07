-- setup anonymous user (guest)
create role guest nologin;
grant guest to appyters;

-- setup api schema
create schema if not exists api;
grant all privileges on schema api to appyters;
grant usage on schema api to guest;