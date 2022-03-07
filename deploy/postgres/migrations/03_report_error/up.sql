create extension if not exists "uuid-ossp";

-- add a table to store errors
create table error_report (
  id uuid default uuid_generate_v4(),
  appyter jsonb,
  url varchar,
  type varchar,
  error jsonb,
  ts timestamp default CURRENT_TIMESTAMP,
  addressed boolean default false,
  primary key (id)
);

-- add public report_error function for reporting an error
create or replace function api.report_error(
  appyter jsonb,
  url varchar,
  type varchar,
  error jsonb
) returns uuid as
$$
declare
  _id uuid;
begin
  insert into error_report (
    appyter,
    url,
    type,
    error
  ) values (
    appyter,
    url,
    type,
    error
  ) returning id
  into _id;
  return _id;
end
$$
language 'plpgsql' security definer;
grant execute on function api.report_error(jsonb, varchar, varchar, jsonb) to guest;
