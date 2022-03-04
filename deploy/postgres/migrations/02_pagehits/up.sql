-- add a table to store pagehits
create table api.pagehits (
  url varchar primary key not null,
  hits bigint not null default 0
);
grant all privileges on api.pagehits to appyters;
grant select on api.pagehits to guest;

-- add public pagehit function for incrementing pagehits
--  create page with 1 hit or or increment page hits
create or replace function api.pagehit(pageurl varchar) returns void as
$$
begin
  insert into api.pagehits (url, hits)
  values (pageurl, 1)
  on conflict (url)
  do update
  set hits = api.pagehits.hits + 1;
end
$$
language 'plpgsql' security definer;
grant all privileges on function api.pagehit(varchar) to appyters;
grant execute on function api.pagehit(varchar) to guest;
