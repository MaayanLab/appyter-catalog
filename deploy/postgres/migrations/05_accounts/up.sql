create role "standard";
grant usage on schema "api" to "standard";

create role "admin";
grant usage on schema "api" to "admin";

-- capture the user in this database
create table "user" (
  "id" varchar default current_setting('request.jwt.claim.email', false),
  "metadata" json,
  "ts" timestamp not null default now(),
  primary key ("id")
);
alter table "user" enable row level security;
create policy "user_policy" on "user"
  for all
  using ("id" = current_setting('request.jwt.claim.email', false))
  with check ("id" = current_setting('request.jwt.claim.email', false));
grant all on "user" to "standard";

-- ensure the user is registered in our database
create or replace function "ensure_user"() returns boolean as
$$
declare
  email varchar;
begin
  select nullif(current_setting('request.jwt.claim.email', true), '')
  into email;

  if email is not null then
    insert into "public"."user" ("id", "metadata", "ts")
    values (
      email,
      json_build_object(
        'name', current_setting('request.jwt.claim.name', true)
      ),
      now()
    )
    on conflict ("id")
    do nothing;
    return true;
  else
    return false;
  end if;
end
$$
language 'plpgsql';
grant execute on function "ensure_user"() to "guest";
grant execute on function "ensure_user"() to "standard";

-- save user configuration
create table "user_config" (
  "user" varchar default current_setting('request.jwt.claim.email', false),
  "config" json not null,
  primary key ("user"),
  foreign key ("user") references "user" ("id")
);
alter table "user_config" enable row level security;
create policy "user_config_policy" on "user_config"
  for all to "standard"
  using ("user" = current_setting('request.jwt.claim.email', true))
  with check ("user" = current_setting('request.jwt.claim.email', true));
grant all on "user_config" to "standard";

-- end point for fetching/writing user config
create or replace function "api"."user_config"(config json) returns json as
$$
declare
  _config json;
begin
  if (select "ensure_user"() is false) then
    raise exception 'Must be logged in to use user_config';
  end if;
  if config is null then
    select uc."config"
    from "user_config" uc
    where uc."user" = current_setting('request.jwt.claim.email', true)
    into _config;

    if _config is null then
      insert into "public"."user_config" ("config")
      values ('{}'::json);
      _config := '{}'::json;
    end if;
  else
    insert into "public"."user_config" ("user", "config")
    values (current_setting('request.jwt.claim.email', true), "config")
    on conflict ("user")
    do update set "config" = EXCLUDED."config"
    returning "public"."user_config"."config"
    into _config;
  end if;
  return _config;
end
$$
language 'plpgsql';
grant execute on function "api"."user_config"(json) to "standard";

create table "user_file" (
  "id" uuid default uuid_generate_v4(),
  "user" varchar default current_setting('request.jwt.claim.email', false),
  "file" varchar,
  "filename" varchar,
  "ts" timestamp default now(),
  "metadata" jsonb,
  primary key ("id"),
  unique ("user", "file", "filename"),
  foreign key ("user") references "user" ("id"),
  foreign key ("file") references "file" ("id")
);
alter table "user_file" enable row level security;
create policy "user_file_policy" on "user_file" for all
  using ("user" = current_setting('request.jwt.claim.email', true))
  with check ("user" = current_setting('request.jwt.claim.email', true));
grant all on "user_file" to "standard";
grant select on "file" to "standard";

create or replace function "api"."add_file"("file" varchar, "filename" varchar, "metadata" jsonb)
returns void as
$$
#variable_conflict use_column
begin
  insert into "public"."file" ("id", "metadata")
  values ($1, $3)
  on conflict ("id")
  do nothing;

  if (select "ensure_user"()) then
    insert into "public"."user_file" ("file", "filename")
    values ($1, $2)
    on conflict ("user", "file", "filename")
    do nothing;
  end if;
end
$$
language 'plpgsql'
security definer;
grant execute on function "api"."add_file"(varchar, varchar, jsonb) to "guest";
grant execute on function "api"."add_file"(varchar, varchar, jsonb) to "standard";

create or replace view "api"."user_file" as
select
  uf.id,
  row_to_json(f.*) as file,
  uf.filename,
  uf.ts,
  uf.metadata
from "user_file" uf
inner join "file" f on f."id" = uf."file"
;
alter view "api"."user_file" owner to "standard";

create or replace function api_user_file_delete()
returns trigger as
$$
begin
  delete from "public"."user_file" uf
  where uf.id = old."id";
  return new;
end
$$
language 'plpgsql';
alter function api_user_file_delete() owner to "standard";

create trigger api_user_file_delete_trigger
instead of delete
on "api"."user_file"
for each row
execute procedure api_user_file_delete();

create table "user_instance" (
  "id" uuid default uuid_generate_v4(),
  "user" varchar default current_setting('request.jwt.claim.email', false),
  "instance" varchar,
  "ts" timestamp default now(),
  "metadata" jsonb,
  primary key ("id"),
  unique ("user", "instance"),
  foreign key ("user") references "user" ("id"),
  foreign key ("instance") references "instance" ("id")
);
alter table "user_instance" enable row level security;
create policy "user_instance_policy" on "user_instance" for all
  using ("user" = current_setting('request.jwt.claim.email', false))
  with check ("user" = current_setting('request.jwt.claim.email', false));
grant all on "user_instance" to "standard";
grant select on "instance" to "standard";

create or replace view "api"."user_instance" as
select
  ui.id,
  row_to_json(i.*) as instance,
  ui.ts,
  ui.metadata
from "user_instance" ui
inner join "instance" i on i."id" = ui."instance"
;
alter view "api"."user_instance" owner to "standard";

create or replace function api_user_instance_delete()
returns trigger as
$$
begin
  delete from "public"."user_instance" ui
  where ui.id = old."id"
  return new;
end
$$
language 'plpgsql';
alter function api_user_instance_delete() owner to "standard";

create trigger api_user_instance_delete_trigger
instead of delete
on "api"."user_instance"
for each row
execute procedure api_user_instance_delete();

create or replace function "api"."add_instance"("instance" varchar, "metadata" jsonb)
returns void as
$$
#variable_conflict use_column
begin
  insert into "public"."instance" ("id", "metadata")
  values ("instance", "metadata")
  on conflict ("id")
  do nothing;

  if (select "ensure_user"()) then
    insert into "public"."user_instance" ("instance")
    values ("add_instance"."instance")
    on conflict ("user", "instance")
    do nothing;
  end if;
end
$$
language 'plpgsql'
security definer;
grant execute on function "api"."add_instance"(varchar, jsonb) to "guest";
grant execute on function "api"."add_instance"(varchar, jsonb) to "standard";

create view "orphaned_file" as
select *
from "file" f
where not exists (
  select 1
  from "user_file" uf
  where uf."file" = f."id"
);

create view "orphaned_instance" as
select *
from "instance" i
where not exists (
  select 1
  from "user_instance" ui
  where ui."instance" = i."id"
);
