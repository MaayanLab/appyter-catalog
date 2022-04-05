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
create policy "user_policy" on "user"
  using ("id" = current_setting('request.jwt.claim.email', true))
  with check ("id" = current_setting('request.jwt.claim.email', true));
grant all on "user" to "standard";

-- ensure the user is registered in our database
create or replace function "ensure_user"() returns void as
$$
begin
  insert into "public"."user" ("id", "metadata", "ts")
  values (
    current_setting('request.jwt.claim.email', true),
    json_build_object(
      'name', current_setting('request.jwt.claim.name', true)
    ),
    now()
  )
  on conflict ("id")
  do nothing;
end
$$
language 'plpgsql';
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
  using ("user" = current_setting('request.jwt.claim.email', true))
  with check ("user" = current_setting('request.jwt.claim.email', true));
grant all on "user_config" to "standard";

-- end point for fetching/writing user config
create or replace function "api"."user_config"(config json) returns json as
$$
declare
  _config json;
begin
  perform "ensure_user"();
  if config is null then
    select uc."config"
    from "user_config" uc
    where uc."user" = current_setting('request.jwt.claim.email', true)
    into _config;
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
  unique ("user", "instance", "filename"),
  foreign key ("user") references "user" ("id"),
  foreign key ("file") references "file" ("id")
);
alter table "user_file" enable row level security;
create policy "user_file_policy" on "user_file"
  using ("user" = current_setting('request.jwt.claim.email', true))
  with check ("user" = current_setting('request.jwt.claim.email', true));
grant all on "user_file" to "standard";

create view "api"."user_file" as
select * from "user_file";
alter view "api"."user_file" owner to "standard";

create table "user_instance" (
  "id" uuid default uuid_generate_v4(),
  "user" varchar default current_setting('request.jwt.claim.email', true),
  "instance" varchar,
  "ts" timestamp default now(),
  "metadata" jsonb,
  primary key ("id"),
  unique ("user", "instance"),
  foreign key ("user") references "user" ("id"),
  foreign key ("instance") references "instance" ("id")
);
alter table "user_instance" enable row level security;
create policy "user_instance_policy" on "user_instance"
  using ("user" = current_setting('request.jwt.claim.email', true))
  with check ("user" = current_setting('request.jwt.claim.email', true));
grant all on "user_instance" to "standard";

create view "api"."user_instance" as
select
  ui.id,
  ui.user,
  row_to_json(i.*) as instance,
  ui.ts,
  ui.metadata
from "user_instance" ui
inner join "instance" i on i."id" = ui."instance"
;
grant all on "api"."user_instance" to "standard";

create or replace function api_user_instance_delete()
returns trigger as
$$
begin
  delete from "public"."user_instance"
  where "public"."user_instance".id = old."id";
  return new;
end
$$
language 'plpgsql';

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

  perform "ensure_user"();

  insert into "public"."user_instance" ("instance")
  values ("add_instance"."instance")
  on conflict ("user", "instance")
  do nothing;
end
$$
language 'plpgsql'
security definer;
grant execute on function "api"."add_instance"(varchar, jsonb) to "standard";
