drop view orphaned_instance;
drop view orphaned_file;
drop function "api"."user_config"(json);
drop function "api"."add_file"(varchar, varchar, jsonb);
drop function "api"."add_instance"(varchar, jsonb);
drop function "ensure_user"();
drop table "user" cascade;
drop role "standard";
drop role "admin";