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
