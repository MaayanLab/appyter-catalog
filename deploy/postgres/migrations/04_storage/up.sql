create table "file" (
  "id" varchar,
  "metadata" jsonb,
  "ts" timestamp default now(),
  primary key ("id")
);

create table "instance" (
  "id" varchar,
  "metadata" jsonb,
  "ts" timestamp default now(),
  primary key (id)
);