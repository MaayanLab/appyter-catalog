# appyter-catalog: A catalog of [appyters](https://github.com/maayanLab/appyter/)

Pull requests encouraged, please refer to the [example](./appyters/example/) for registering your own appyter.

## Deployment

Currently, because this application deals with several independent appyters, we construct Dockerfiles independently for each and facilitate deployment with docker-compose.yml. In the future this can be extended to automatically generating a kubernetes deployment or simply using docker-in-docker, but for now a simple Makefile will do the trick of hosting the docker-compose on a single system.

```bash
# Download the catalog locally
git clone git@github.com:MaayanLab/appyter-catalog.git

# Build+run everything
make deploy

# Build+run specific components
make appyters/example/.deploy app/.deploy

# Publish specific component
make appyters/example/.publish
```

## Details

The appter-catalog does several things to permit integration of several independent appyters with their own dependencies while permitting various modifications performed at the entire application level.

1. Submit pull request with new appyter added to `appyters` directory.
2. `.github/workflows/validate_merge.yml` instructs github to execute `validate/validate_merge.py`
3. `validate/validate_merge.py` executes, validating the structure of the directory including
    1. Asserting that `appyter.json` is formatted according to the `schema/appyter-validator.json` json-schema validator
    2. Asserting that other relevant files are present
    3. Uses `compose/build_dockerfile.py` to construct and build a Dockerfile the same way it would be done in production
4. PR is accepted if and only if the validation and manual review is passed
5. `Makefile` can be used to facilitate the remaining steps
6. Run `compose/build_dockerfile.py` for each appyter to inject `override`s, `catalog_helper`, and construct a Dockerfile for the `appyter`
    1. When built, the files in `deploy/override` will be merged (using `compose/catalog_helper.py`) with the appyter's own `appyter` overrides
7. Run `compose/build_appyters.py` to build a unified `appyters.json` file, containing information about each appyter for the `app`
8. Run `cd app && npm i && npm run build` to build the `app` (written in nodejs) with the most recently rendered `appyters.json`
9. Run `compose/build_compose.py` to build a application wide `docker-compose.yml` which includes a unified proxy for serving all apps on one endpoint
10. Run `docker-compose build` to build all Dockerfiles for the `appyters` and the `app`
    1. Variables in `.env` are automatically loaded by `docker-compose`
    2. `appyter_version` in `.env` is used as a `Dockerfile arg` permitting easy updates to the version used by all `appyters`
    3. A `postgres` database is used through `postgrest` for `app` state.
        1. `postgres/migrations` contains the `postgres` schema of that database, which are applied at database initialization in `postgres/Dockerfile`
11. Run `docker-compose up -d` to start all docker containers in the application.
    1. Variables in `.env` are automatically loaded by `docker-compose`
    2. `maayanlab/proxy` is used to proxy different paths to the respective containers and set up `https` with `letsencrypt`
    3. `postgrest` exposes `postgres` tables, views, and functions on the `api` schema with the `guest` role over HTTP at `/postgrest`
    4. `app` facilitates showing all `appyters` and navigating users to the mount location of the actual `appyter` container.
    5. `data/<container_name>` contains all application data split up by container mounted from the host.

