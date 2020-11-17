# Publishing Appyters

To publish your Appyter, please submit pull request on GitHub for the Appyters Catalog GitHub repository (<https://github.com/MaayanLab/appyter-catalog>) following these instructions: 

1. Submit pull request with new appyter added to appyters directory.
2. `.github/workflows/validate_merge.yml` instructs github to execute `validate/validate_merge.py`
3. `validate/validate_merge.py` executes, validating the structure of the directory including
    1. Asserting that `appyter.json` is formatted according to the `schema/appyter-validator.json` json-schema validator
    2. Asserting that other relevant files are present
    3. Uses `compose/build_dockerfile.py` to construct and build a Dockerfile the same way it would be done in production
4. PR is accepted if and only if the validation and manual review is passed
5. Makefile can be used to facilitate the remaining steps
6. Run `compose/build_dockerfile.py` for each appyter to inject overrides, `merge_j2`, and construct a Dockerfile for the appyter
    1. When built, the files in override will be merged (using `compose/merge_j2.py`) with the appyter's own appyter overrides
7. Run `compose/build_appyters.py` to build a unified `appyters.json` file, containing information about each appyter for the app
8. Run `cd app && npm i && npm run build` to build the app (written in nodejs) with the most recently rendered `appyters.json`
9. Run `compose/build_compose.py` to build a application wide `docker-compose.yml` which includes a unified proxy for serving all apps on one endpoint
10. Run docker-compose build to build all `Dockerfile`s for the appyters and the app
    1. Variables in `.env` are automatically loaded by docker-compose
    2. `appyter_version` in `.env` is used as a `Dockerfile` arg permitting easy updates to the version used by all appyters
    3. A postgres database is used through postgrest for app state.
        1. `postgres/migrations` contains the postgres schema of that database, which are applied at database initialization in `postgres/Dockerfile`
11. Run `docker-compose up -d` to start all docker containers in the application.
    1. Variables in `.env` are automatically loaded by docker-compose
    2. `maayanlab/proxy` is used to proxy different paths to the respective containers and set up https with letsencrypt
    3. postgrest exposes postgres tables, views, and functions on the api schema with the guest role over HTTP at `/postgrest`
    4. app facilitates showing all appyters and navigating users to the mount location of the actual appyter container.
    5. `data/<container_name>` contains all application data split up by container mounted from the host.
