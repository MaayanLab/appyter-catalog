# jupyter-template-catalog

A catalog of jupyter templates.

Pull requests encouraged, please refer to the [example](./templates/example/) for registering your own template.


## Deployment

Currently, because this application deals with several independent templates, we construct Dockerfiles independently for each and facilitate deployment with docker-compose.yml. In the future this can be extended to automatically generating a kubernetes deployment or simply using docker-in-docker, but for now a simple Makefile will do the trick of hosting the docker-compose on a single system.

```bash
# Download the catalog locally
git clone git@github.com:MaayanLab/jupyter-template-catalog.git

# Start the server (this has several dependent steps including constructing dockerfiles, docker-compose and building it all)
make start

# Update the server with latest github (WARNING: this will force delete anything not tracked in your current directory, then restart the application)
make update
```