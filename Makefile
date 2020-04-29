# These are recipes for getting all template docker-files setup in a unified docker-compose

# Get all template directories
TEMPLATES=$(shell ls -d templates/*)

# The Dockerfile in each template directory
DOCKERFILES=$(foreach template, $(TEMPLATES), $(template)/Dockerfile)

# You can make a dockerfile for a given template with compose/build_dockerfile.py
$(DOCKERFILES):
	python3 compose/build_dockerfile.py $(shell basename $(shell dirname $@)) > $@

# You can make a docker-compose given that all dockerfiles are present with compose/build_compose.py
docker-compose.yml: $(DOCKERFILES)
	python3 compose/build_compose.py > $@


# These are just convenient wrappers for docker-compose depending on the fact that the docker-compose exists

.PHONY: build
build: docker-compose.yml
	docker-compose build

.PHONY: start
start: docker-compose.yml
	docker-compose up -d

.PHONY: stop
stop: docker-compose.yml
	docker-compose down


# These are convenient wrappers for deployment management

# clean this directory of any un-tracked files
.PHONY: clean
clean:
	git reset --hard

# upgrade the project to the latest git commit
.PHONY: upgrade
upgrade: clean
	git pull --rebase

# update the deployment, upgrading the source and restarting
.PHONY: update
update:
	make upgrade && make build && make start
