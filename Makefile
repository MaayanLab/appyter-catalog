PYTHON ?= python3

# Get all template directories
TEMPLATES = $(shell ls -d templates/*)
TEMPLATEFILES = $(foreach template, $(TEMPLATES), $(template)/template.json)
DOCKERFILES = $(foreach template, $(TEMPLATES), $(template)/Dockerfile)
BUILDTEMPLATES = $(foreach template, $(TEMPLATES), $(template)/.build)

.SECONDEXPANSION:
templates/%/Dockerfile: compose/build_dockerfile.py $$(shell find $$(@D) -type f ! \( -name Dockerfile -o -name .build \))
	$(PYTHON) compose/build_dockerfile.py $(shell basename $(shell dirname $@)) > $@

docker-compose.yml: app/Dockerfile $(DOCKERFILES)
	$(PYTHON) compose/build_compose.py > $@

app/public/templates.json: $(TEMPLATEFILES)
	cat $^ | jq -s '.' > $@

app/.build: app/public/templates.json app/package.json $$(shell find app/public -type f)
	cd app && npm i && npm run build && cd .. && docker-compose build app && touch $@

templates/%/.build: templates/%/Dockerfile
	docker-compose build $(shell basename $(shell dirname $@ | awk '{print tolower($$0)}')) && touch $@

.build: docker-compose.yml app/.build $(BUILDTEMPLATES)

.PHONY: build
build: .build

.env: .env.example
	test -f .env || ( echo "Warning: Using .env.example, please update .env as required" && cp .env.example .env )

.PHONY: start
start: build .env
	docker-compose up -d

.PHONY: stop
stop: docker-compose.yml
	docker-compose down

.PHONY: clean
clean:
	git reset --hard

.PHONY: upgrade
upgrade: clean
	git pull --rebase

.PHONY: update
update:
	make upgrade && make start
