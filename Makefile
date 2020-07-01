PYTHON ?= python3

# Get all appyter directories
APPYTERS = $(shell ls -d appyters/*)
APPYTERFILES = $(foreach appyter, $(APPYTERS), $(appyter)/appyter.json)
DOCKERFILES = $(foreach appyter, $(APPYTERS), $(appyter)/Dockerfile)
BUILDAPPYTERS = $(foreach appyter, $(APPYTERS), $(appyter)/.build)

s+ = $(subst \ ,+,$1)
+s = $(subst +,\ ,$1)

.SECONDEXPANSION:
appyters/%/Dockerfile: compose/build_dockerfile.py compose/templates/Dockerfile.j2 compose/merge_j2.py $$(call +s,$$(shell find $$(@D) -type f ! \( -name Dockerfile -o -name .build \) | sed 's/ /+/g'))
	$(PYTHON) compose/build_dockerfile.py $(shell basename $(shell dirname $@)) > $@

docker-compose.yml: compose/build_compose.py compose/templates/docker-compose.yml.j2 app/Dockerfile $(DOCKERFILES)
	$(PYTHON) compose/build_compose.py > $@

app/public/appyters.json: $(APPYTERFILES)
	$(PYTHON) compose/build_appyters.py > $@

app/.build: app/public/appyters.json app/package.json $$(call +s,$$(shell find app/public -type f | sed 's/ /+/g'))
	cd app && npm i && npm run build && cd .. && docker-compose build app && touch $@

appyters/%/.build: appyters/%/Dockerfile
	docker-compose build $(shell basename $(shell dirname $@ | awk '{print tolower($$0)}')) && touch $@

.build: app/.build $(BUILDAPPYTERS)

.env: .env.example
	test -f .env || ( echo "Warning: Using .env.example, please update .env as required" && cp .env.example .env )

.PHONY: prepare
prepare: .env docker-compose.yml $(DOCKERFILES)

.PHONY: build
build: .build

.PHONY: push
push: build
	docker-compose push

.PHONY: pull
pull: prepare
	docker-compose pull

.PHONY: start
start: prepare
	( make pull || make build ) && docker-compose up -d

.PHONY: stop
stop: prepare
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
