PYTHON ?= python3
COMPOSE_ARGS ?= 

# Get all appyter directories
APPYTERS = $(shell find appyters -name appyter.json -exec sh -c 'realpath --relative-to=appyters $$(dirname {})' \;)
APPYTER_FILES = $(foreach appyter, $(APPYTERS), appyters/$(appyter)/appyter.json)
DOCKERFILES = $(foreach appyter, $(APPYTERS), appyters/$(appyter)/Dockerfile)
BUILD_APPYTERS = $(foreach appyter, $(APPYTERS), appyters/$(appyter)/.build)
PUBLISH_APPYTERS = $(foreach appyter, $(APPYTERS), appyters/$(appyter)/.publish)
DEPLOY_APPYTERS = $(foreach appyter, $(APPYTERS), appyters/$(appyter)/.deploy)

s+ = $(subst \ ,+,$1)
+s = $(subst +,\ ,$1)

.SECONDEXPANSION:
compose/.build: $$(shell find $$(@D) -type f ! \( -name Dockerfile -o -name .build \) | sed 's/ /+/g')
	touch $@

.SECONDEXPANSION:
$(DOCKERFILES): compose/.build $$(call +s,$$(shell find $$(@D) -type f ! \( -name Dockerfile -o -name .build -o -name .deploy -o -name .publish \) | sed 's/ /+/g'))
	$(PYTHON) compose/build_dockerfile.py $(shell basename $(shell dirname $@)) > $@

.SECONDEXPANSION:
$(BUILD_APPYTERS): docker-compose.yml $$(@D)/Dockerfile
	docker-compose build appyter-$(shell basename $(shell dirname $@ | awk '{print tolower($$0)}')) && touch $@

.SECONDEXPANSION:
$(PUBLISH_APPYTERS): docker-compose.yml $$(@D)/.build
	docker-compose push appyter-$(shell basename $(shell dirname $@ | awk '{print tolower($$0)}')) && touch $@

.SECONDEXPANSION:
$(DEPLOY_APPYTERS): docker-compose.yml .env $$(@D)/.build
	docker-compose up -d appyter-$(shell basename $(shell dirname $@ | awk '{print tolower($$0)}')) && touch $@

docker-compose.yml: compose/.build $(DOCKERFILES)
	$(PYTHON) compose/build_compose.py $(COMPOSE_ARGS) > $@

app/public/appyters.json: compose/.build .env $(APPYTER_FILES)
	$(PYTHON) compose/build_appyters.py > $@

app/.build: app/public/appyters.json app/package.json $$(call +s,$$(shell find app/public -type f | sed 's/ /+/g'))
	cd app && npm i && npm run build && cd .. && docker-compose build appyters-catalog && touch $@

app/.publish: app/.build
	docker-compose push appyters-catalog && touch $@

app/.deploy: app/.build .env
	docker-compose up -d appyters-catalog && touch $@

.build: app/.build $(BUILD_APPYTERS)
.publish: app/.publish $(PUBLISH_APPYTERS)
.deploy: app/.deploy $(DEPLOY_APPYTERS)

.env: .env.example
	test -f .env || ( echo "Warning: Using .env.example, please update .env as required" && cp .env.example .env )

.PHONY: build
build: .build

.PHONY: publish
publish: .publish

.PHONY: deploy
deploy: .deploy
