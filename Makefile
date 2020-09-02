PYTHON ?= python3

# Get all appyter directories
APPYTERS = $(shell find appyters -name appyter.json -exec sh -c 'realpath --relative-to=appyters $$(dirname {})' \;)
APPYTERFILES = $(foreach appyter, $(APPYTERS), appyters/$(appyter)/appyter.json)
DOCKERFILES = $(foreach appyter, $(APPYTERS), appyters/$(appyter)/Dockerfile)
BUILDAPPYTERS = $(foreach appyter, $(APPYTERS), appyters/$(appyter)/.build)

s+ = $(subst \ ,+,$1)
+s = $(subst +,\ ,$1)

.SECONDEXPANSION:
compose/.build: $$(shell find $$(@D) -type f ! \( -name Dockerfile -o -name .build \) | sed 's/ /+/g')
	touch $@

.SECONDEXPANSION:
$(DOCKERFILES): compose/.build $$(call +s,$$(shell find $$(@D) -type f ! \( -name Dockerfile -o -name .build \) | sed 's/ /+/g'))
	$(PYTHON) compose/build_dockerfile.py $(shell basename $(shell dirname $@)) > $@

.SECONDEXPANSION:
$(BUILDAPPYTERS): docker-compose.yml $$(@D)/Dockerfile
	docker-compose build appyter-$(shell basename $(shell dirname $@ | awk '{print tolower($$0)}')) && touch $@

docker-compose.yml: compose/.build .env $(DOCKERFILES)
	$(PYTHON) compose/build_compose.py > $@

app/public/appyters.json: $(APPYTERFILES)
	$(PYTHON) compose/build_appyters.py > $@

app/.build: app/public/appyters.json app/package.json $$(call +s,$$(shell find app/public -type f | sed 's/ /+/g'))
	cd app && npm i && npm run build && cd .. && docker-compose build appyters-catalog && touch $@

.build: app/.build $(BUILDAPPYTERS)

.env: .env.example
	test -f .env || ( echo "Warning: Using .env.example, please update .env as required" && cp .env.example .env )

.PHONY: build
build: .build
