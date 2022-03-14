VERSION=v1.7.1

crds/cert-manager.crds.yaml:
	curl -L https://github.com/cert-manager/cert-manager/releases/download/$(VERSION)/cert-manager.crds.yaml -o $@

templates/cert-manager.yaml:
	helm template cert-manager jetstack/cert-manager --version $(VERSION) > $@

.PHONY: build
build: crds/cert-manager.crds.yaml templates/cert-manager.yaml

.PHONY: clean
clean:
	rm crds/cert-manager.crds.yaml templates/cert-manager.yaml
