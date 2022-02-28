# Provision Kubernetes Cluster

## Locally
```bash
HOST_ARGS=$(getent hosts dev.appyter.maayanlab.cloud | python -c "
import sys
for argv in map(str.split, map(str.strip, sys.stdin)):
  print(f\"--host-alias {argv[0]}:{','.join(argv[1:])}\")
")
sh -c "k3d cluster create appyter -p 80:80@loadbalancer -p 443:443@loadbalancer ${HOST_ARGS}"
k3d kubeconfig merge appyter
kubectl get nodes
```

# Deploy appyter-catalog

## Install base cluster infrastructure
Install `appyter-catalog-base` cluster infrastructure which includes:
- secret-generator (auto-secrets)
- cert-manager (auto-tls)
- db-operator (auto-db provisioning)
- keycloak-operator (auth infrastructure)

```bash
helm dependency build charts/appyter-catalog-base
helm install appyter-catalog-base charts/appyter-catalog-base
```

## Install appyter-catalog
Utilizes the operators initialized in the base cluster infrastructure to auto-provision everything this application needs.
- create a postgres database and provision postgres-dbinstance

```bash
helm install appyter-catalog charts/appyter-catalog
```
