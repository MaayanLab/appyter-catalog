''' Produce acceptable.yaml
'''
import sys
import json
import yaml
import pandas as pd
import fsspec

acceptable_schemes = set()
acceptable_prefixes = set()

# IANA
iana_schemes = pd.read_csv('https://www.iana.org/assignments/uri-schemes/uri-schemes-1.csv', index_col=0)
acceptable_iana_schemes = iana_schemes[iana_schemes.Status.isin(['Permanent', 'Provisional'])]
acceptable_schemes.update(acceptable_iana_schemes.index)

# identifiers.org
with fsspec.open('https://registry.api.identifiers.org/resolutionApi/getResolverDataset') as fr:
  identifiers_org = pd.DataFrame.from_dict(json.load(fr)['payload']['namespaces'])

acceptable_prefixes.update(identifiers_org.prefix)

# N2T
with fsspec.open('https://n2t.net/e/n2t_full_prefixes.yaml') as fr:
  n2t = pd.DataFrame.from_dict(yaml.safe_load(fr), orient='index')

acceptable_schemes.update(n2t[n2t.type=='scheme'].index)
acceptable_prefixes.update(n2t[n2t.type=='commonspfx'].index)

# DRS
acceptable_drs = {'drs'}
acceptable_schemes.update(acceptable_drs)

# export
acceptable = dict(
  scheme=sorted(list(acceptable_schemes)),
  prefix=sorted(list(acceptable_prefixes)),
)

yaml.dump(acceptable, sys.stdout)
