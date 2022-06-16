import requests

def enrichr_libraries():
  return sorted([stat['libraryName'] for stat in requests.get('https://maayanlab.cloud/Enrichr/datasetStatistics').json()['statistics']])
