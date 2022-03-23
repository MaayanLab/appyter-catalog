import sh
import re
import os
import json
import glob
from dotenv import load_dotenv
load_dotenv(os.path.join(os.path.dirname(os.path.dirname(__file__)), '.env'))

root_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
docker_org = os.environ['DOCKER_REGISTRY']
library_version = os.environ['appyter_version']

appyter_path = os.path.join(root_dir, 'appyters')
def get_appyters(appyter_path):
  for path in map(os.path.dirname, glob.glob(os.path.join(appyter_path, '*', 'appyter.json'))):
    appyter = json.load(open(os.path.join(path, 'appyter.json'), 'r'))
    yield dict(
      appyter,
      path=path,
      long_description=open(os.path.join(path, 'README.md'), 'r').read(),
      # find the oldest commit containing the appyter's appyter.json (follow for detecting renames)
      creation_timestamp=str(sh.tail(
        sh.git.log(
          '--follow', r'--pretty=format:%aI', '--', os.path.join(path, 'appyter.json'),
          _tty_out=False,
        ),
        '-n1'
      )).strip(),
      # find the most recent commit containing the appyter's directory
      update_timestamp=str(sh.head(
        sh.git.log(
          r'--pretty=format:%aI', '--', path,
          _tty_out=False,
        ),
        '-n1'
      )).strip(),
    )

appyters = list(get_appyters(appyter_path))

config = json.load(open(os.path.join(os.path.dirname(__file__), 'templates', 'appyters.json'), 'r'))
config['appyters'] = appyters
config['docker_org'] = docker_org
config['library_version'] = os.environ['appyter_version']
if os.environ.get('keycloak_url'):
  config['keycloak'] = dict(
    url=os.environ['keycloak_url'],
    realm=os.environ['keycloak_realm'],
    clientId=os.environ['keycloak_client_id'],
  )
if os.environ.get('ga_id'):
  config['ga_id'] = os.environ['ga_id']

if __name__ == '__main__':
  print(json.dumps(config))
