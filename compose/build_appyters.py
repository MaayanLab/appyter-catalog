import sh
import os
import json
import glob
from dotenv import load_dotenv
load_dotenv(os.path.join(os.path.dirname(os.path.dirname(__file__)), '.env'))

root_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
docker_org = os.environ['DOCKER_REGISTRY']
library_version = os.environ['LIBRARY_VERSION']

appyter_path = os.path.join(root_dir, 'appyters')
def get_appyters(appyter_path):
  for path in map(os.path.dirname, glob.glob(os.path.join(appyter_path, '*', 'appyter.json'))):
    appyter = json.load(open(os.path.join(path, 'appyter.json'), 'r'))
    yield dict(
      appyter,
      path=path,
      long_description=open(os.path.join(path, 'README.md'), 'r').read(),
      # find the oldest commit containing the appyter's README (follow for detecting renames)
      creation_timestamp=str(sh.tail(
        '-n1',
        _in=sh.git.log(
          '--follow', r'--pretty=format:%aI', '--', os.path.join(path, 'README.md'),
          _tty_out=False,
        ),
      )).strip(),
      # find the most recent commit containing the appyter's directory
      update_timestamp=str(sh.head(
        '-n1',
        _in=sh.git.log(
          r'--pretty=format:%aI', '--', path,
          _tty_out=False,
        ),
      )).strip(),
    )

appyters = list(get_appyters(appyter_path))

config = json.load(open(os.path.join(os.path.dirname(__file__), 'templates', 'appyters.json'), 'r'))
config['appyters'] = appyters
config['docker_org'] = docker_org
config['library_version'] = os.environ['LIBRARY_VERSION']
if os.environ.get('KEYCLOAK_URL'):
  config['keycloak'] = dict(
    url=os.environ['KEYCLOAK_URL'],
    realm=os.environ['KEYCLOAK_REALM'],
    clientId=os.environ['KEYCLOAK_CLIENT_ID'],
  )
if os.environ.get('GA_ID'):
  config['ga_id'] = os.environ['GA_ID']

if __name__ == '__main__':
  print(json.dumps(config))
