import os
import json
from jinja2 import Environment, FileSystemLoader

def build_dockerfile(appyter_path, config):
  env = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), 'templates')))
  template = env.get_template('Dockerfile.j2')
  dockerfile = template.render(
    appyter_path=appyter_path,
    config=config,
    os=os,
  )
  return dockerfile

def prepare_appyter(appyter_path, config):
  import shutil
  override_path = os.path.join(appyter_path, 'override')
  if os.path.exists(override_path):
    shutil.rmtree(override_path)
  shutil.copytree(
    os.path.join(os.path.dirname(__file__), '..', 'override'),
    os.path.join(appyter_path, 'override'),
  )
  shutil.copy(
    os.path.join(os.path.dirname(__file__), 'merge_j2.py'),
    os.path.join(appyter_path, 'merge_j2.py')
  )
  return build_dockerfile(appyter_path, config)

if __name__ == '__main__':
  import sys
  appyter = sys.argv[1]
  appyter_path = os.path.join(os.path.dirname(__file__), '..', 'appyters', appyter)
  config = json.load(open(os.path.join(appyter_path, 'appyter.json'), 'r'))
  print(prepare_appyter(appyter_path, config))
