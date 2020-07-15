import os
import json
import glob
from math import log10
from itertools import count
from jinja2 import Environment, FileSystemLoader

version = '0.0.1'
root_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
appyter_path = os.path.join(root_dir, 'appyters')
appyters = [
  dict(
    path=os.path.dirname(path),
    **json.load(open(path, 'r')),
  )
  for path in glob.glob(os.path.join(appyter_path, '*', 'appyter.json'))
]
env = Environment(loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), 'templates')))
template = env.get_template('docker-compose.yml.j2')
docker_compose = template.render(
  appyters=appyters,
  count=count,
  enumerate=enumerate,
  int=int,
  iter=iter,
  len=len,
  log10=log10,
  next=next,
  os=os,
  root_dir=root_dir,
  version=version,
)
print(docker_compose)
