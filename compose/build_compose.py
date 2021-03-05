import click

@click.command(help='Build the docker-compose.yml file')
@click.option('--tls', default=False, type=bool, is_flag=True, help='Whether or not to build the docker-compose.yml with tls support')
def build_compose(tls):
  import os
  import json
  import glob
  from math import log10
  from itertools import count
  from jinja2 import Environment, FileSystemLoader
  #
  root_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
  version = open(os.path.join(root_dir, 'VERSION'), 'r').read().strip()
  appyter_path = os.path.join(root_dir, 'appyters')
  appyters = [
    dict(
      path=os.path.dirname(path),
      **json.load(open(path, 'r')),
    )
    for path in glob.glob(os.path.join(appyter_path, '*', 'appyter.json'))
  ]
  env = Environment(
    extensions=['jinja2.ext.with_'],
    loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), 'templates')),
  )
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
    tls=tls,
  )
  print(docker_compose)


if __name__ == '__main__':
  build_compose()
