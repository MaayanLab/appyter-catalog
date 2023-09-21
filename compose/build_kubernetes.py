import click

@click.command(help='Build the kubernetes.yml file')
@click.option('--tls', default=False, type=bool, is_flag=True, help='Whether or not to build the kubernetes.yml with tls support')
@click.option('--minio', default=False, type=bool, is_flag=True, help='Whether or not to build the kubernetes.yml with a minio deployment')
@click.option('--auth', default=False, type=bool, is_flag=True, help='Whether or not to build the kubernetes.yml with a authentication support')
@click.option('--aws-proxy', default=False, type=bool, is_flag=True, help='An additional proxy to aws for legacy reasons, maps /storage => aws bucket')
def build_kubernetes(tls, minio, auth, aws_proxy):
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
    loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), 'templates')),
  )
  template = env.get_template('kubernetes.yml.j2')
  kubernetes = template.render(
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
    auth=auth,
    minio=minio,
    aws_proxy=aws_proxy,
  )
  print(kubernetes)


if __name__ == '__main__':
  from dotenv import load_dotenv
  load_dotenv()
  build_kubernetes()
