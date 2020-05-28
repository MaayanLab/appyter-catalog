import re
import os
import json
import glob

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

proxy_environment = '\n'.join(f"""
      - nginx_proxy_{n:03}=/{appyter['name']}(/.*) http://{appyter['name'].lower()}:80/{appyter['name']}$$1
""".strip('\n') for n, appyter in enumerate(appyters)).strip('\n')

proxy_service = f"""
  proxy:
    image: maayanlab/proxy:1.2.0
    environment:
{proxy_environment}
      - nginx_server_name=${{nginx_server_name}}
      - nginx_ssl=${{nginx_ssl}}
      - nginx_ssl_letsencrypt=${{nginx_ssl}}
      - letsencrypt_email=${{letsencrypt_email}}
      - nginx_proxy_{len(appyters):03}=(/.*) http://app:80$$1
    ports:
      - 80:80
      - 443:443
    volumes:
      - ./data/letsencrypt/:/etc/letsencrypt/
""".strip('\n')

docker_compose_services = '\n'.join(f"""
  {appyter['name'].lower()}:
    build:
      context: {os.path.relpath(appyter['path'], root_dir)}
      dockerfile: Dockerfile
      args:
        - appyter_version=${{appyter_version}}
    image: maayanlab/appyter-{appyter['name'].lower()}:{appyter['version']}
    environment:
      - PREFIX=/{appyter['name']}/
    volumes:
      - ./data/{appyter['name'].lower()}/:/app/data
""".strip('\n') for appyter in appyters)

docker_compose = f"""
version: '3'
services:
{proxy_service}
  app:
    build: app
    image: maayanlab/appyters:{version}
{docker_compose_services}
""".strip('\n')

print(docker_compose)
