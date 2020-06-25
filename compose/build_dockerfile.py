import os
from textwrap import dedent

def build_dockerfile(appyter_path, config):
  dockerfile_parts = ['''
    FROM ubuntu
  ''', '''
    ENV DEBIAN_FRONTEND "noninteractive"
    ENV TZ "America/New_York"
  ''', '''
    RUN set -x \\
        && echo "Preparing system..." \\
        && apt-get -y update \\
        && apt-get -y install git python3-pip python3-dev \\
        && rm -rf /var/lib/apt/lists/* \\
        && pip3 install --upgrade pip
  ''', '''
    RUN set -x \\
      && echo "Installing jupyter kernel..." \\
      && pip3 install ipykernel \\
      && python3 -m ipykernel install
  ''']
  if os.path.isfile(os.path.join(appyter_path, 'deps.txt')):
    dockerfile_parts.append('''
      ADD deps.txt /app/deps.txt
      RUN set -x \\
        && echo "Installing system dependencies from deps.txt..." \\
        && apt-get -y update \\
        && apt-get -y install $(grep -v '^#' /app/deps.txt) \\
        && rm -rf /var/lib/apt/lists/* \\
        && rm /app/deps.txt
    ''')
  if os.path.isfile(os.path.join(appyter_path, 'setup.R')):
    dockerfile_parts.append('''
      ADD setup.R /app/setup.R
      RUN set -x \\
        && echo "Installing R..." \\
        && apt-get -y update \\
        && apt-get -y install r-base \\
        && rm -rf /var/lib/apt/lists/* \\
        && echo "Setting up R with setup.R..." \\
        && R -e "source('/app/setup.R')" \\
        && rm /app/setup.R
    ''')
  if os.path.isfile(os.path.join(appyter_path, 'requirements.txt')):
    dockerfile_parts.append('''
      ADD requirements.txt /app/requirements.txt
      RUN set -x \\
        && echo "Installing python dependencies from requirements.txt..." \\
        && pip3 install -Ivr /app/requirements.txt \\
        && rm /app/requirements.txt
    ''')
  dockerfile_parts.append('''
    ARG appyter_version=git+git://github.com/Maayanlab/appyter.git
    RUN set -x \\
      && echo "Installing appyter..." \\
      && pip3 install -Iv ${appyter_version}
  ''')
  dockerfile_parts.append('''
    WORKDIR /app
    EXPOSE 80
    VOLUME /app/data
  ''')
  dockerfile_parts.append('''
    ENV PREFIX "/"
    ENV HOST "0.0.0.0"
    ENV PORT "80"
    ENV DEBUG "false"
  ''')
  dockerfile_parts.append('''
    COPY . /app
  ''')
  dockerfile_parts.append('''
    RUN set -x \\
      && echo "Overriding appyter templates..." \\
      && python3 /app/merge_j2.py /app /app/override /app
  ''')
  dockerfile_parts.append(f'''
    CMD [ "appyter", "--profile={config['appyter'].get('profile', 'default')}", "{config['appyter']['file']}" ]
  ''')
  return '\n\n'.join(map(str.strip, map(dedent, dockerfile_parts)))

if __name__ == '__main__':
  import sys
  import json
  import shutil
  appyter = sys.argv[1]
  appyter_path = os.path.join(os.path.dirname(__file__), '..', 'appyters', appyter)
  config = json.load(open(os.path.join(appyter_path, 'appyter.json'), 'r'))
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
  print(build_dockerfile(appyter_path, config))
