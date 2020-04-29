import os
from textwrap import dedent

def build_dockerfile(template_path, config):
  dockerfile_parts = ['FROM ubuntu']
  if os.path.isfile(os.path.join(template_path, 'deps.txt')):
    dockerfile_parts.append('''
      ADD deps.txt /app/deps.txt
      RUN set -x \\
        && echo "Installing system dependencies from deps.txt..." \\
        && apt-get -y update \\
        && apt-get -y install git $(grep -v '^#' /app/deps.txt) \\
        && rm /app/deps.txt
    ''')
  if os.path.isfile(os.path.join(template_path, 'setup.R')):
    dockerfile_parts.append('''
      ADD setup.R /app/setup.R
      RUN set -x \\
        && echo "Installing R..." \\
        && apt-get install -y r-base \\
        && echo "Installing R dependencies from setup.R..." \\
        && R -e "source('/app/setup.R')" \\
        && rm /app/setup.R
    ''')
  if os.path.isfile(os.path.join(template_path, 'requirements.txt')):
    dockerfile_parts.append('''
      ADD requirements.txt /app/requirements.txt
      RUN set -x \\
        && echo "Installing python..." \\
        && apt-get install -y  python3-pip python3-dev \\
        && pip3 install --upgrade pip \\
        && echo "Installing python dependencies from requirements.txt..." \\
        && pip3 install -Ivr /app/requirements.txt \\
        && rm /app/requirements.txt
    ''')
  dockerfile_parts.append('''
    WORKDIR /app
    EXPOSE 80
  ''')
  dockerfile_parts.append('''
    ENV PREFIX "/"
    ENV HOST "0.0.0.0"
    ENV PORT "80"
  ''')
  dockerfile_parts.append('''
    COPY . /app
  ''')
  dockerfile_parts.append(f'''
    CMD [ "jupyter-template", "--profile={config['template'].get('profile', 'default')}", "{config['template']['file']}" ]
  ''')
  return '\n\n'.join(map(str.strip, map(dedent, dockerfile_parts)))

if __name__ == '__main__':
  import sys
  import json
  template = sys.argv[1]
  template_path = os.path.join(os.path.dirname(__file__), '..', 'templates', template)
  config = json.load(open(os.path.join(template_path, 'template.json'), 'r'))
  print(build_dockerfile(template_path, config))
