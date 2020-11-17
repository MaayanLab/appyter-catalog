import os
import re
import json
import glob
import yaml
from io import StringIO
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
env = Environment(
  extensions=['jinja2.ext.with_'],
  loader=FileSystemLoader(os.path.join(os.path.dirname(__file__), 'templates')),
)
template = env.get_template('Chart.yaml.j2')
chart_files_spec = template.render(
  appyters=appyters,
  version=version,
)
chart_files = re.compile(r'\n+---', re.MULTILINE).split(chart_files_spec)
# load first two files (Chart.yaml & questions.yaml)
chart_desc, chart_questions = yaml.load(StringIO(chart_files[0])), yaml.load(StringIO(chart_files[1]))
chart_root = os.path.join(root_dir, 'charts', chart_desc['name'], chart_desc['version'])
# create values.yaml from questions.yaml
os.makedirs(chart_root, exist_ok=True)
with open(os.path.join(chart_root, 'values.yaml'), 'w') as fw:
  for question in chart_questions['questions']:
    print(
      f"# {question['description']} ({question['group']})",
      f"{'' if question['required'] else '#'}{question['variable']}: \"{question['default']}\"",
      sep='\n',
      file=fw,
    )
# create templates
for chart_file in chart_files:
  m = re.compile(r'\n*# Source: ([^\n]+)\n(.+)', re.DOTALL | re.MULTILINE).match(chart_file)
  chart_file_path, chart_file_content = m.group(1), m.group(2)
  chart_file_full_path = os.path.join(chart_root, chart_file_path)
  os.makedirs(os.path.dirname(chart_file_full_path), exist_ok=True)
  with open(chart_file_full_path, 'w') as fw:
    fw.write(chart_file_content)
