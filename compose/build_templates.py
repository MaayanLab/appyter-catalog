import re
import os
import json
import glob

root_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
template_path = os.path.join(root_dir, 'templates')
templates = [
  dict(
    path=path,
    long_description=open(os.path.join(path, 'README.md'), 'r').read(),
    **json.load(open(os.path.join(path, 'template.json'), 'r')),
  )
  for path in map(os.path.dirname, glob.glob(os.path.join(template_path, '*', 'template.json')))
]

if __name__ == '__main__':
  print(json.dumps(templates))
