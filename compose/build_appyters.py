import re
import os
import json
import glob

root_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
appyter_path = os.path.join(root_dir, 'appyters')
appyters = [
  dict(
    path=path,
    long_description=open(os.path.join(path, 'README.md'), 'r').read(),
    **json.load(open(os.path.join(path, 'appyter.json'), 'r')),
  )
  for path in map(os.path.dirname, glob.glob(os.path.join(appyter_path, '*', 'appyter.json')))
]

if __name__ == '__main__':
  print(json.dumps(appyters))
