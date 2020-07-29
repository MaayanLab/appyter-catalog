import sh
import re
import os
import json
import glob

root_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), '..'))
appyter_path = os.path.join(root_dir, 'appyters')
def get_appyters(appyter_path):
  for path in map(os.path.dirname, glob.glob(os.path.join(appyter_path, '*', 'appyter.json'))):
    appyter = json.load(open(os.path.join(path, 'appyter.json'), 'r'))
    if appyter.get('public') == False:
      continue
    yield dict(
      appyter,
      path=path,
      long_description=open(os.path.join(path, 'README.md'), 'r').read(),
      # find the oldest commit containing the appyter's appyter.json (follow for detecting renames)
      creation_timestamp=str(sh.tail(
        sh.git.log(
          '--follow', r'--pretty=format:%aI', '--', os.path.join(path, 'appyter.json'),
          _tty_out=False,
        ),
        '-n1'
      )).strip(),
      # find the most recent commit containing the appyter's directory
      update_timestamp=str(sh.head(
        sh.git.log(
          r'--pretty=format:%aI', '--', path,
          _tty_out=False,
        ),
        '-n1'
      )).strip(),
    )

appyters = list(get_appyters(appyter_path))

if __name__ == '__main__':
  print(json.dumps(appyters))
