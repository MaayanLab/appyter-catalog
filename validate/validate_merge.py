import sys
import json
import click
import pathlib
import logging
import traceback
from appyter import __version__ as appyter_library_version
from subprocess import Popen, PIPE

sys.path.insert(0, str(pathlib.Path(__file__).parent)); from validate import validate_appyter
root_dir = pathlib.Path(__file__).parent.parent
appyters_dir = root_dir/'appyters'

class BufferedLog:
  ''' A logger wrapper which saves a buffer and dumps debugging messages in the buffer when necessary
  (after an exception, we get more detailed information leading up to the exception)
  '''
  def __init__(self, logger=logging.getLogger(), buffer_size=25):
    self.logger = logger
    self.buffer = []
    self.buffer_size = buffer_size
  
  def log(self, level, msg):
    self.logger.log(level, msg)
    self.buffer.append((level, msg))
    while len(self.buffer) >= self.buffer_size:
      self.buffer.pop(0)

  def debug(self, msg):
    self.log(logging.DEBUG, msg)
  def info(self, msg):
    self.log(logging.INFO, msg)
  def warning(self, msg):
    self.log(logging.WARNING, msg)
  def error(self, msg):
    self.log(logging.ERROR, msg)

  def getChild(self, child):
    return BufferedLog(logger=self.logger.getChild(child), buffer_size=self.buffer_size)

  def replay_buffer(self, replay_level=logging.DEBUG):
    current_level = self.logger.getEffectiveLevel()
    self.logger.setLevel(replay_level)
    for level, msg in self.buffer:
      if level < current_level and level >= replay_level:
        self.logger.log(level, msg)
    self.logger.setLevel(current_level)

def get_changed_appyters_gh_action():
  ops = {
    'added': 'A',
    'renamed': 'R',
    'modified': 'M',
    'removed': 'D',
  }
  for record in json.load(sys.stdin):
    yield dict(
      filename=record['filename'],
      status=ops[record['status']],
      prev=record.get('previous_filename'),
    )

def get_changed_appyters_git():
  with Popen(['git', 'diff', '--name-status', 'origin/main'], stdout=PIPE, stderr=sys.stderr) as p:
    for line in filter(None, map(str.strip, map(bytes.decode, p.stdout))):
      op, filename, *rest = line.split()
      prev = rest[-1] if rest else None
      yield dict(
        filename=filename,
        status=op[0],
        prev=prev,
      )

def dict_groupby(iterator, grouper, keyfunc):
  groups = {}
  for element in iterator:
    group = grouper(element)
    if group not in groups:
      groups[group] = {}
    key = keyfunc(element)
    groups[group][key] = element
  return groups

def get_changed_appyters(github_action):
  changed_files = list(get_changed_appyters_gh_action() if github_action else get_changed_appyters_git())
  current = dict_groupby(
    [
      record
      for record in changed_files
      if record['filename'].startswith('appyters/')
    ] + [
      dict(
        filename=record['prev'],
        status='D',
        prev=None,
      )
      for record in changed_files
      if record.get('prev') and record['prev'].startswith('appyters/')
    ],
    grouper=lambda record: record['filename'].split('/', maxsplit=3)[1],
    keyfunc=lambda record: record['filename'].split('/', maxsplit=3)[2],
  )
  logging.debug(f"{current=}")
  for appyter, changes in reversed(list(current.items())):
    logging.getLogger(appyter).info(f"Changes detected")
    appyterJson = changes.get('appyter.json')
    assert appyterJson, f"Expected appyter.json modification for {appyter}"
    yield appyter

@click.command(help='Performs validation tests for all appyters that were changed when diffing against origin/main')
@click.option('-v', '--verbose', count=True, default=0, help='How verbose this should be, more -v = more verbose')
@click.option('--github-action', default=False, type=bool, is_flag=True, help='Use for receiving json on stdin from github actions')
@click.option('--library-version', envvar='LIBRARY_VERSION', default=appyter_library_version, type=str, help='The appyter library version to use')
def cli(github_action=False, verbose=0, library_version=appyter_library_version):
  logging.basicConfig(level=30 - (verbose*10))
  valid = True
  for appyter in get_changed_appyters(github_action):
    logger = BufferedLog(logging.getLogger(appyter))
    if not (appyters_dir/appyter).exists():
      logger.info(f"{appyter} directory no longer exists, ignoring")
      continue
    elif not (appyters_dir/appyter).is_dir():
      logger.info(f"{appyter} is not a directory, ignoring")
      continue
    #
    try:
      validate_appyter(appyter, library_version=library_version, logger=logger)
    except Exception as e:
      logger.replay_buffer()
      logger.error(str(e))
      logger.error(traceback.format_exc())
      valid = False
  #
  if valid:
    sys.exit(0)
  else:
    sys.exit(1)

if __name__ == '__main__':
  try:
    from dotenv import load_dotenv
    load_dotenv()
  except ImportError:
    logging.warn("Install dotenv to load env from .env")
  #
  cli()
