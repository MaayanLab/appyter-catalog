import os
import re
import sys
import json
import click
import shutil
import logging
import nbformat as nbf
import traceback
import jsonschema
import urllib.request
from fsspec.core import url_to_fs
from appyter import __version__ as appyter_library_version
from appyter.ext.urllib import URI
from PIL import Image
from subprocess import Popen, PIPE

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

def try_json_loads(s):
  import json
  try:
    return json.loads(s)
  except:
    s

def ensure_list(v):
  if v is None: return []
  elif type(v) == list: return v
  else: return [v]

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

def validate_appyter(appyter, library_version=appyter_library_version, logger=logging.getLogger()):
  logger.info("Preparing temporary directory...")
  tmp_directory = os.path.realpath('.tmp')
  os.makedirs(tmp_directory, exist_ok=True)
  #
  logger.info("Checking for existing of files...")
  assert os.path.isfile(os.path.join('appyters', appyter, 'README.md')), f"Missing appyters/{appyter}/README.md"
  assert os.path.isfile(os.path.join('appyters', appyter, 'appyter.json')), f"Missing appyters/{appyter}/appyter.json"
  #
  logger.info("Validating `{appyter}/appyter.json`...")
  config = json.load(open(os.path.join('appyters', appyter, 'appyter.json'), 'r'))
  validator = jsonschema.Draft7Validator({
    '$ref': f"file:///{os.path.realpath(os.path.join(os.path.dirname(__file__), '..', 'schema', 'appyter-validator.json'))}",
  })
  errors = [error.message for error in validator.iter_errors(config)]
  assert errors == [], '\n'.join(errors)
  #
  name = config['name']
  assert name == appyter, f"The directory should be named like `name`"
  #
  if 'image' in config:
    image = config['image']
    if re.match(r'^https?://', image):
      image_name = os.path.basename(image)
      image_path = os.path.join(tmp_directory, image_name)
      logger.warning("It is recommended to use a relative path instead of a url")
      _, response = urllib.request.urlretrieve(config['image'], filename=os.path.join(tmp_directory, image_name))
      assert response.get_content_maintype() == 'image', 'Expected image content'
    else:
      image_path = f"appyters/{appyter}/static/{image}"
    #
    with Image.open(image_path, 'r') as img:
      assert img.size == (1280, 720), "Image should be 1280x720 px"
  else:
    logger.warning(f"`{appyter}/appyter.json` should have an 'image' defined...")
  #
  nbfile = config['appyter']['file']
  nbpath = os.path.join('appyters', appyter, nbfile)
  #
  logger.info(f"Checking notebook for issues..")
  nb = nbf.read(open(nbpath, 'r'), as_version=4)
  for cell in nb.cells:
    if cell['cell_type'] == 'code':
      assert not cell.get('execution_count'), "Please clear all notebook output & metadata"
      assert not cell.get('metadata'), "Please clear all notebook output & metadata"
      assert not cell.get('outputs'), "Please clear all notebook output & metadata"
  assert not nb['metadata'].get('widgets'), "Please clear all notebook output & metadata"
  assert not nb['metadata'].get('execution_info'), "Please clear all notebook output & metadata"
  #
  logger.info(f"Preparing docker to run `{nbfile}`...")
  assert os.path.isfile(nbpath), f"Missing {nbpath}"
  try:
    json.load(open(nbpath, 'r'))
  except Exception as e:
    logger.error(f"{nbfile} is not valid json")
    raise e
  #
  assert not os.path.isfile(os.path.join('appyters', appyter, 'Dockerfile')), 'Custom Dockerfiles are no longer supported'
  logger.info("Creating Dockerfile...")
  import sys; sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
  from compose.build_dockerfile import prepare_appyter
  with open(os.path.join('appyters', appyter, 'Dockerfile'), 'w') as fw:
    print(prepare_appyter(os.path.join('appyters', appyter), config), file=fw)
  #
  logger.info("Building Dockerfile...")
  appyter_tag = f"maayanlab/appyter-{config['name'].lower()}:{config['version']}-{library_version}"
  with Popen([
    'docker', 'build',
    '--build-arg', f"appyter_version=appyter[production]@git+https://github.com/Maayanlab/appyter@v{library_version}",
    '-t', appyter_tag,
    '.',
  ], cwd=os.path.join('appyters', appyter), stdout=PIPE, stderr=sys.stderr) as p:
    for line in filter(None, map(str.strip, map(bytes.decode, p.stdout))):
      logger.debug(f"`docker build .`: {line}")
    assert p.wait() == 0, '`docker build .` command failed'
  #
  logger.info("Inspecting appyter...")
  with Popen([
    'docker', 'run',
    '-e', 'APPYTER_PREFIX=', # hotfix because prefix is baked-in and necessary at production initialization time, but should be empty here
    appyter_tag,
    'appyter', 'nbinspect',
    nbfile,
  ], stdout=PIPE, stderr=sys.stderr) as p:
    nbinspect_output = p.stdout.read().decode().strip()
    logger.info(f"`appyter nbinspect {nbfile}`: {nbinspect_output})")
    assert p.wait() == 0, f"`appyter nbinspect {nbfile}` command failed"
  #
  inspect = json.loads(nbinspect_output)
  field_args = {
    field['args']['name']: field['args']
    for field in inspect
  }
  assert len(field_args) == len(inspect), "Some of your fields weren't captured, there might be duplicate `name`s"
  #
  logger.info("Preparing defaults...")
  default_args = {
    field_name: field.get('default')
    for field_name, field in field_args.items()
  }
  file_fields = {
    field['args']['name']
    for field in inspect
    if field['field'] in {'MultiFileField', 'FileField'}
  }
  early_stopping = False
  for file_field in file_fields:
    field_examples = field_args[file_field].get('examples', {})
    default_files = ensure_list(default_args[file_field])
    if default_files:
      for default_file in default_files:
        if default_file in field_examples:
          default_file_path = field_examples[default_file].lstrip('/')
          if os.path.exists(os.path.join('appyters', appyter, default_file_path)):
            logger.info(f"Copying example file {default_file} from {default_file_path}...")
            shutil.copyfile(os.path.join('appyters', appyter, default_file_path), os.path.join(tmp_directory, default_file))
          else:
            logger.info(f"Using example file {default_file} from {default_file_path}...")
            try:
              fs, fs_path = url_to_fs(default_file_path)
              assert fs.exists(fs_path), "file does not exist"
            except AssertionError as e:
              logger.warning(f"example file {default_file} from {default_file_path} resulted in error {str(e)}.")
              early_stopping = True
            else:
              uri_parsed = URI(default_file_path)
              default_args[file_field] = str(uri_parsed.with_fragment(default_file))
        else:
          logger.warning(f"default file isn't in examples, we won't know how to get it if it isn't available in the image")
    else:
      logger.warning(f"no default file is provided")
  #
  if early_stopping:
    logger.warning(f"Stopping early as a download requires manual intervention.")
    return
  logger.info(f"Fixing permissions...")
  assert Popen(['chmod', '-R', '777', tmp_directory]).wait() == 0, f"ERROR: Changing permissions failed"
  logger.info(f"Constructing default notebook from appyter...")
  with Popen([
    'docker', 'run',
    '-v', f"{tmp_directory}:/data",
    "-i", appyter_tag,
    'appyter', 'nbconstruct',
    f"--output=/data/{nbfile}",
    nbfile,
  ], stdin=PIPE, stdout=PIPE, stderr=sys.stderr) as p:
    procLogger = logger.getChild(f"appyter nbconstruct {nbfile}")
    procLogger.info(f"`< {default_args}")
    stdout, _ = p.communicate(json.dumps(default_args).encode())
    for line in filter(None, map(str.strip, stdout.decode().splitlines())):
      procLogger.debug(f"{line}")
    assert p.wait() == 0, f"`appyter nbconstruct {nbfile}` command failed"
    assert os.path.exists(os.path.join(tmp_directory, config['appyter']['file'])), f"nbconstruct output was not created"
  #
  logger.info(f"Executing default notebook with appyter...")
  with Popen([
    'docker', 'run',
    '-v', f"{tmp_directory}:/data",
    '-e', 'PYTHONPATH=/app',
    '--device', '/dev/fuse',
    '--cap-add', 'SYS_ADMIN',
    '--security-opt', 'apparmor:unconfined',
    appyter_tag,
    'appyter', 'nbexecute',
    '--fuse=true',
    f"/data/{nbfile}",
  ], stdout=PIPE, stderr=sys.stderr) as p:
    procLogger = logger.getChild(f"appyter nbexecute {nbfile}")
    last_msg = None
    for msg in map(try_json_loads, p.stdout):
      if type(msg) == dict and msg['type'] == 'error':
        procLogger.error(f"{json.dumps(last_msg)}")
        procLogger.error(f"{json.dumps(msg)}")
        raise Exception(f"error {msg.get('data')}")
      else:
        procLogger.debug(f"{json.dumps(msg)}")
        last_msg = msg
    assert p.wait() == 0, f"`appyter nbexecute {nbfile}` command failed"
  #
  logger.info(f"Success!")

@click.command(help='Performs validation tests for all appyters that were changed when diffing against origin/main')
@click.option('-v', '--verbose', count=True, default=0, help='How verbose this should be, more -v = more verbose')
@click.option('--github-action', default=False, type=bool, is_flag=True, help='Use for receiving json on stdin from github actions')
@click.option('--library-version', envvar='LIBRARY_VERSION', default=appyter_library_version, type=str, help='The appyter library version to use')
def validate_merge(github_action=False, verbose=0, library_version=appyter_library_version):
  logging.basicConfig(level=30 - (verbose*10))
  valid = True
  for appyter in get_changed_appyters(github_action):
    logger = BufferedLog(logging.getLogger(appyter))
    if not os.path.exists(os.path.join('appyters', appyter)):
      logger.info(f"{appyter} directory no longer exists, ignoring")
      continue
    elif not os.path.isdir(os.path.join('appyters', appyter)):
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
  validate_merge()
