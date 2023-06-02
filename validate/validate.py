import re
import json
import click
import shutil
import logging
import pathlib
import tempfile
import nbformat as nbf
import jsonschema
import urllib.request
from fsspec.core import url_to_fs
from appyter import __version__ as appyter_library_version
from appyter.ext.urllib import URI
from subprocess import Popen, PIPE

root_directory = pathlib.Path(__file__).parent.parent
appyters_directory = root_directory/'appyters'

def ensure_list(v):
  if v is None: return []
  elif type(v) == list: return v
  else: return [v]

def try_json_loads(s):
  try:
    return json.loads(s)
  except:
    return s

def validate_appyter(appyter, library_version=appyter_library_version, logger=logging.getLogger()):
  appyter_directory = appyters_directory/appyter
  logger.info("Checking for existing of files...")
  assert (appyter_directory).is_dir(), f"Appyter not found, missing appyters/{appyter} directory"
  assert (appyter_directory/'README.md').is_file(), f"Missing appyters/{appyter}/README.md"
  assert (appyter_directory/'appyter.json').is_file(), f"Missing appyters/{appyter}/appyter.json"
  #
  logger.info("Preparing temporary directory...")
  tmp_root = root_directory/'.tmp'/appyter
  tmp_root.mkdir(exist_ok=True, parents=True)
  tmp_directory = pathlib.Path(tempfile.mkdtemp(dir=tmp_root))
  #
  logger.info("Validating `{appyter}/appyter.json`...")
  with (appyter_directory/'appyter.json').open('r') as fr:
    config = json.load(fr)
  validator = jsonschema.Draft7Validator({
    '$ref': (root_directory/'schema'/'appyter-validator.json').absolute().as_uri(),
  })
  errors = [error.message for error in validator.iter_errors(config)]
  assert errors == [], '\n'.join(errors)
  #
  name = config['name']
  assert name == appyter, f"The directory should be named like appyter.json:.name"
  #
  if 'image' in config:
    image = config['image']
    if re.match(r'^https?://', image):
      _, _, image_name = image.rpartition('/')
      image_path = tmp_directory/image_name
      logger.warning("It is recommended to use a relative path instead of a url")
      _, response = urllib.request.urlretrieve(config['image'], filename=tmp_directory/image_name)
      assert response.get_content_maintype() == 'image', 'Expected image content'
    else:
      image_path = appyter_directory/'static'/image
    #
    from PIL import Image
    with Image.open(image_path, 'r') as img:
      assert img.size == (1280, 720), "Image should be 1280x720 px"
  else:
    logger.warning(f"`{appyter}/appyter.json` should have an 'image' defined...")
  #
  nbfile = config['appyter']['file']
  nbpath = appyter_directory/nbfile
  assert nbpath.is_file(), f"Missing {nbfile}"
  #
  logger.info(f"Checking notebook for issues..")
  try:
    with nbpath.open('r') as fr:
      nb = nbf.read(fr, as_version=4)
  except Exception as e:
    logger.error(f"{nbfile} is not valid json")
    raise e
  for cell in nb.cells:
    if cell['cell_type'] == 'code':
      assert not cell.get('execution_count'), "Please clear all notebook output & metadata"
      assert not cell.get('metadata'), "Please clear all notebook output & metadata"
      assert not cell.get('outputs'), "Please clear all notebook output & metadata"
  assert not nb['metadata'].get('widgets'), "Please clear all notebook output & metadata"
  assert not nb['metadata'].get('execution_info'), "Please clear all notebook output & metadata"
  #
  logger.info("Creating Dockerfile...")
  import sys; sys.path.insert(0, str(root_directory))
  from compose.build_dockerfile import prepare_appyter
  with (appyter_directory/'Dockerfile').open('w') as fw:
    print(prepare_appyter(str(appyter_directory), config), file=fw)
  #
  appyter_tag = f"maayanlab/appyter-{config['name'].lower()}:{config['version']}-{library_version}"
  logger.info("Building Dockerfile...")
  with Popen([
    'docker', 'build',
    '--build-arg', f"appyter_version=appyter[production]@git+https://github.com/Maayanlab/appyter@v{library_version}",
    '-t', appyter_tag,
    '.',
  ], cwd=appyter_directory, stdout=PIPE, stderr=sys.stderr) as p:
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
          if (appyter_directory/default_file_path).is_file():
            logger.info(f"Copying example file {default_file} from {default_file_path}...")
            shutil.copyfile(appyter_directory/default_file_path, tmp_directory/default_file)
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
    assert (tmp_directory/config['appyter']['file']).is_file(), f"nbconstruct output was not created"
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
  shutil.rmtree(tmp_root)
  #
  logger.info(f"Success!")

@click.command(help='Performs validation test for a single appyter')
@click.option('-v', '--verbose', count=True, default=0, help='How verbose this should be, more -v = more verbose')
@click.option('--library-version', envvar='LIBRARY_VERSION', default=appyter_library_version, type=str, help='The appyter library version to use')
@click.argument('appyter', type=str)
def cli(appyter, verbose=0, library_version=appyter_library_version):
  logging.basicConfig(level=30 - (verbose*10))
  logger = logging.getLogger(appyter)
  validate_appyter(appyter, library_version=library_version, logger=logger)

if __name__ == '__main__':
  try:
    from dotenv import load_dotenv
    load_dotenv()
  except ImportError:
    logging.warn("Install dotenv to load env from .env")
  #
  cli()
