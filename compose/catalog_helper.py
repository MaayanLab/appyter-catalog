import os
import re
import click
import shutil

@click.group()
def cli():
  pass

@cli.command(name='merge-j2')
@click.argument('primary', type=click.Path(exists=True, dir_okay=True, file_okay=False))
@click.argument('override', type=click.Path(exists=True, dir_okay=True, file_okay=False))
@click.argument('merged', type=click.Path(dir_okay=True, file_okay=False))
def merge_j2_cli(primary, override, merged):
  merge_j2_directories(primary, override, merged)

@cli.command(name='insert-info')
@click.option('-i', '--info', type=click.File('r'), default='-', help='File containing info json')
@click.argument('ipynb', type=click.Path(file_okay=True, dir_okay=False))
def insert_info_cli(ipynb, info):
  import json
  insert_info(ipynb, json.load(info))

def insert_info(ipynb, info):
  ''' Given an ipynb, insert { 'metadata': { 'appyter': 'info': info } } }
  '''
  import nbformat as nbf
  with open(ipynb, 'r') as fr:
    nb = nbf.read(fr, as_version=4)
  if 'appyter' not in nb.metadata:
    nb.metadata['appyter'] = {}
  nb.metadata['appyter']['info'] = info
  with open(ipynb, 'w') as fw:
    nbf.write(nb, fw)

def merge_j2(*j2s):
  ''' Given a set of independent jinja2 templates, under certain conditions, we can merge the two together into one template.
  '''
  assert len(j2s) > 0, 'Nothing to merge'
  # ensure extends are the same and present
  extends = set()
  for j2 in j2s:
    m = re.search(r'\{% *extends (?P<extend>.+?) *%\}', j2)
    if m:
      extends.add(m.group('extend'))
    else:
      extends.add(None)
  assert len(extends) == 1, f"Expected extends to be the same, got {extends}"
  extend = next(iter(extends))
  assert extend != None, f"Can only merge extended j2s"
  # locate blocks
  blocks = {}
  for ind, j2 in enumerate(j2s):
    # WARNING: this would fail on nested blocks
    for black_match in re.finditer(r'\{% *block *(?P<block_name>.+?) *%\}(?P<block_content>(.|\r|\n)+?)\{% *endblock *%\}', j2, re.MULTILINE):
      block_name = black_match.group('block_name')
      block_content = black_match.group('block_content')
      if block_name not in blocks:
        blocks[block_name] = {}
      # capture location of super
      super_match = re.search(r'\{\{ *super\(\) *\}\}', block_content)
      if super_match:
        start_super, end_super = super_match.span()
        blocks[block_name][ind] = {
          'pre_super': block_content[:start_super],
          'post_super': block_content[end_super:],
        }
      else:
        blocks[block_name][ind] = {
          'no_super': block_content,
        }
  # merge blocks, combining contents of blocks relative to super
  merged_blocks = {}
  for block_name, block_contents in blocks.items():
    for ind, block_content in block_contents.items():
      if block_name not in merged_blocks:
        merged_blocks[block_name] = {}
      merged_blocks[block_name]['pre_super'] = '\n'.join(filter(None, [
        merged_blocks[block_name].get('pre_super', '').strip(),
        block_content.get('pre_super', '').strip(),
      ])).strip()
      merged_blocks[block_name]['no_super'] = '\n'.join(filter(None, [
        merged_blocks[block_name].get('no_super', '').strip(),
        block_content.get('no_super', '').strip(),
      ])).strip()
      merged_blocks[block_name]['post_super'] = '\n'.join(filter(None, [
        merged_blocks[block_name].get('post_super', '').strip(),
        block_content.get('post_super', '').strip(),
      ])).strip()
  # merge into final template
  merged = '\n\n'.join([
    f"{{% extends {extend} %}}",
    *[
      '\n'.join(filter(None, [
        f"{{% block {block_name} %}}",
        block_content['pre_super'].strip(),
        block_content['no_super'].strip() if block_content.get('no_super') else '{{ super() }}',
        block_content['post_super'].strip(),
        f"{{% endblock %}}"
      ]))
      for block_name, block_content in merged_blocks.items()
    ],
  ])
  return merged


def merge_j2_directories(primary_dir, override_dir, merged_dir):
  ''' Given a primary directory and an override directory, recursively
  copy over or merge and overrides
  '''
  if primary_dir != merged_dir:
    for dirpath, dirnames, filenames in os.walk(primary_dir):
      for filename in filenames:
        # handle each override file
        file_dir = os.path.relpath(dirpath, primary_dir)
        file_path = os.path.join(file_dir, filename)
        # copy primary into directory
        os.makedirs(os.path.join(merged_dir, file_dir), exist_ok=True)
        shutil.copyfile(
          os.path.join(primary_dir, file_path),
          os.path.join(merged_dir, file_path)
        )
  #
  for dirpath, dirnames, filenames in os.walk(override_dir):
    for filename in filenames:
      # handle each override file
      file_dir = os.path.relpath(dirpath, override_dir)
      file_path = os.path.join(file_dir, filename)
      if os.path.exists(os.path.join(primary_dir, file_path)):
        # join j2 and update
        if os.path.splitext(filename)[1] == '.j2':
          merged = merge_j2(
            open(os.path.join(primary_dir, file_path), 'r').read(),
            open(os.path.join(override_dir, file_path), 'r').read(),
          )
          os.makedirs(os.path.join(merged_dir, file_dir), exist_ok=True)
          with open(os.path.join(merged_dir, file_path), 'w') as fw:
            fw.write(merged)
        else:
          raise Exception(f"Unable to merge {file_path}")
      else:
        # copy override into directory
        os.makedirs(os.path.join(merged_dir, file_dir), exist_ok=True)
        shutil.copyfile(
          os.path.join(override_dir, file_path),
          os.path.join(merged_dir, file_path)
        )

if __name__ == '__main__':
  cli()
