import os
import shutil

def test_merge_j2():
  from compose.merge_j2 import merge_j2
  test_0 = open(os.path.join(os.path.dirname(__file__), 'test_merge_j2_0.j2'), 'r').read()
  test_1 = open(os.path.join(os.path.dirname(__file__), 'test_merge_j2_1.j2'), 'r').read()
  expectation = open(os.path.join(os.path.dirname(__file__), 'test_merge_j2_01.j2'), 'r').read()
  result = merge_j2(test_0, test_1)
  print(result)
  assert expectation == result

def test_merge_j2_directories():
  from compose.merge_j2 import merge_j2_directories
  test_primary = os.path.join(os.path.dirname(__file__), 'primary')
  test_override = os.path.join(os.path.dirname(__file__), 'override')
  test_merged = os.path.join(os.path.dirname(__file__), 'merged')
  shutil.rmtree(test_merged)
  merge_j2_directories(test_primary, test_override, test_merged)
