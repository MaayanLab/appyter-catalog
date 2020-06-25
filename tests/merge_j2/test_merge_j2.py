def test_merge_j2():
  from compose.merge_j2 import merge_j2
  test_0 = open(os.path.join(os.path.dirname(__file__), 'test_merge_j2_0.j2'), 'r').read()
  test_1 = open(os.path.join(os.path.dirname(__file__), 'test_merge_j2_1.j2'), 'r').read()
  expectation = open(os.path.join(os.path.dirname(__file__), 'test_merge_j2_01.j2'), 'r').read()
  result = merge_j2(test_0, test_1)
  print(result)
  assert expectation == result
