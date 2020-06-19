def re_match(target, expr):
  '''
  Usage:
    ```j2
    {% set a, b = re_match("hello world", "^(\w+) (\w+)$") %}
    print("{{ a }} {{ b }}")
    ```
  '''
  import re
  return re.match(expr, str(target)).groups()
