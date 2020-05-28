def do_op(op, a, b):
  return {
      'add': lambda a, b: a + b,
      'subtract': lambda a, b: a - b,
      'multiply': lambda a, b: a * b,
      'divide': lambda a, b: a / b,
      'power': lambda a, b: a ** b,
  }[op.raw_value](a.raw_value, b.raw_value)
