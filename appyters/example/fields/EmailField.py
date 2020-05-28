from appyter.fields import Field

class EmailField(Field):
  def __init__(self, **kwargs):
    super(EmailField, self).__init__(**kwargs)

  def constraint(self):
    import re
    return self.raw_value is not None and re.match(r'\w[\w\d_-]+@\w[\w\d_-]+\.\w[\w\d_-]+', self.raw_value)
