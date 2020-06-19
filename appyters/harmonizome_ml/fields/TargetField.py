from appyter.fields import Field

class TargetField(Field):
  '''
  Usage:
    ```j2
    {{ TargetField(
      choices={
        "A": StringField(
          name="A",
          default="a"
        ),
        "B": StringField(
          name="B",
          default="b"
        )
      }
    ) }}
    ```
  '''
  def __init__(self, choices={}, **kwargs):
    super(TargetField, self).__init__(choices=choices, **kwargs)
