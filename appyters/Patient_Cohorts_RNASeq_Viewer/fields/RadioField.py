from appyter.fields import Field

class RadioField(Field):
    def __init__(self, **kwargs):
    # if we wanted to inject something/edit the math here, we could do it this way:
    # choices = [ {**choice, "math": "Hello"} for choice in kwargs["choices"] ]
    # kwargs = {**kwargs, "choices": choices}
        super().__init__(**kwargs)

    def constraint(self):
      return self.raw_value is not None
