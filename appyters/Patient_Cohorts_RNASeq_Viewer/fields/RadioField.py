from appyter.fields import Field

class RadioField(Field):
    def __init__(self, **kwargs):
    # if we wanted to inject something/edit the math here, we could do it this way:
    # choices = [ {**choice, "math": "Hello"} for choice in kwargs["choices"] ]
    # kwargs = {**kwargs, "choices": choices}
        super().__init__(**kwargs)

    def constraint(self):
        return self.raw_value is not None

    def to_cwl(self):
        _choices = self.args['choices']
        self.args['choices'] = [choice['value'] for choice in _choices]
        cwl = super().to_cwl()
        self.args['choices'] = _choices
        return cwl
