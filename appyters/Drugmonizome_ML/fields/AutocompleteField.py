import re
from appyter.fields import Field
import json, requests

class AutocompleteField(Field):
  ''' Representing a field that accepts a string with autocomplete.
  Be careful with this field, consider defining a constraint regex.
  :param name: (str) A name that will be used to refer to the object as a variable and in the HTML form.
  :param label: (str) A human readable label for the field for the HTML form
  :param description: (Optional[str]) A long human readable description for the field for the HTML form
  :param constraint: A regular expression for validating the file name.
  :param hint: A hint to put in the field prior to content.
  :param choices: (Union[List[str], Set[str], Dict[str, str]]) A set of choices that are available for this field or lookup table mapping from choice label to resulting value
  :param default: (str) A default value as an example and for use during prototyping
  :param section: (Optional[str]) The name of a SectionField for which to nest this field under, defaults to a root SectionField
  :param value: (INTERNAL Any) The raw value of the field (from the form for instance)
  :param \**kwargs: Additional keyword arguments used by other fields
  '''
  def __init__(self, constraint=r'.*', hint=None, **kwargs):
    super().__init__(
      constraint=constraint,
      hint=hint,
      **kwargs,
    )

  def constraint(self):
    return (self.raw_value is not None and
            re.match(self.args['constraint'], self.raw_value) and
            (self.args['choices'] is None or len(self.args['choices']) == 0 or self.raw_value in self.args['choices']))