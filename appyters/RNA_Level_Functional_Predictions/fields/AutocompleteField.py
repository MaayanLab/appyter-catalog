from appyter.fields import Field 

class AutocompleteField(Field):
    """
    Representing field that allows search for thousands of possible RNA genes.
    :param name: (str) A name that will be used to refer to the object as a variable and in the HTML form.
    :param label: (str) A human readable label for the field for the HTML form
    :description: (Optional[str]) A long human readable description for the field for the HTML form
    :param default: (float) A default value as an example and for use during prototyping
    :param section: (Optional[str]) The name of a SectionField for which to nest this field under, defaults to a root SectionField
    :param value: (INTERNAL Any) The raw value of the field (from the form for instance)
    """
    def __init__(self, **kwargs):
        kwargs["choices"] = None 
        super().__init__(**kwargs)  
    
    def constraint(self):
        import re 
        return self.raw_value is not None and self.raw_value