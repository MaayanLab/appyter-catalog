''' pronto has proven very inconsistent for our needs,
so here's our own ontology parsers.

All parsers construct a dictionary of the form:

{
  'ONTOID:id': {
    id: 'ONTOID:id',
    synonyms: [synonym1, synonym2, ...],
    ... # ontology specific
  },
  ...
}
'''

import re
import xml.etree.ElementTree as ET

def ensure_list(L):
  if type(L) == list: return L
  else: return [L]

class Ontology(dict):
  def reversed(self):
    __reversed = getattr(self, '__reversed', None)
    if __reversed is None:
      __reversed = {
        node['name']: node_id
        for node_id, node in self.items()
      }
      setattr(self, '__reversed', __reversed)
    return __reversed

  def reversed_synonyms(self):
    __reversed_synonyms = getattr(self, '__reversed_synonyms', None)
    if __reversed_synonyms is None:
      __reversed_synonyms = {
        name: node_id
        for node_id, node in self.items()
        for name in [node['name']] + node['synonyms']
      }
      setattr(self, '__reversed_synonyms', __reversed_synonyms)
    return __reversed_synonyms


class OBOOntology(Ontology):
  @classmethod
  def parse(cls, source):
    with open(source, 'r') as fr:
      return cls({
        block['id']: dict(block, synonyms=set(ensure_list(block.get('synonym', []))))
        for block in cls._walk_obo(fr)
        if block['_type'] == 'Term'
      })

  _section_re = re.compile(r'^\[([^\]]+)\]$')
  _kv_re = re.compile(r'^([^:]+):\s*(.+)$')

  @staticmethod
  def _walk_obo(fr):
    _type = None
    buf = []
    blocks = []
    for line in map(str.strip, fr):
      if not line: continue
      if line.startswith('!'): continue
      m = OBOOntology._section_re.match(line)
      if m is None:
        try:
          k, v = OBOOntology._kv_re.match(line).groups()
          buf.append((k, v))
        except:
          print(line)
          raise
        continue
      else:
        yield OBOOntology._prepare_block(_type, buf)
        buf = []
        _type = m.group(1)
    if buf: yield OBOOntology._prepare_block(_type, buf)

  @staticmethod
  def _prepare_block(_type, block_buf):
    # we've buffered a bunch of parsed lines and hit a new section
    # create a jsonld-style dictionary, add the section type, and add the block
    as_dict = {}
    for k, v in block_buf:
      if k not in as_dict: as_dict[k] = [v]
      else: as_dict[k].append(v)
    for k in list(as_dict.keys()):
      if len(as_dict[k]) == 1: as_dict[k] = as_dict[k][0]
    return dict(as_dict, _type=_type)

class CellosaurusOntology(Ontology):
  @classmethod
  def parse(cls, source):
    root = ET.parse(source).getroot()
    return cls({
      accession.text.replace('_', ':'): {
        'id': accession.text.replace('_', ':'),
        'name': cell_line.find('name-list').find("name[@type='identifier']").text,
        'synonyms': {
          synonym.text
          for synonym in cell_line.find('name-list').iterfind("name[@type='synonym']")
        }
      }
      for cell_line in root.find('cell-line-list').iterfind('cell-line')
      for accession in cell_line.find('accession-list').iterfind('accession')
    })
