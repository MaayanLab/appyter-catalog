import json
from pathlib import Path

def load_static(path):
    with (Path(__file__).parent.parent / 'static' / path).open('r') as fr:
        return json.load(fr)
