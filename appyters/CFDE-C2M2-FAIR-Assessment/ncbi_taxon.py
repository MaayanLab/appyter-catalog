''' This constructs a thread-safe NCBI Taxon client that
works in its own thread and caches results across any threads that access it.
'''
import os
import time
import traceback
import urllib.request
import xml.etree.ElementTree
import shelve
from threading import Thread, Event, Lock
from queue import Queue, Empty
from contextlib import contextmanager

@contextmanager
def acquire_with_timeout(lock):
  assert lock.acquire(timeout=10), 'Timeout'
  yield lock
  lock.release()

class NCBITaxonClient(Thread):
  def __init__(self, cachedir='.cached'):
    super().__init__()
    self._queue = Queue()
    self._data = {}
    self._data_lock = Lock()
    self._join = False
    os.makedirs(cachedir, exist_ok=True)
    self._local_data = shelve.open(os.path.join(cachedir, 'ncbitaxon.cache'))
  #
  def run(self):
    # when _join is set, we'll exit the loop when the queue is empty
    while not self._join:
      # if nothing new is added to the queue in a second, we'll check if
      #  join was called
      try:
        id = self._queue.get(timeout=1)
        if id is None:
          self._queue.task_done()
          break
        # process an item on the queue
        try:
          T = xml.etree.ElementTree.parse(
            urllib.request.urlopen(
              f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={id}"
            )
          )
          time.sleep(1)
          t = T.find('Taxon')
          ret = {
            'taxid': t.find('TaxId').text,
            'name': t.find('ScientificName').text,
            'rank': t.find('Rank').text,
          }
        except:
          traceback.print_exc()
          ret = None
        # update data and trigger event
        with acquire_with_timeout(self._data_lock):
          evt = self._data[id]
          self._local_data[id] = self._data[id] = ret
          self._local_data.sync()
          evt.set()
        # move on
        self._queue.task_done()
      except Empty:
        pass
  #
  def join(self, timeout=None):
    self._join = True
    self._queue.put(None)
    super().join(timeout=timeout)
    self._local_data.close()
  #
  def fetch(self, id):
    with acquire_with_timeout(self._data_lock):
      if id in self._data:
        value = self._data[id]
        if not isinstance(value, Event):
          return value
      elif id in self._local_data:
        value = self._data[id] = self._local_data[id]
        return value
      else:
        value = self._data[id] = Event()
        self._queue.put(id)
    #
    value.wait(timeout=10)
    return self.fetch(id)

@contextmanager
def create_ncbi_taxon_client(cachedir='.cached'):
  ncbi_taxon_client = NCBITaxonClient(cachedir=cachedir)
  ncbi_taxon_client.start()
  try:
    yield ncbi_taxon_client
  finally:
    ncbi_taxon_client.join()
