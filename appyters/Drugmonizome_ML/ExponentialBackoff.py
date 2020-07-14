from threading import Lock

class ExponentialBackoff(object):
    """
    Threadsafe class designed for implementing exponential backoff when making API requests across threads.
    """

    def __init__(self, value=0.5, min_value=1e-6):
        self.val = value
        self.min_val = min_value
        self.lock = Lock()

    def double(self):
        with self.lock:
            self.val *= 2
    
    def halve(self):
        with self.lock:
            if self.val > self.min_val:
                self.val /= 2

    def value(self):
        with self.lock:
            return self.val