import os

class logger(object):
    def __init__(self, path=None):
        self._path = path

    def new(self, times=None):
        with open(self._path, 'w'):
            os.utime(self._path, times)
        return

    def append(self, lines=None):
        with open(self._path, 'a') as o:
            o.write(lines+'\n')
        return
