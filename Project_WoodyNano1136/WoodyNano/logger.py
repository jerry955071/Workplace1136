import os

class logger(object):
    def __init__(self, path:str) -> None:
        self.__path = path

    def new(self, times=None):
        with open(self.__path, 'w'):
            os.utime(self.__path, times)
        return

    def append(self, lines=None):
        with open(self.__path, 'a') as o:
            o.write(lines+'\n')
        return
