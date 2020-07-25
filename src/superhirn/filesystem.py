import os
import pathlib

def find_subdirs(dirpath):
    return [f'{dirpath}/{dir}' for dir in os.listdir(dirpath)]

def contains_format(dir, format):
        for file in os.listdir(dir):
            if pathlib.Path(file).suffix == format.value:
                return True

        return False

