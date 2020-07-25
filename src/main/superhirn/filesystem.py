import os
from pathlib import Path
from enum import Enum


class Format(Enum):
    RAW = '.raw'
    MZXML = '.mzXML'


def contains_format(directory: Path, file_format: Format) -> bool:
    if not directory.is_dir():
        return False

    for file in directory.iterdir():
        if file.suffix == file_format.value:
            return True

    return False


if __name__ == '__main__':
    path = '../../../input/mzXML/'