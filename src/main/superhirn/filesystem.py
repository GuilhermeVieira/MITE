import os
import shutil
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


def change_directory(directory: Path):
    if not directory.is_dir():
        raise FileNotFoundError(f'{str(directory)} is not a valid directory')
    os.chdir(str(directory))


def copy_file(source_file: Path, destination_file: Path, create_parents=False):
    if not source_file.exists():
        raise FileNotFoundError(f'{str(source_file)} does not exist')

    destination_directory = destination_file.parent
    if create_parents and not destination_directory.is_dir():
        destination_directory.mkdir(parents=True)

    if not destination_directory.is_dir():
        raise FileNotFoundError(f'{str(destination_directory)} is not a valid destination directory')

    shutil.copyfile(source_file, destination_file)


def remove_all_folders_starting_with(prefix: str, root_directory: Path):
    if not root_directory.is_dir():
        raise FileNotFoundError(f'{str(root_directory)} is not a valid directory')

    for item in root_directory.iterdir():
        if item.is_dir() and str(item.stem).startswith(prefix):
            shutil.rmtree(item)


if __name__ == '__main__':
    path = '../../../input/mzXML/'


