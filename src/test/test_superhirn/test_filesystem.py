from pathlib import Path
import os
import pytest
import filesystem as filesystem
from filesystem import Format


# Test contains_format()
def test_contains_format_should_find_mzxml(tmp_path):
    tmp_path.joinpath('file.mzXML').touch()
    assert filesystem.contains_format(tmp_path, Format.MZXML) is True


def test_contains_format_should_find_raw(tmp_path):
    tmp_path.joinpath('file.raw').touch()
    assert filesystem.contains_format(tmp_path, Format.RAW) is True


def test_contains_format_with_invalid_path():
    assert filesystem.contains_format(Path("invalid path"), Format.MZXML) is False


def test_contains_format_with_file_instead_of_directory():
    assert filesystem.contains_format(Path('../test_superhirn/test_filesystem.py'), Format.RAW) is False


# Test change_directory()
def test_should_change_directory(tmp_path):
    filesystem.change_directory(tmp_path)
    assert os.getcwd() == str(tmp_path)


def test_should_raise_file_not_found_exception():
    with pytest.raises(FileNotFoundError):
        path = Path('invalid/path')
        filesystem.change_directory(path)


# Test copy_file()
def test_should_copy_file(tmp_path):
    file_name = 'file.txt'
    source_directory = tmp_path.joinpath('source')
    destination_directory = tmp_path.joinpath('destination')
    source_file = source_directory.joinpath(file_name)
    destination_file = destination_directory.joinpath(file_name)
    source_directory.mkdir()
    destination_directory.mkdir()
    source_file.touch()
    filesystem.copy_file(source_file, destination_file)
    assert source_file.exists() and destination_file.exists()


def test_should_not_copy_file_if_it_does_not_exist(tmp_path):
    with pytest.raises(FileNotFoundError):
        file_name = 'file.txt'
        source_file = tmp_path.joinpath(f'source/{file_name}')
        destination_file = tmp_path.joinpath(f'destination/{file_name}')
        filesystem.copy_file(source_file, destination_file)


def test_should_not_copy_file_if_destination_parent_folder_does_not_exist(tmp_path):
    with pytest.raises(FileNotFoundError):
        file_name = 'file.txt'
        source_directory = tmp_path.joinpath('source')
        destination_directory = tmp_path.joinpath('destination')
        source_file = source_directory.joinpath(file_name)
        destination_file = destination_directory.joinpath(file_name)
        source_directory.mkdir()
        source_file.touch()
        filesystem.copy_file(source_file, destination_file)


def test_should_copy_file_if_destination_parent_folder_does_not_exist_and_create_parents_is_true(tmp_path):
    file_name = 'file.txt'
    source_directory = tmp_path.joinpath('source')
    destination_directory = tmp_path.joinpath('destination')
    source_file = source_directory.joinpath(file_name)
    destination_file = destination_directory.joinpath(file_name)
    source_directory.mkdir()
    source_file.touch()
    filesystem.copy_file(source_file, destination_file, create_parents=True)
    assert source_file.exists() and destination_file.exists()


def test_should_copy_file_if_destination_parent_folders_do_not_exist_and_create_parents_is_true(tmp_path):
    file_name = 'file.txt'
    source_directory = tmp_path.joinpath('source')
    destination_directory = tmp_path.joinpath('destination/folder1/folder2')
    source_file = source_directory.joinpath(file_name)
    destination_file = destination_directory.joinpath(file_name)
    source_directory.mkdir()
    source_file.touch()
    filesystem.copy_file(source_file, destination_file, create_parents=True)
    assert source_file.exists() and destination_file.exists()


# Test remove_all_folders_starting_with()
def test_should_remove_all_folders_starting_with(tmp_path):
    folder1 = tmp_path.joinpath('ANALYSIS_1')
    folder2 = tmp_path.joinpath('ANALYSIS_2')
    folder1.mkdir()
    folder2.mkdir()
    filesystem.remove_all_folders_starting_with("ANALYSIS_", tmp_path)
    assert not folder1.exists() and not folder2.exists()



