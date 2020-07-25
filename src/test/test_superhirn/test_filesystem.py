from pathlib import Path
import superhirn.filesystem as fs
from superhirn.filesystem import Format


def test_contains_format_should_find_mzXML():
    assert fs.contains_format(Path('../resources/format-supported-files'), Format.MZXML) is True


def test_contains_format_should_find_raw():
    assert fs.contains_format(Path('../resources/format-supported-files'), Format.RAW) is True


def test_contains_format_with_invalid_path():
    assert fs.contains_format(Path("invalid path"), Format.MZXML) is False


def test_contains_format_with_file_instead_of_directory():
    assert fs.contains_format(Path('../test_superhirn/test_filesystem.py'), Format.RAW) is False
