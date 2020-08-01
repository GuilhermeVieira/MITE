import pytest
from pathlib import Path
from paramwriter import ParamWriter


def test_should_create_param_file_without_custom_ms_settings(tmp_path):
    pw = ParamWriter(tmp_path)
    project_name = 'TestParamWriter'
    mzxml_directory = Path('/foo/bar')
    param_file_path = pw.create_param_file(project_name, mzxml_directory)
    expected_file_path = str(tmp_path.joinpath('param.def'))

    assert expected_file_path == str(param_file_path), f'param.def is {str(param_file_path)}, but should be {expected_file_path} '
    assert param_file_path.exists(), f'param.def was not properly created'
    assert param_file_path.is_file(), 'param.def was created but it is not a file'
    assert param_file_path.read_text() == f'MY PROJECT NAME={project_name}\nMZXML DIRECTORY=/foo/bar\nROOT PARAMETER FILE=./ROOT_PARAM.def\n', 'there were wrong values written'


def test_should_create_param_file_with_custom_ms_settings(tmp_path):
    custom_ms_settings = {
        "a": 2.0,
        "b": 1.0
    }
    pw = ParamWriter(tmp_path, custom_ms_settings)
    project_name = 'TestParamWriter'
    mzxml_directory = Path('/foo/bar')
    param_file_path = pw.create_param_file(project_name, mzxml_directory)
    expected_file_path = str(tmp_path.joinpath('param.def'))

    assert expected_file_path == str(param_file_path), f'param.def is {str(param_file_path)}, but should be {expected_file_path} '
    assert param_file_path.exists(), f'param.def was not properly created'
    assert param_file_path.read_text() == f'MY PROJECT NAME={project_name}\nMZXML DIRECTORY=/foo/bar\nROOT PARAMETER FILE=./ROOT_PARAM.def\na=2.0\nb=1.0\n'
    assert param_file_path.is_file(), 'param.def was created but it is not a file'


def test_should_not_create_param_path_in_invalid_location():
    invalid_path: Path = Path('invalid/path')
    pw = ParamWriter(invalid_path)
    project_name = 'TestParamWriter'
    mzxml_directory = Path('/foo/bar')

    with pytest.raises(FileNotFoundError):
        pw.create_param_file(project_name, mzxml_directory)


def test_should_not_create_param_file_with_invalid_superhirn_location():
    invalid_superhirn_location = Path('invalid/path')
    pw = ParamWriter(invalid_superhirn_location)
    project_name = 'TestParamWriter'
    mzxml_directory = Path('/foo/bar')
    with pytest.raises(FileNotFoundError):
        pw.create_param_file(project_name, mzxml_directory)


def test_should_delete_param(tmp_path):
    tmp_path.joinpath('param.def').touch()
    pw = ParamWriter(tmp_path)
    pw.delete_param_file()

    assert not tmp_path.joinpath('param.def').exists(), 'the param.def file was not properly deleted'


def test_delete_param_files_should_ignore_if_param_file_does_not_exist(tmp_path):
    pw = ParamWriter(tmp_path)
    pw.delete_param_file()

    assert not tmp_path.joinpath('param.def').exists()


def test_should_create_param_file_and_delete_as_context_manager(tmp_path):
    with ParamWriter(tmp_path) as pw:
        param_file_path = pw.create_param_file('TestParamWriter', Path('foo/bar'))
        assert param_file_path.is_file()
        assert str(param_file_path) == str(tmp_path.joinpath('param.def'))

    assert not param_file_path.exists()




