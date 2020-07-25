import pytest

from superhirn.settingsloader import *
import json


def test_should_load_json_file(tmp_path):
    settings = {"a": 2.0, "b": 3.0}
    json_file_path = tmp_path.joinpath('custom_settings.json')
    with json_file_path.open('w') as json_file:
        json.dump(settings, json_file)

    assert load_superhirn_ms_settings(json_file_path) == settings


def test_should_not_load_if_file_does_not_exist():
    invalid_path = Path('invalid/path')
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        load_superhirn_ms_settings(invalid_path)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == f'ERROR: Could not load mass spectrometer custom settings. ' \
                                          f'File {invalid_path} could not be located.'


def test_should_not_load_non_json_file(tmp_path):
    non_json_file_path = tmp_path.joinpath('custom_settings.ini')
    non_json_file_path.touch()
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        load_superhirn_ms_settings(non_json_file_path)
    assert pytest_wrapped_e.type == SystemExit
    assert pytest_wrapped_e.value.code == f'ERROR: Could not load mass spectrometer custom settings. ' \
                                          f'{str(non_json_file_path)} is not a json file.'
