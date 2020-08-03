from pathlib import Path
from typing import Dict, Optional


class ParamWriter:

    def __init__(self, superhirn_path: Path, custom_ms_settings: Optional[Dict[str, float]] = None):
        self.__superhirn_path: Path = superhirn_path
        self.__custom_ms_settings: Optional[Dict[str, float]] = custom_ms_settings

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.delete_param_file()
    
    def create_param_file(self, project_name: str, mzxml_directory: Path) -> Path:
        param_file_path: Path = self.__superhirn_path.joinpath('param.def')
        with param_file_path.open('w') as f:
            f.write(f'MY PROJECT NAME={project_name}\n')
            f.write(f'MZXML DIRECTORY={str(mzxml_directory)}\n')
            f.write('ROOT PARAMETER FILE=./ROOT_PARAM.def\n')

            if self.__custom_ms_settings:
                for key, value in self.__custom_ms_settings.items():
                    f.write(f'{key}={value}\n')

        return param_file_path

    def delete_param_file(self):
        param_file_path: Path = self.__superhirn_path.joinpath('param.def')
        if param_file_path.is_file():
            param_file_path.unlink()
