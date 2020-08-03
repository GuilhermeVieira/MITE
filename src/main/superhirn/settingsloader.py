import sys
import json
from pathlib import Path
from typing import Dict


def load_superhirn_ms_settings(settings_file_path: Path) -> Dict[str, float]:
    try:
        with settings_file_path.open() as json_file:
            settings = json.load(json_file)
    except FileNotFoundError:
        sys.exit(f'ERROR: Could not load mass spectrometer custom settings. '
                 f'File {str(settings_file_path)} could not be located.')
    except json.JSONDecodeError:
        sys.exit(f'ERROR: Could not load mass spectrometer custom settings. '
                 f'{str(settings_file_path)} is not a json file.')

    return settings
