from pathlib import Path
import filesystem as filesystem
from filesystem import Format
from typing import List, Optional, Dict
import settingsloader as loader
from paramwriter import ParamWriter
import subprocess


class Taxon:

    def __init__(self, path: Path):
        self._path: Path = path
        self._name: str = path.stem

    @property
    def path(self):
        return self._path

    @property
    def name(self):
        return self._name


class SuperHirnManager:

    def __init__(self, superhirn_path: Path, ms_settings_file_path: Optional[Path] = None):
        self.__NUM_MIN_TAXA: int = 2
        self.__detected_taxa: List[Taxon] = []
        self.__superhirn_path: Path = superhirn_path
        self.__ms_settings_file_path: Dict[str, float] = \
            loader.load_superhirn_ms_settings(ms_settings_file_path)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.clean()

    def process(self, directory_path: Path, raw: bool = False):
        file_format: Format = Format.RAW if raw else Format.MZXML

        if raw:
            raise NotImplementedError("SuperHirnManager does not support raw files")

        self.__detect_taxa(directory_path, file_format)
        self.__run_pipeline()

    def clean(self):
        filesystem.remove_all_folders_starting_with("ANALYSIS_", self.__superhirn_path)

    def __detect_taxa(self, dirpath: Path, file_format: Format):
        detected_taxa: List[Taxon] = []

        for dir in dirpath.iterdir():
            if filesystem.contains_format(dir, file_format):
                detected_taxa.append(Taxon(dir))

        if len(detected_taxa) < self.__NUM_MIN_TAXA:
            raise RuntimeError(
                f'''Invalid number of taxa. Found {len(detected_taxa)} in {dirpath}, 
                but the minimum is {self.__NUM_MIN_TAXA}'''
            )

        print("The following taxa were detected:")
        for taxon in detected_taxa:
            print(taxon.name)

        self.__detected_taxa = detected_taxa

    def __run_pipeline(self):
        for taxon in self.__detected_taxa:
            self.__do_feature_extraction(taxon.name, taxon.path)
            self.__build_alignment_tree(taxon.name, taxon.path)
            self.__multiple_lcms_alignmenmt(taxon.name, taxon.path)
            self.__mastermap_intensity_normalization(taxon.name, taxon.path)

        self.__multiple_alignment_between_runs()

    def __do_feature_extraction(self, run_name: str, mzxml_directory: Path):
        self.__call_superhirn('-FE', run_name, mzxml_directory)

    def __build_alignment_tree(self, run_name: str, mzxml_directory: Path):
        self.__call_superhirn('-BT', run_name, mzxml_directory)

    def __multiple_lcms_alignmenmt(self, run_name: str, mzxml_directory: Path):
        self.__call_superhirn('-CM', run_name, mzxml_directory)

    def __mastermap_intensity_normalization(self, run_name: str, mzxml_directory: Path):
        self.__call_superhirn('-IN', run_name, mzxml_directory)

    def __multiple_alignment_between_runs(self):
        for taxon in self.__detected_taxa:
            filesystem.copy_file(
                self.__superhirn_path.joinpath(f'ANALYSIS_{taxon.name}/NORMALIZED_{taxon.name}.xml'),
                self.__superhirn_path.joinpath(f'ANALYSIS_MITE/LC_MS_RUNS/NORMALIZED_{taxon.name}.xml'),
                create_parents=True
            )

        self.__call_superhirn('-AR', "MITE", Path())

        for taxon in self.__detected_taxa:
            filesystem.copy_file(
                self.__superhirn_path.joinpath(f'ANALYSIS_MITE/LC_MS_RUNS_ALIGNED/{taxon.name}.xml'),
                Path(f'/user/src/mite/input/xml/{taxon.name}.xml'),
                create_parents=True
            )

    def __call_superhirn(self, superhirn_param: str, run_name: str, mzxml_directory: Path):
        filesystem.change_directory(self.__superhirn_path)
        with ParamWriter(self.__superhirn_path, self.__ms_settings_file_path) as pw:
            pw.create_param_file(run_name, mzxml_directory)
            superhirn: str = str(self.__superhirn_path.joinpath("SuperHirnv03"))
            subprocess.call([superhirn, superhirn_param])


if __name__ == '__main__':
    path = Path('/user/src/mite/input/mzXML/')
    superhirn_directory = Path('/apps/SuperHirn/SuperHirnv03/make/')
    settings_file_path = Path('/user/src/mite/input/config/superhirn_params.json')

    with SuperHirnManager(superhirn_directory, settings_file_path) as manager:
        manager.process(Path(path))
