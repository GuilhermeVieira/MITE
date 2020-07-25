class ParamWriter:

    def __init__(self, superhirn_directory, custom_ms_settings):
        self.superhirn_directory = superhirn_directory
        self.custom_ms_settings = custom_ms_settings 
    
    def create_param_file(self, project_name, mzxml_directory):

        with open(f'{self.superhirn_directory}param.def', 'w') as f:
            f.write(f'MY PROJECT NAME = {project_name.upper()}\n')
            f.write(f'MZXML DIRECTORY = {mzxml_directory}\n')

            for key, value in self.custom_ms_settings.items():
                f.write(f'{key} = {value}\n')

    def delete_param_file(self):
        

if __name__ == '__main__':

    print("ParamWriter")

    superhirn_directory = '/apps/SuperHirn/SuperHirnv03/make/'
    dicio = {'a': 1, 'b': 2}
    pw = ParamWriter(superhirn_directory, dicio)
    pw.create_param_file('test', "sla/mzxml")