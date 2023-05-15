import pandas as pd
import openpyxl
class ExcelRW:
    accuracyFile = ''
    systemUsageFile = ''
    def writeGRch37VarDataInExcel(self, grch37vardata_dict, fname):
        df = pd.DataFrame.from_dict(grch37vardata_dict, orient ='index')
        df.to_excel(fname)
    def writeGRch38VarDataInExcel(self, grch38vardata_dict, fname):
        #df = pd.DataFrame(list(grch38vardata_dict.items()), columns=['rsId', 'Chromosome Coordinates'])
        df = pd.DataFrame.from_dict(grch38vardata_dict, orient='index')
        df.to_excel(fname)
