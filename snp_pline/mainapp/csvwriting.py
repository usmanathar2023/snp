class CSVFileWriting:
    accuracyFile=''
    systemUsageFile = ''
    def writeGRch37VarDataCSV(self, data, fname):
        header = ['rsId', 'Chromosome_Coordinates']
        count = 0
        with open(fname, 'w') as f:
            for key in data.keys():
                coordinates=''.join(data[key])
                coordinates=' '+coordinates+'\n'
                f.write("%s:%s" %(key,coordinates))

    def writeGRch38VarDataCSV(self, data, fname):
        header = ['rsId', 'Chromosome_Coordinates']
        count = 0
        with open(fname, 'w') as f:
            for key in data.keys():
                coordinates = ''.join(data[key])
                coordinates = ' ' + coordinates + '\n'
                f.write("%s:%s" % (key, coordinates))