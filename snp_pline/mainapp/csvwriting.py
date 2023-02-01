import csv
class CSVFileWriting:
    accuracyFile=''
    systemUsageFile = ''
    def writeGRch38VarDataCSV(self, chr_coord_list, fname):
        fieldnames = ['rsId','HGVS Id', 'Chromosome_Coordinates']
        print('data', chr_coord_list)
        count = 0
        with open(fname, 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(fieldnames)
            writer.writerows(chr_coord_list)

    def writeGRch37VarDataCSV(self, chr_coord_list, fname):
        fieldnames = ['rsId', 'HGVS Id','Chromosome_Coordinates']
        print('data', chr_coord_list)
        count = 0
        with open(fname, 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(fieldnames)
            writer.writerows(chr_coord_list)

    def writeAnnotationDateCSV(self, AnnoData, fname,fieldnames):
        #print('fname', fname)
        with open(fname, 'w', encoding='utf-8',newline='') as f:
            writer = csv.writer(f)
            writer.writerow(fieldnames)
            writer.writerows(AnnoData)



