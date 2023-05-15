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
        print('fname', fname)
        count = 0
        with open(fname, 'w', encoding='UTF8', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(fieldnames)
            writer.writerows(chr_coord_list)

    def writeAnnotationDateCSV(self, AnnoData, fname,fieldnames):
        print('fname sift=== ', fname)
        print('AnnoData=== ', AnnoData)
        with open(fname, 'w', encoding='utf-8',newline='') as f:
            writer = csv.writer(f)
            writer.writerow(fieldnames)
            writer.writerows(AnnoData)

    def readCSVFile(self, fname):
        rowData=[]
        allData=[]

        with open(fname) as csv_file:

            csv_reader = csv.reader(csv_file, delimiter=',')
            line_count = 0
            for row in csv_reader:

                if line_count == 0:
                   line_count += 1
                else:
                   # print('row...', row)
                    allData.append(row)
                    line_count += 1
            return allData;




