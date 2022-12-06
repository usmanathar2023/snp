import csv
import re
class CSVFileWriting:
    accuracyFile=''
    systemUsageFile = ''
    def writeCSV(self, files):
        header = ['Test File','SP-Score', 'Modeler', 'SPFN', 'SPFP','Compression(naive)','Compression', 'TC']
        count=0
        self.accuracyFile='media2/accuracy.csv'
        csvFile=open(self.accuracyFile, 'w', encoding='UTF8', newline='')
        writer = csv.writer(csvFile)
        # write the header
        writer.writerow(header)
        data2=[]
        for file in files:
            data = []
            lastSlashIndex=str(file).rindex('/')
            lastUnderScoreIndex=str(file).rindex('_')
            fName=str(file)[lastSlashIndex+1:lastUnderScoreIndex]
            data.append(fName)
            accContentFile = open(file, "r")
            while True:
                # Get next line from file
                newLine = accContentFile.readline()
                # if line is empty
                # end of file is reached
                if not newLine:
                    break

                lastSpaceIndex = newLine.rindex(' ')
                lastSubString = newLine[lastSpaceIndex:]
                data.append(lastSubString.strip())
            accContentFile.close()
            data2.append(data)
        # write multiple rows
        writer.writerows(data2)

    def writeAlignmentMemTimeUsageCSV(self, data):
        header = ['File', 'Time', 'Memory', 'CPU']
        count = 0
        self.systemUsageFile = 'sequences/systemusage.csv'
        csvFile = open(self.systemUsageFile, 'w', encoding='UTF8', newline='')
        writer = csv.writer(csvFile)
        # write the header
        writer.writerow(header)
        writer.writerows(data)