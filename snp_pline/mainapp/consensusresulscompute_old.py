from Bio import Entrez
import myvariant
from .dictionaryIO import DictionaryIO
import uuid
from .csvwriting import CSVFileWriting
import requests, sys, json
from .dbsnpvarientdataretrieval import DBSnpVarientDataRetrieval

from biothings_client import get_client
from .grch37variantdata import Grch37VariantData
class ConsensusResultsCompute:
    dataNotFound=[]
    toolAnnonotationNotFound={}
    varDataNotFound=[]
    isPathogenicBySelectedTools=0
    consensusResults=[]
    chekboxValues=[]
    isSift = False; isPolyphen2HVAR=False; isPolyphen2HDIV=False; isPrimateAI=False; isRevel=False;
    islisLS2=False; isFATHMM=False; isVest4=False; isSNPEff=False;  isCadd=False;
    pathogenicVarFile='';variantNotFoundFile='';
    resultsFile = '';
    siftPVars=0;polyphen2hvarPVars=0;polyphen2hdivPVars=0;primateaiPVars=0; revelPVars=0;
    ls2PVars=0;fathmmPVars=0; vest4PVars=0; snpeffPVars=0; caddPVars=0; damagingCount=0;
    numberOfToolsFound=0
    totalSNPs=0;
    def parseMyVarinats(self, request,fabricatedTerm):
        self.pathogenicVariants.clear()
        mv = myvariant.MyVariantInfo();
        self.chekboxValues.clear()
        mvVar = ''
        Entrez.email = "usman.athar@gmail.com"
        handle = Entrez.esearch(db="snp", term=fabricatedTerm, retmax=5)
        variantData = Entrez.read(handle)
        self.totalSNPs = variantData['Count']
        varids = variantData["IdList"]
        dictIO = DictionaryIO()
        csvfw = CSVFileWriting();
        self.isSift = False; self.isPolyphen2HVAR = False; self.isPolyphen2HDIV = False; self.isPrimateAI = False; self.isRevel = False;
        self.islisLS2 = False; self.isFATHMM = False; self.isVest4 = False; self.isSNPEff = False; self.isCadd = False;
        self.resultsFile='media/' + str(uuid.uuid4()) +'_ipsnp_results' + '.csv'

        self.siftPVars = 0; self.polyphen2hvarPVars = 0; self.polyphen2hdivPVars = 0; self.primateaiPVars = 0; self.revelPVars = 0;
        self.ls2PVars = 0; self.fathmmPVars = 0; self.vest4PVars = 0; self.snpeffPVars = 0; self.caddPVars = 0;self.damagingCount=0;
        self.numberOfToolsFound=0;
        loopCount = 0;
        print('chekboxValues=== ', self.chekboxValues)
        fabricatedFields = []; listOfTuplesSift = []; annotationDataRowSift = []; finalAnnotationDataSift = [];  annotationDataRowSift4g = [];
        finalAnnotationDataSift4g = [];  annotationDataRowProvean = []; finalAnnotationDataProvean = []; annotationDataRowMVP = []; finalAnnotationDataMVP = [];
        annotationDataRowPolyphen2hvar = []; finalAnnotationDataPolyphen2hvar = []; annotationDataRowPolyphen2hdiv = []; finalAnnotationDataPolyphen2hdiv = [];
        annotationDataRowPrimateai = [];  finalAnnotationDataPrimateai = [];  annotationDataRowRevel = []; finalAnnotationDataRevel = []; annotationDataRowMPC = [];
        finalAnnotationDataMPC = []; annotationDataRowMutpred = []; finalAnnotationDataMutpred = []; annotationDataRowMtaster = []; finalAnnotationDataMtaster = [];
        annotationDataRowMassessor = []; finalAnnotationDataMassessor = []; annotationDataRowMRNN = []; finalAnnotationDataMRNN = []; annotationDataRowMSVM = []
        finalAnnotationDataMSVM = [];  annotationDataRowMLR = []; finalAnnotationDataMLR = []; annotationDataRowMCAP = []; finalAnnotationDataMCAP = [];
        annotationDataRowLS2 = []; finalAnnotationDataLS2 = []; annotationDataRowFathmm = []; finalAnnotationDataFathmm = []; annotationDataRowFXF = [];
        finalAnnotationDataFXF = [];  annotationDataRowFMKL = []; finalAnnotationDataFMKL = [];  annotationDataRowBdeladdaf = []; finalAnnotationDataBdeladdaf = [];
        annotationDataRowBdelnoaf = []; finalAnnotationDataBdelnoaf = []; annotationDataRowVest4 = []; finalAnnotationDataVest4 = []; annotationDataRowDann = [];
        finalAnnotationDataDann = [];  annotationDataRowEigen = []; finalAnnotationDataEigen = [];  annotationDataRowEigenpc = []; finalAnnotationDataEigenpc = [];
        annotationDataRowDoegen2 = [];  finalAnnotationDataDoegen2 = [];  annotationDataRowGenocanyon = [];  finalAnnotationDataGenocanyon = []
        fabricatedFields.clear(); self.toolAnnonotationNotFound.clear(); self.varDataNotFound.clear();
        hgvs37Ids=self.getHGVS37Ids(varids)
        #print("hgvs37Ids",hgvs37Ids)
        varDataNotFound_local=[]
        consensusResults_local=[]
        for vars in hgvs37Ids: #varids: for (a, b, c) in itertools.zip_longest(num, color, value, fillvalue=-1):
            id = vars[0]
            hgvs37Id = vars[1] #'rs' + str(varid)
            self.isPathogenicBySelectedTools = 0
            varDataNotFound_local.clear()
            consensusResults_local.clear()
            #id='rs104894072'
            try:
                r = requests.get('http://myvariant.info/v1/variant/'+hgvs37Id)  # , headers={"Content-Type": "application/json"})
            except:
                print("Error: There was no data for" + id + ' and ' + hgvs37Id)
            #print('r.status_code= ', r.status_code)
            if r.status_code==200:
                varDic = r.json()
                if isinstance(varDic, list):
                    varDicStr=str(varDic)
                    index1=varDicStr.find('dbnsfp')
                    if index1 < 0:
                        varDataNotFound_local.append(id)
                        varDataNotFound_local.append(hgvs37Id)
                        self.varDataNotFound.append(varDataNotFound_local.copy())
                        continue
                    varDicStr=varDicStr[index1+9:]
                    index2 = varDicStr.index('dbsnp')

                    #print('index2= ', index2)
                    varDicStr=varDicStr[:index2-3]
                    varDic=eval(varDicStr)
                    dbnsfpVarDic = varDic
                else:
                    if varDic.get('dbnsfp') == None:
                        varDataNotFound_local.append(id)
                        varDataNotFound_local.append(hgvs37Id)
                        self.varDataNotFound.append(varDataNotFound_local.copy())
                        continue
                    dbnsfpVarDic = varDic.get('dbnsfp')
                #print('varDic= ', type(varDic))
                #print('varDic= ', varDic)
                #print('r.status_code= ', r.status_code)
                #print('dbnsfpVarDic= '+str(loopCount)+'  =',dbnsfpVarDic)

                if 'sift' in dbnsfpVarDic:
                    self.numberOfToolsFound +=1;
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('sift'), 'sift' ):
                        #print(x)
                        #if len(x)>0: #if x[0] == 'sift':
                        self.annotateBySift(x, annotationDataRowSift, finalAnnotationDataSift, id,hgvs37Id)


                if dbnsfpVarDic.__contains__('polyphen2'):
                    self.numberOfToolsFound += 1;
                    tempDic=dbnsfpVarDic.get('polyphen2')
                    #print('dbnsfpVarDic.get(hvar)==', dbnsfpVarDic.get('hvar'))
                    for x in dictIO.nested_dict_pairs_iterator(tempDic.get('hvar'), 'polyphen2hvar'): #Polyphen2hdiv
                        #print(x)
                        #if len(x)>0:# 'hvar':
                        self.annotateByPolyphen2hvar(x,annotationDataRowPolyphen2hvar,finalAnnotationDataPolyphen2hvar,id,hgvs37Id)

                if dbnsfpVarDic.__contains__('polyphen2'):
                    self.numberOfToolsFound += 1;
                    tempDic = dbnsfpVarDic.get('polyphen2')
                    #print('dbnsfpVarDic.get(hvar)==', tempDic.get('hdiv'))
                    for x in dictIO.nested_dict_pairs_iterator(tempDic.get('hdiv'),'polyphen2hdiv'): #Polyphen2hdiv
                       # print(x)
                        #if len(x)>0:#'hdiv':
                        self.annotateByPolyphen2hdiv(x,annotationDataRowPolyphen2hdiv,finalAnnotationDataPolyphen2hdiv,id,hgvs37Id)

                if 'primateai' in dbnsfpVarDic:
                   self.numberOfToolsFound += 1;
                   for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('primateai'),'primateai'):
                        #print(x)
                        #if len(x)>0:#if x[0] == 'primateai':
                        self.annotateByPrimateai(x,annotationDataRowPrimateai,finalAnnotationDataPrimateai,id,hgvs37Id)

                if 'revel' in dbnsfpVarDic:
                   self.numberOfToolsFound += 1;
                   for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('revel'),'revel'):
                        #print(x)
                        #if len(x)>0:#if x[0] == 'primateai':
                        self.annotateByRevel(x,annotationDataRowRevel,finalAnnotationDataRevel,id,hgvs37Id)

                if 'list-s2' in dbnsfpVarDic:
                    self.numberOfToolsFound += 1;
                    #fabricatedFields.append('dbnsfp.list-s2.pred, dbnsfp.list-s2.rankscore, dbnsfp.list-s2.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('list-s2'),'ls2'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByLS2(x, annotationDataRowLS2, finalAnnotationDataLS2, id,hgvs37Id)
                if 'fathmm' in dbnsfpVarDic:
                    self.numberOfToolsFound += 1;
                    #fabricatedFields.append('dbnsfp.fathmm.pred, dbnsfp.fathmm.rankscore, dbnsfp.fathmm.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('fathmm'),'fathmm'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByFathmm(x, annotationDataRowFathmm, finalAnnotationDataFathmm, id,hgvs37Id)
                    if loopCount==0:
                        self.fathmmFile='media/' + str(uuid.uuid4()) +'_fathmm' + '.csv'

                if 'vest4' in dbnsfpVarDic:
                    self.numberOfToolsFound += 1;
                    #fabricatedFields.append('dbnsfp.vest4.rankscore,dbnsfp.vest4.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('vest4'),'vest4'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByVest4(x, annotationDataRowVest4, finalAnnotationDataVest4, id,hgvs37Id)


            else:
                varDataNotFound_local.append(id)
                varDataNotFound_local.append(hgvs37Id)
                self.varDataNotFound.append(varDataNotFound_local.copy())
                continue
            if self.damagingCount>self.numberOfToolsFound:
                consensusResults_local.append(id)
                consensusResults_local.append(hgvs37Id)
                consensusResults_local.append('pathogenic')
                self.pconsensusResults.append(consensusResults_local.copy())
            else:
                consensusResults_local.append(id)
                consensusResults_local.append(hgvs37Id)
                consensusResults_local.append('benign')
                self.pconsensusResults.append(consensusResults_local.copy())
        ###ends outmost for loop

        fieldnames = ['rsid','HGVSId', 'Prediction']
        csvfw.writeAnnotationDateCSV(self.consensusResults, self.resultsFile,fieldnames)


    ###ends startParsing method

    def annotateBySift(self, x, annotationDataRowSift,finalAnnotationDataSift,id,hgvs37Id):
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDataRowSift.insert(2, 'Tolerant')

            else:
                annotationDataRowSift.insert(2, 'Damaging')
                self.siftPVars+=1
                self.damagingCount = self.damagingCount + 1

        if x[0] == 'converted_rankscore':
            annotationDataRowSift.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowSift.insert(1, x[1])
        if len(annotationDataRowSift) == 3:
            annotationDataRowSift.insert(0, id)
            annotationDataRowSift.insert(1, hgvs37Id)
            finalAnnotationDataSift.append(annotationDataRowSift.copy())
            annotationDataRowSift.clear()
   #Polyphen2hdiv
    def annotateByPolyphen2hvar(self, x, annotationDataRowPolyphen2hvar, finalAnnotationDataPolyphen2hvar, id,hgvs37Id):
        if x[0] == 'pred':
            if x[1] == 'B':
                annotationDataRowPolyphen2hvar.insert(2, 'Tolerant')

            else:
                annotationDataRowPolyphen2hvar.insert(2, 'Damaging')
                self.damagingCount = self.damagingCount + 1
                self.polyphen2hvarPVars+=1

        if x[0] == 'rankscore':
            annotationDataRowPolyphen2hvar.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowPolyphen2hvar.insert(1, x[1])
        if len(annotationDataRowPolyphen2hvar) == 3:
            annotationDataRowPolyphen2hvar.insert(0, id)
            annotationDataRowPolyphen2hvar.insert(1, hgvs37Id)
            finalAnnotationDataPolyphen2hvar.append(annotationDataRowPolyphen2hvar.copy())
            annotationDataRowPolyphen2hvar.clear()

    def annotateByPolyphen2hdiv(self, x, annotationDataRowPolyphen2hdiv, finalAnnotationDataPolyphen2hdiv, id,hgvs37Id):
        if x[0] == 'pred':
            if x[1] == 'B':
                annotationDataRowPolyphen2hdiv.insert(2, 'Tolerant')

            else:
                annotationDataRowPolyphen2hdiv.insert(2, 'Damaging')
                self.damagingCount = self.damagingCount + 1
                self.polyphen2hdivPVars+=1

        if x[0] == 'rankscore':
            annotationDataRowPolyphen2hdiv.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowPolyphen2hdiv.insert(1, x[1])
        if len(annotationDataRowPolyphen2hdiv) == 3:
            annotationDataRowPolyphen2hdiv.insert(0, id)
            annotationDataRowPolyphen2hdiv.insert(1, hgvs37Id)
            finalAnnotationDataPolyphen2hdiv.append(annotationDataRowPolyphen2hdiv.copy())
            annotationDataRowPolyphen2hdiv.clear()
    def annotateByPrimateai(self, x, annotationDataRowPrimateai, finalAnnotationDataPrimateai, id,hgvs37Id):
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDataRowPrimateai.insert(2, 'Tolerant')

            else:
                annotationDataRowPrimateai.insert(2, 'Damaging')
                self.damagingCount = self.damagingCount + 1
                self.primateaiPVars+=1

        if x[0] == 'rankscore':
            annotationDataRowPrimateai.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowPrimateai.insert(1, x[1])
        if len(annotationDataRowPrimateai) == 3:
            annotationDataRowPrimateai.insert(0, id)
            annotationDataRowPrimateai.insert(1, hgvs37Id)
            finalAnnotationDataPrimateai.append(annotationDataRowPrimateai.copy())
            annotationDataRowPrimateai.clear()
#annotationDataRowRevel
    def annotateByRevel(self, x, annotationDataRowRevel, finalAnnotationDataRevel, id,hgvs37Id):
        if x[0] == 'rankscore':
            annotationDataRowRevel.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowRevel.insert(1, x[1])
            if isinstance(x[1], list):
                tempX=x[1][0]
            else:
                tempX=x[1]
            if tempX<0.5:
                annotationDataRowRevel.insert(2, 'Tolerant')

            else:
                annotationDataRowRevel.insert(2, 'Damaging')
                self.damagingCount = self.damagingCount + 1
                self.revelPVars+=1

        if len(annotationDataRowRevel) == 3:
            annotationDataRowRevel.insert(0, id)
            annotationDataRowRevel.insert(1, hgvs37Id)
            finalAnnotationDataRevel.append(annotationDataRowRevel.copy())
            annotationDataRowRevel.clear()

    def annotateByFathmm(self, x, annotationDatRowFathmm, finalAnnotationDataFathmm, id,hgvs37Id):
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDatRowFathmm.insert(2, 'Tolerant')

            else:
                annotationDatRowFathmm.insert(2, 'Damaging')
                self.damagingCount = self.damagingCount + 1
                self.fathmmPVars+=1

        if x[0] == 'converted_rankscore':
            annotationDatRowFathmm.insert(0, x[1])
        if x[0] == 'score':
            annotationDatRowFathmm.insert(1, x[1])
        if len(annotationDatRowFathmm) == 3:
            annotationDatRowFathmm.insert(0, id)
            annotationDatRowFathmm.insert(1, hgvs37Id)
            finalAnnotationDataFathmm.append(annotationDatRowFathmm.copy())
            annotationDatRowFathmm.clear()

    def annotateByVest4(self, x, annotationDataRowVest4, finalAnnotationDataVest4, id,hgvs37Id):
        if x[0] == 'rankscore':
            annotationDataRowVest4.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowVest4.insert(1, x[1])
            if isinstance(x[1], list):
                tempX = x[1][0]
            else:
                tempX = x[1]
            if tempX<0.5:
                annotationDataRowVest4.insert(2, 'Tolerant')

            else:
                annotationDataRowVest4.insert(2, 'Damaging')
                self.damagingCount = self.damagingCount + 1
                self.vest4PVars+=1

        if len(annotationDataRowVest4) == 3:
            annotationDataRowVest4.insert(0, id)
            annotationDataRowVest4.insert(1, hgvs37Id)
            finalAnnotationDataVest4.append(annotationDataRowVest4.copy())
            annotationDataRowVest4.clear()
   
    def getHGVS37Ids(self, varids):
        dbsnpvardata = DBSnpVarientDataRetrieval()
        grch37vardatainstance = Grch37VariantData()
        for varid in varids:
           dbsnpvardata_str = dbsnpvardata.getvariantdata(varid)
           grch37vardatainstance.parsevardatabystring(dbsnpvardata_str, varid)
        hgvs37Ids=grch37vardatainstance.hgvs37Ids
        return hgvs37Ids