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
    caddDataNotFound=[]
    snpeffDataNotFound=[]
    isPathogenicBySelectedTools=0 #this variable is for computing a variant as pathogenic by all selected tools i.e. finding common pathogenic
    pathogenicVariants=[]
    chekboxValues=[]
    isSift = False; isSift4G=False; isProvean=False; isMVP=False; isPolyphen2HVAR=False; isPolyphen2HDIV=False; isPrimateAI=False; isRevel=False;
    isMPC=False; isMutPred=False; isMTaster=False; isMAssessor=False; isMRNN=False; isMSVM=False; isMLR=False; isMCAP=False; isLS2=False; isFATHMM=False;
    isFXF=False; isFMKL=False; isBDelAddAF=False; isBDelNoAF=False; isVest4=False; isDANN=False; isEigen=False; isEigenPC=False; isDeogen2=False;
    isGenoCanyon=False; isCadd=False; isSneEff=False; pathogenicVarFile='';variantNotFoundFile='';

    resultsFile = '';
    totalSNPs=0;
    damagingCount=0;
    numberOfToolsFound = 0;
    consensusResults =[];
    def parseMyVarinats(self, request,fabricatedTerm,rsids):
        self.pathogenicVariants.clear()
        mv = myvariant.MyVariantInfo();
        self.chekboxValues.clear()
        mvVar = ''
        varids = []
        if len(rsids) > 0: #and rsids != NULL:
            rsids=rsids.split(",")
            for rsid in rsids:
                varids.append(rsid.removeprefix("rs").strip())
        else:
            Entrez.email = "usman.athar@gmail.com"
            handle = Entrez.esearch(db="snp", term=fabricatedTerm, retmax=10)
            variantData = Entrez.read(handle)
            #self.totalSNPs = variantData['Count']
            varids = variantData["IdList"]
        print("varids", varids)
        dictIO = DictionaryIO()
        csvfw = CSVFileWriting();
        self.isSift = False; self.isSift4G = False; self.isProvean = False; self.isMVP = False; self.isPolyphen2HVAR = False; self.isPolyphen2HDIV = False;
        self.isPrimateAI = False; self.isRevel = False; self.isMPC = False; self.isMutPred = False; self.isMTaster = False; self.isMAssessor = False;
        self.isMRNN = False; self.isMSVM = False; self.isMLR = False; self.isMCAP = False; self.isLS2 = False; self.isFATHMM = False; self.isFXF = False;
        self.isFMKL = False; self.isBDelAddAF = False; self.isBDelNoAF = False; self.isVest4 = False; self.isDANN = False; self.isEigen = False; self.isEigenPC = False;
        self.isDeogen2 = False;  self.isGenoCanyon = False; self.isCadd=False; self.isSneEff=False;
        self.chekboxValues = ((request.POST.getlist('chekboxValues', '')))
        print('self.chekboxValues === ', self.chekboxValues)
        if self.chekboxValues.__contains__('on'):
            self.numberOfToolsFound-=1

        self.numberOfToolsFound += len(self.chekboxValues);
        print('Selected Tools === ', self.numberOfToolsFound)

        loopCount = 0;

        fabricatedFields = []; listOfTuplesSift = [];
        annotationDataRowSift = []; finalAnnotationDataSift = [];  annotationDataRowSift4g = [];
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
        self.totalSNPs = len(hgvs37Ids)
        print('self.totalSNPs inside consus tools==', self.totalSNPs)
        #print("hgvs37Ids",hgvs37Ids)
        varDataNotFound_local=[]
        consensusResults_local=[]
        self.consensusResults=[]
        caddStr='';
        snpeffStr='';
        caddSourceFound = 0; dbnsfpSourceFound = 0; snpeffSourceFound = 0;
        self.resultsFile = 'media/' + str(uuid.uuid4()) + '_ipsnp_results' + '.csv'
        #print('hgvs37Ids== ',hgvs37Ids)
        #print('varids== ', varids)
        for vars in hgvs37Ids: #varids: for (a, b, c) in itertools.zip_longest(num, color, value, fillvalue=-1):
            id = vars[0]
            hgvs37Id = vars[1] #'rs' + str(varid)
            self.isPathogenicBySelectedTools = 0
            varDataNotFound_local.clear()
            consensusResults_local.clear()
            caddDataNotFound_local=[]
            snpeffDataNotFound_local=[]
            self.damagingCount=0;
            caddSourceFound = 0; dbnsfpSourceFound = 0; snpeffSourceFound = 0;
            print('hgvs37Id== ',hgvs37Id)
            print('id== ', id)
            try:
                r = requests.get('http://myvariant.info/v1/variant/'+hgvs37Id)  # , headers={"Content-Type": "application/json"})
            except:
                print("Error: There was no data for" + id + ' and ' + hgvs37Id)
            #print('r.status_code= ', r.status_code)
            if r.status_code==200:
                varDic = r.json()
                if isinstance(varDic, list):
                    varDicStr=str(varDic)
                    ##extracting cadd string
                    index1 = varDicStr.find('cadd')

                    if index1 < 0:
                        caddDataNotFound_local.append(id)
                        caddDataNotFound_local.append(hgvs37Id)
                        self.caddDataNotFound.append(caddDataNotFound_local.copy())
                        caddSourceFound=0
                    else:

                        caddStr=varDicStr[index1:]
                        index2 = caddStr.index('dbnsfp')
                        #print('index2= ', index2)
                        caddStr=caddStr[:index2-3]
                        caddSourceFound = 1
                    #varDic=eval(varDicStr)
                   # extracting dbnsfp dictionary
                    index1=varDicStr.find('dbnsfp')
                    if index1 < 0:
                        varDataNotFound_local.append(id)
                        varDataNotFound_local.append(hgvs37Id)
                        self.varDataNotFound.append(varDataNotFound_local.copy())
                        dbnsfpSourceFound=0
                    else:
                        varDicStr=varDicStr[index1+9:]
                        index2 = varDicStr.index('dbsnp')
                        #print('index2= ', index2)
                        varDicStr=varDicStr[:index2-3]
                        varDic=eval(varDicStr)
                        dbnsfpVarDic = varDic
                        dbnsfpSourceFound = 1
                    # extracting snpeff string
                    index1 = varDicStr.find('snpeff')
                    if index1 < 0:
                        snpeffDataNotFound_local.append(id)
                        snpeffDataNotFound_local.append(hgvs37Id)
                        self.snpeffDataNotFound.append(snpeffDataNotFound_local.copy())
                        snpeffSourceFound = 0
                    else:
                        varDicStr = varDicStr[index1 + 10:]
                        index2 = varDicStr.index('vcf-1')
                        # print('index2= ', index2)
                        snpeffStr = varDicStr[:index2 - 3]
                        snpeffSourceFound = 1
                else:
                    if varDic.get('cadd') == None:
                        caddDataNotFound_local.append(id)
                        caddDataNotFound_local.append(hgvs37Id)
                        self.caddDataNotFound.append(caddDataNotFound_local.copy())
                        caddSourceFound = 0
                    else:
                        caddStr= str(varDic.get('cadd'))
                        caddSourceFound = 1
                        #print('caddStr= ', caddStr)
                    if varDic.get('dbnsfp') == None:
                        varDataNotFound_local.append(id)
                        varDataNotFound_local.append(hgvs37Id)
                        self.varDataNotFound.append(varDataNotFound_local.copy())
                        dbnsfpSourceFound = 0
                    else:
                        dbnsfpVarDic = varDic.get('dbnsfp')
                        dbnsfpSourceFound = 1
                    if varDic.get('snpeff') == None:
                        snpeffDataNotFound_local.append(id)
                        snpeffDataNotFound_local.append(hgvs37Id)
                        self.snpeffDataNotFound.append(snpeffDataNotFound_local.copy())
                        snpeffSourceFound = 0
                    else:
                        snpeffStr = str(varDic.get('snpeff'))
                        snpeffSourceFound = 1 #
                if caddSourceFound==0 and dbnsfpSourceFound==0 and snpeffSourceFound==0:
                    continue

                #print('caddStr= ', caddStr)
                #print('dbnsfpVarDic= ', dbnsfpVarDic)
                #print('snpeffStr= ', snpeffStr)
                #print('dbnsfpVarDic= '+str(loopCount)+'  =',dbnsfpVarDic)

                if self.chekboxValues.__contains__('sft') and  'sift' in dbnsfpVarDic:
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('sift'), 'sift' ):
                        #print(x)
                        #if len(x)>0: #if x[0] == 'sift':
                        self.annotateBySift(x, annotationDataRowSift, finalAnnotationDataSift, id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['sift' + str(loopCount)] = id

                if self.chekboxValues.__contains__('sift4g') and  'sift4g' in dbnsfpVarDic:
                   for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('sift4g'), 'sift4g'):
                        #print(x)
                        #if len(x)>0:#if x[0] == 'sift4g':
                        self.annotateBySift4g(x,annotationDataRowSift4g,finalAnnotationDataSift4g,id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['sift4g' + str(loopCount)] = id

                if self.chekboxValues.__contains__('provean') and  'provean' in dbnsfpVarDic:
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('provean'), 'provean'):
                        #print(x)
                        #if len(x)>0:#if x[0] == 'provean':
                            self.annotateByProvean(x,annotationDataRowProvean,finalAnnotationDataProvean,id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['provean' + str(loopCount)] = id

                if self.chekboxValues.__contains__('mvp') and  'mvp' in dbnsfpVarDic:
                   for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('mvp'), 'mvp'):
                        #print(x)
                        #if len(x)>0:#if x[0] == 'mvp':
                        self.annotateByMVP(x,annotationDataRowMVP,finalAnnotationDataMVP,id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['mvp' + str(loopCount)] = id

                if self.chekboxValues.__contains__('polyphen2hvar') and dbnsfpVarDic.__contains__('polyphen2'):
                    tempDic=dbnsfpVarDic.get('polyphen2')
                    #print('dbnsfpVarDic.get(hvar)==', dbnsfpVarDic.get('hvar'))
                    for x in dictIO.nested_dict_pairs_iterator(tempDic.get('hvar'), 'polyphen2hvar'): #Polyphen2hdiv
                        #print(x)
                        #if len(x)>0:# 'hvar':
                        self.annotateByPolyphen2hvar(x,annotationDataRowPolyphen2hvar,finalAnnotationDataPolyphen2hvar,id,hgvs37Id)



                else:
                    self.toolAnnonotationNotFound['polyphen2hvar' + str(loopCount)] = id

                if self.chekboxValues.__contains__('polyphen2hdiv') and dbnsfpVarDic.__contains__('polyphen2'):
                    tempDic = dbnsfpVarDic.get('polyphen2')
                    #print('dbnsfpVarDic.get(hvar)==', tempDic.get('hdiv'))
                    for x in dictIO.nested_dict_pairs_iterator(tempDic.get('hdiv'),'polyphen2hdiv'): #Polyphen2hdiv
                       # print(x)
                        #if len(x)>0:#'hdiv':
                        self.annotateByPolyphen2hdiv(x,annotationDataRowPolyphen2hdiv,finalAnnotationDataPolyphen2hdiv,id,hgvs37Id)

                else:
                     self.toolAnnonotationNotFound['polyphen2hdiv' + str(loopCount)] = id

                if self.chekboxValues.__contains__('primateai') and  'primateai' in dbnsfpVarDic:
                   for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('primateai'),'primateai'):
                        #print(x)
                        #if len(x)>0:#if x[0] == 'primateai':
                        self.annotateByPrimateai(x,annotationDataRowPrimateai,finalAnnotationDataPrimateai,id,hgvs37Id)

                else:
                    self.toolAnnonotationNotFound['primateai' + str(loopCount)] = id

                if self.chekboxValues.__contains__('revel') and  'revel' in dbnsfpVarDic:
                   for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('revel'),'revel'):
                        #print(x)
                        #if len(x)>0:#if x[0] == 'primateai':
                        self.annotateByRevel(x,annotationDataRowRevel,finalAnnotationDataRevel,id,hgvs37Id)

                else:
                    self.toolAnnonotationNotFound['revel' + str(loopCount)] = id

                if self.chekboxValues.__contains__('mpc') and  'mpc' in dbnsfpVarDic:
                   for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('mpc'),'mpc'):
                        #print(x)
                        #if len(x)>0:#if x[0] == 'primateai':
                        self.annotateByMPC(x,annotationDataRowMPC,finalAnnotationDataMPC,id,hgvs37Id)

                else:
                    self.toolAnnonotationNotFound['mpc' + str(loopCount)] = id

                if self.chekboxValues.__contains__('mutpred') and  'mutpred' in dbnsfpVarDic:
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('mutpred'),'mutpred'):
                        #print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByMutpred(x, annotationDataRowMutpred, finalAnnotationDataMutpred, id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['mutpred' + str(loopCount)] = id

                if self.chekboxValues.__contains__('mtaster') and  'mutationtaster' in dbnsfpVarDic:

                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('mutationtaster'),'mtaster'):

                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByMtaster(x, annotationDataRowMtaster, finalAnnotationDataMtaster, id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['mtaster' + str(loopCount)] = id

                if self.chekboxValues.__contains__('massessor') and  'mutationassessor' in dbnsfpVarDic:
                   for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('mutationassessor'),'massessor'):
                        #print(x)
                        if len(x) > 0:  # if x[0] == 'primateai':
                            self.annotateByMassessor(x, annotationDataRowMassessor, finalAnnotationDataMassessor, id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['massessor' + str(loopCount)] = id

                if self.chekboxValues.__contains__('mrnn') and  'metarnn' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.metarnn.pred, dbnsfp.metarnn.rankscore, dbnsfp.metarnn.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('metarnn'),'mrnn'):
                        #print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                            self.annotateByMRNN(x, annotationDataRowMRNN, finalAnnotationDataMRNN, id,hgvs37Id)



                else:
                    self.toolAnnonotationNotFound['mrnn' + str(loopCount)] = id

                if self.chekboxValues.__contains__('msvm') and  'metasvm' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.metasvm.pred, dbnsfp.metasvm.rankscore, dbnsfp.metasvm.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('metasvm'),'msvm'):
                        #print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByMRNN(x, annotationDataRowMSVM, finalAnnotationDataMSVM, id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['msvm' + str(loopCount)] = id

                if self.chekboxValues.__contains__('mlr') and  'metalr' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.metalr.pred, dbnsfp.metalr.rankscore, dbnsfp.metalr.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('metalr'),'mlr'):
                        #print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByMLR(x, annotationDataRowMLR, finalAnnotationDataMLR, id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['mlr' + str(loopCount)] = id

                if self.chekboxValues.__contains__('mcap') and  'm-cap' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.m-cap.pred, dbnsfp.m-cap.rankscore, dbnsfp.m-cap.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('m-cap'),'mcap'):
                        #print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByMCAP(x, annotationDataRowMCAP, finalAnnotationDataMCAP, id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['mcap' + str(loopCount)] = id

                if self.chekboxValues.__contains__('ls2') and  'list-s2' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.list-s2.pred, dbnsfp.list-s2.rankscore, dbnsfp.list-s2.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('list-s2'),'ls2'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByLS2(x, annotationDataRowLS2, finalAnnotationDataLS2, id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['ls2' + str(loopCount)] = id

                if self.chekboxValues.__contains__('fathmm') and  'fathmm' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.fathmm.pred, dbnsfp.fathmm.rankscore, dbnsfp.fathmm.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('fathmm'),'fathmm'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByFathmm(x, annotationDataRowFathmm, finalAnnotationDataFathmm, id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['fathmm' + str(loopCount)] = id
                if self.chekboxValues.__contains__('fxf') and  'fathmm-xf' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.fathmm-xf.coding_pred,dbnsfp.fathmm-xf.coding_rankscore, dbnsfp.fathmm-xf.coding_score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('fathmm-xf'),'fathmm-xf'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByFXF(x, annotationDataRowFXF, finalAnnotationDataFXF, id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['fxf' + str(loopCount)] = id

                if self.chekboxValues.__contains__('fmkl') and  'fathmm-mkl' in dbnsfpVarDic:
                    #fabricatedFields.append('fathmm-mkl.coding_pred,dbnsfp.fathmm-mkl.coding_rankscore,dbnsfp.fathmm-mkl.coding_score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('fathmm-mkl'),'fathmm-mkl'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByFMKL(x, annotationDataRowFMKL, finalAnnotationDataFMKL, id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['fmkl' + str(loopCount)] = id

                if self.chekboxValues.__contains__('bdeladdaf') and  'bayesdel' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.bayesdel.add_af.pred, dbnsfp.bayesdel.add_af.rankscore,dbnsfp.bayesdel.add_af.score')
                    tempDic = dbnsfpVarDic.get('bayesdel')
                    #print('tempDic===', tempDic)
                    for x in dictIO.nested_dict_pairs_iterator(tempDic.get('add_af'),'bayesdel-addaf'):
                        #print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByBdeladdaf(x, annotationDataRowBdeladdaf, finalAnnotationDataBdeladdaf, id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['bdeladdaf' + str(loopCount)] = id

                if self.chekboxValues.__contains__('bdelnoaf') and  'bayesdel' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.bayesdel.no_af.pred, dbnsfp.bayesdel.no_af.rankscore,dbnsfp.bayesdel.no_af.score')
                    tempDic = dbnsfpVarDic.get('bayesdel')
                    #print('tempDic===', tempDic)
                    for x in dictIO.nested_dict_pairs_iterator(tempDic.get('no_af'),'bayesdel-noaf'):
                        #print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByBdelnoaf(x, annotationDataRowBdelnoaf, finalAnnotationDataBdelnoaf, id,hgvs37Id)

                else:
                    self.toolAnnonotationNotFound['bdelnoaf' + str(loopCount)] = id

                if self.chekboxValues.__contains__('vest4') and  'vest4' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.vest4.rankscore,dbnsfp.vest4.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('vest4'),'vest4'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByVest4(x, annotationDataRowVest4, finalAnnotationDataVest4, id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['vest4' + str(loopCount)] = id

                if self.chekboxValues.__contains__('dann') and  'dann' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.dann.rankscore, dbnsfp.dann.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('dann'),'dann'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByDann(x, annotationDataRowDann, finalAnnotationDataDann, id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['dann' + str(loopCount)] = id

                if self.chekboxValues.__contains__('eign') and  'eigen' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.eigen.raw_coding, dbnsfp.eigen.raw_coding_rankscore')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('eigen'),'eigen'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByEigen(x, annotationDataRowEigen, finalAnnotationDataEigen, id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['eigen' + str(loopCount)] = id

                if self.chekboxValues.__contains__('eigenpc') and  'eigen-pc' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.eigen-pc.raw_coding, dbnsfp.eigen-pc.raw_coding_rankscore')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('eigen-pc'),'eigenpc'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByEigenpc(x, annotationDataRowEigenpc, finalAnnotationDataEigenpc, id,hgvs37Id)

                else:
                    self.toolAnnonotationNotFound['eigenpc' + str(loopCount)] = id

                if self.chekboxValues.__contains__('deogen2') and  'deogen2' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.deogen2.pred,dbnsfp.deogen2.rankscore, dbnsfp.deogen2.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('deogen2'),'doegen2'):
                        # print(x)
                        if len(x) > 0:  # if x[0] == 'primateai':
                            self.annotateByDoegen2(x, annotationDataRowDoegen2, finalAnnotationDataDoegen2, id,hgvs37Id)


                else:
                    self.toolAnnonotationNotFound['deogen2' + str(loopCount)] = id

                if self.chekboxValues.__contains__('genocanyon') and  'genocanyon' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.genocanyon.rankscore,dbnsfp.genocanyon.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('genocanyon'),'genocanyon'):
                        #print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByGenocanyon(x, annotationDataRowGenocanyon, finalAnnotationDataGenocanyon, id,hgvs37Id)

                else:
                    self.toolAnnonotationNotFound['genocanyon' + str(loopCount)] = id
                if self.chekboxValues.__contains__('cadd') and caddStr.find('http://bit.ly/2TIuab9') > 0:

                    index1 = caddStr.find('phred')

                    if index1 < 0:
                        caddScoreNotFound=0
                    else:
                        caddPhredStr = caddStr[index1 + 7:]
                        #print('caddPhredStr: ', caddPhredStr)
                        index2 = caddPhredStr.find(',')
                        # print('index2= ', index2)
                        caddPgredScore = caddPhredStr[:index2]
                        caddPgredScore = caddPgredScore.strip()
                        #print('caddPgredScore: ', caddPgredScore)
                        caddPgredScoreFloat=float(caddPgredScore)
                        caddPgredScoreInt=round(caddPgredScoreFloat)
                        if caddPgredScoreInt >=20:
                            self.damagingCount += 1
                        #print('caddPgredScoreInt: ', caddPgredScoreInt)

                else:
                    self.toolAnnonotationNotFound['cadd' + str(loopCount)] = id
                if self.chekboxValues.__contains__('snpeff') and snpeffStr.find('http://bit.ly/2suyRKt')> 0:
                    index1 = snpeffStr.find('putative_impact')

                    if index1 < 0:
                        snpeeScoreNotFound = 0
                    else:
                        snpeffpredictionstr = snpeffStr[index1 + 19:]
                        #print('snpeffpredictionstr: ', snpeffpredictionstr)
                        index2 = snpeffpredictionstr.find(',')
                        # print('index2= ', index2)
                        snpeffpredition = snpeffpredictionstr[:index2-1]
                        snpeffpredition=snpeffpredition.strip()
                        if snpeffpredition=='HIGH':
                            self.damagingCount +=1
                        #print('snpeffpredition: ', snpeffpredition)


                else:
                    self.toolAnnonotationNotFound['snpeff' + str(loopCount)] = id
                loopCount +=1
            else:
                varDataNotFound_local.append(id)
                varDataNotFound_local.append(hgvs37Id)
                self.varDataNotFound.append(varDataNotFound_local.copy())
                continue
            print('self.numberOfToolsFound/2==',self.numberOfToolsFound/2)
            print('damagingCount==', self.damagingCount)
            if self.damagingCount>(self.numberOfToolsFound/2):
                consensusResults_local.append(id)
                consensusResults_local.append(hgvs37Id)
                consensusResults_local.append('pathogenic')
                consensusResults_local.append(self.damagingCount)
                self.consensusResults.append(consensusResults_local.copy())
            else:
                consensusResults_local.append(id)
                consensusResults_local.append(hgvs37Id)
                consensusResults_local.append('benign')
                consensusResults_local.append(self.numberOfToolsFound-self.damagingCount)
                self.consensusResults.append(consensusResults_local.copy())
        ###ends outmost for loop
        fieldnames = ['rsid', 'HGVSId', 'Prediction','Tools support']
        csvfw.writeAnnotationDateCSV(self.consensusResults, self.resultsFile, fieldnames)
    ###ends startParsing method

    def annotateBySift(self, x, annotationDataRowSift,finalAnnotationDataSift,id,hgvs37Id):
        x2 = ''
        if x[0] == 'pred':
            if isinstance(x[1],list):
                x2=list(x[1])
            if x[1] == 'T' or x2.__contains__('T'):

                annotationDataRowSift.insert(2, 'Tolerant')

            else:
                annotationDataRowSift.insert(2, 'Damaging')
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

    def annotateBySift4g(self, x,annotationDataRowSift4g,finalAnnotationDataSift4g, id,hgvs37Id):
        x2 = ''
        if x[0] == 'pred':
            if isinstance(x[1], list):
                x2 = list(x[1])
            if x[1] == 'T' or x2.__contains__('T'):
                annotationDataRowSift4g.insert(2, 'Tolerant')

            else:
                annotationDataRowSift4g.insert(2, 'Damaging')
                self.damagingCount = self.damagingCount + 1


        if x[0] == 'converted_rankscore':
            annotationDataRowSift4g.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowSift4g.insert(1, x[1])
        if len(annotationDataRowSift4g) == 3:
            annotationDataRowSift4g.insert(0, id)
            annotationDataRowSift4g.insert(1, hgvs37Id)
            finalAnnotationDataSift4g.append(annotationDataRowSift4g.copy())
            annotationDataRowSift4g.clear()

    def annotateByProvean(self, x, annotationDataRowProvean, finalAnnotationDataProvean, id,hgvs37Id):
        x2=''
        if x[0] == 'pred':
            print('provean x1', x[1])
            if isinstance(x[1],list):
                x2=list(x[1])
            if x[1] == 'N' or x2.__contains__('N'):
                annotationDataRowProvean.insert(2, 'Tolerant')

            else:
                annotationDataRowProvean.insert(2, 'Damaging')
                self.damagingCount = self.damagingCount + 1


        if x[0] == 'converted_rankscore':
            annotationDataRowProvean.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowProvean.insert(1, x[1])
        if len(annotationDataRowProvean) == 3:
            annotationDataRowProvean.insert(0, id)
            annotationDataRowProvean.insert(1, hgvs37Id)
            finalAnnotationDataProvean.append(annotationDataRowProvean.copy())
            annotationDataRowProvean.clear()

    def annotateByMVP(self, x, annotationDataRowMVP, finalAnnotationDataMVP, id,hgvs37Id):
        tempX=0
        #print('x[0]==', x[0])
        if x[0] == 'rankscore':
            annotationDataRowMVP.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowMVP.insert(1, x[1])
            #print('x[1]==', x[1])
            #print('isinstance(x[1], list)==', isinstance(x[1], list))
            if isinstance(x[1], list):
                tempX=x[1][0]
                #print('tempX==', tempX)
            else:
                tempX=x[1]
                #print('tempX==', tempX)
            if tempX <= 0.76:
                annotationDataRowMVP.insert(2, 'Tolerant')

            else:
                annotationDataRowMVP.insert(2, 'Damaging')
                self.damagingCount = self.damagingCount + 1


        if len(annotationDataRowMVP) == 3:
            annotationDataRowMVP.insert(0, id)
            annotationDataRowMVP.insert(1, hgvs37Id)
            finalAnnotationDataMVP.append(annotationDataRowMVP.copy())
            annotationDataRowMVP.clear()
    #Polyphen2hdiv
    def annotateByPolyphen2hvar(self, x, annotationDataRowPolyphen2hvar, finalAnnotationDataPolyphen2hvar, id,hgvs37Id):
        x2 = ''
        if x[0] == 'pred':
            if isinstance(x[1], list):
                x2 = list(x[1])
            if x[1] == 'B' or x2.__contains__('B'):
                annotationDataRowPolyphen2hvar.insert(2, 'Tolerant')

            else:
                annotationDataRowPolyphen2hvar.insert(2, 'Damaging')
                self.damagingCount = self.damagingCount + 1


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
        x2 = ''
        if x[0] == 'pred':
            if isinstance(x[1], list):
                x2 = list(x[1])
            if x[1] == 'B' or x2.__contains__('B'):
                annotationDataRowPolyphen2hdiv.insert(2, 'Tolerant')

            else:
                annotationDataRowPolyphen2hdiv.insert(2, 'Damaging')
                self.damagingCount = self.damagingCount + 1


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
        x2 = ''
        if x[0] == 'pred':
            if isinstance(x[1], list):
                x2 = list(x[1])
            if x[1] == 'T' or x2.__contains__('T'):
                annotationDataRowPrimateai.insert(2, 'Tolerant')

            else:
                annotationDataRowPrimateai.insert(2, 'Damaging')
                self.damagingCount = self.damagingCount + 1


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


        if len(annotationDataRowRevel) == 3:
            annotationDataRowRevel.insert(0, id)
            annotationDataRowRevel.insert(1, hgvs37Id)
            finalAnnotationDataRevel.append(annotationDataRowRevel.copy())
            annotationDataRowRevel.clear()

    def annotateByMPC(self, x, annotationDataRowMPC, finalAnnotationDataMPC, id,hgvs37Id):
        if x[0] == 'rankscore':
            annotationDataRowMPC.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowMPC.insert(1, x[1])
            if isinstance(x[1], list):
                tempX = x[1][0]
            else:
                tempX = x[1]
            if tempX<0.5:
                annotationDataRowMPC.insert(2, 'Tolerant')

            else:
                annotationDataRowMPC.insert(2, 'Damaging')
                self.damagingCount = self.damagingCount + 1


        if len(annotationDataRowMPC) == 3:
            annotationDataRowMPC.insert(0, id)
            annotationDataRowMPC.insert(1, hgvs37Id)
            finalAnnotationDataMPC.append(annotationDataRowMPC.copy())
            annotationDataRowMPC.clear()

    def annotateByMutpred(self, x, annotationDataRowMutpred, finalAnnotationDataMutpred, id,hgvs37Id):
        if x[0] == 'rankscore':
            annotationDataRowMutpred.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowMutpred.insert(1, x[1])
            if isinstance(x[1], list):
                tempX = x[1][0]
            else:
                tempX = x[1]
            if tempX<0.5:
                annotationDataRowMutpred.insert(2, 'Tolerant')

            else:
                annotationDataRowMutpred.insert(2, 'Damaging')
                self.damagingCount = self.damagingCount + 1


        if len(annotationDataRowMutpred) == 3:
            annotationDataRowMutpred.insert(0, id)
            annotationDataRowMutpred.insert(1, hgvs37Id)
            finalAnnotationDataMutpred.append(annotationDataRowMutpred.copy())
            annotationDataRowMutpred.clear()

    def annotateByMtaster(self, x, annotationDataRowMtaster, finalAnnotationDataMtaster, id,hgvs37Id):
        x2 = ''
        if x[0] == 'pred':
            print('mtaster x1', x[1])
            if isinstance(x[1], list):
                x2 = list(x[1])
            if x[1] == 'N' or x2.__contains__('N'):
                annotationDataRowMtaster.insert(2, 'Tolerant')

            else:
                annotationDataRowMtaster.insert(2, 'Damaging')
                self.damagingCount = self.damagingCount + 1


        if x[0] == 'converted_rankscore':
            annotationDataRowMtaster.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowMtaster.insert(1, x[1])
        if len(annotationDataRowMtaster) == 3:
            annotationDataRowMtaster.insert(0, id)
            annotationDataRowMtaster.insert(1, hgvs37Id)
            finalAnnotationDataMtaster.append(annotationDataRowMtaster.copy())
            annotationDataRowMtaster.clear()  # annotationDataRowMtaster

    def annotateByMassessor(self, x, annotationDataRowMassessor, finalAnnotationDataMassessor, id,hgvs37Id):
        x2 = ''
        if x[0] == 'pred':
            print('massessor x1', x[1])
            if isinstance(x[1], list):
                x2 = list(x[1])
            if x[1] == 'T' or x2.__contains__('T'):
                annotationDataRowMassessor.insert(2, 'Tolerant')

            else:
                annotationDataRowMassessor.insert(2, 'Damaging')
                self.damagingCount = self.damagingCount + 1


        if x[0] == 'rankscore':
            annotationDataRowMassessor.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowMassessor.insert(1, x[1])
        if len(annotationDataRowMassessor) == 3:
            annotationDataRowMassessor.insert(0, id)
            annotationDataRowMassessor.insert(1, hgvs37Id)
            finalAnnotationDataMassessor.append(annotationDataRowMassessor.copy())
            annotationDataRowMassessor.clear()

    def annotateByMRNN(self, x, annotationDatRowMRNN, finalAnnotationDataMRNN, id,hgvs37Id):
        x2 = ''
        if x[0] == 'pred':
            if isinstance(x[1], list):
                x2 = list(x[1])
            if x[1] == 'T' or x2.__contains__('T'):
                annotationDatRowMRNN.insert(2, 'Tolerant')

            else:
                annotationDatRowMRNN.insert(2, 'Damaging')

                self.damagingCount = self.damagingCount + 1
        if x[0] == 'rankscore':
            annotationDatRowMRNN.insert(0, x[1])
        if x[0] == 'score':
            annotationDatRowMRNN.insert(1, x[1])
        if len(annotationDatRowMRNN) == 3:
            annotationDatRowMRNN.insert(0, id)
            annotationDatRowMRNN.insert(1, hgvs37Id)
            finalAnnotationDataMRNN.append(annotationDatRowMRNN.copy())
            annotationDatRowMRNN.clear()

    def annotateByMSVM(self, x, annotationDatRowMSVM, finalAnnotationDataMSVM, id,hgvs37Id):
        x2 = ''
        if x[0] == 'pred':
            if isinstance(x[1], list):
                x2 = list(x[1])
            if x[1] == 'T' or x2.__contains__('T'):
                annotationDatRowMSVM.insert(2, 'Tolerant')

            else:
                annotationDatRowMSVM.insert(2, 'Damaging')

                self.damagingCount = self.damagingCount + 1
        if x[0] == 'rankscore':
            annotationDatRowMSVM.insert(0, x[1])
        if x[0] == 'score':
            annotationDatRowMSVM.insert(1, x[1])
        if len(annotationDatRowMSVM) == 3:
            annotationDatRowMSVM.insert(0, id)
            annotationDatRowMSVM.insert(1, hgvs37Id)
            finalAnnotationDataMSVM.append(annotationDatRowMSVM.copy())
            annotationDatRowMSVM.clear() #annotationDataRowMLR
    def annotateByMLR(self, x, annotationDatRowMLR, finalAnnotationDataMLR, id,hgvs37Id):
        x2 = ''
        if x[0] == 'pred':
            if isinstance(x[1], list):
                x2 = list(x[1])
            if x[1] == 'T' or x2.__contains__('T'):
                annotationDatRowMLR.insert(2, 'Tolerant')

            else:
                annotationDatRowMLR.insert(2, 'Damaging')

                self.damagingCount = self.damagingCount + 1
        if x[0] == 'rankscore':
            annotationDatRowMLR.insert(0, x[1])
        if x[0] == 'score':
            annotationDatRowMLR.insert(1, x[1])
        if len(annotationDatRowMLR) == 3:
            annotationDatRowMLR.insert(0, id)
            annotationDatRowMLR.insert(1, hgvs37Id)
            finalAnnotationDataMLR.append(annotationDatRowMLR.copy())
            annotationDatRowMLR.clear()

    def annotateByMCAP(self, x, annotationDatRowMCAP, finalAnnotationDataMCAP, id,hgvs37Id):
        x2 = ''
        if x[0] == 'pred':
            if isinstance(x[1], list):
                x2 = list(x[1])
            if x[1] == 'T' or x2.__contains__('T'):
                annotationDatRowMCAP.insert(2, 'Tolerant')

            else:
                annotationDatRowMCAP.insert(2, 'Damaging')

                self.damagingCount = self.damagingCount + 1
        if x[0] == 'rankscore':
            annotationDatRowMCAP.insert(0, x[1])
        if x[0] == 'score':
            annotationDatRowMCAP.insert(1, x[1])
        if len(annotationDatRowMCAP) == 3:
            annotationDatRowMCAP.insert(0, id)
            annotationDatRowMCAP.insert(1, hgvs37Id)
            finalAnnotationDataMCAP.append(annotationDatRowMCAP.copy())
            annotationDatRowMCAP.clear()
    def annotateByLS2(self, x, annotationDatRowLS2, finalAnnotationDataLS2, id,hgvs37Id):
        x2 = ''
        if x[0] == 'pred':
            if isinstance(x[1], list):
                x2 = list(x[1])
            if x[1] == 'T' or x2.__contains__('T'):
                annotationDatRowLS2.insert(2, 'Tolerant')

            else:
                annotationDatRowLS2.insert(2, 'Damaging')

                self.damagingCount = self.damagingCount + 1
        if x[0] == 'rankscore':
            annotationDatRowLS2.insert(0, x[1])
        if x[0] == 'score':
            annotationDatRowLS2.insert(1, x[1])
        if len(annotationDatRowLS2) == 3:
            annotationDatRowLS2.insert(0, id)
            annotationDatRowLS2.insert(1, hgvs37Id)
            finalAnnotationDataLS2.append(annotationDatRowLS2.copy())
            annotationDatRowLS2.clear()
    def annotateByFathmm(self, x, annotationDatRowFathmm, finalAnnotationDataFathmm, id,hgvs37Id):
        x2 = ''
        if x[0] == 'pred':
            print('fathmm x1', x[1])
            if isinstance(x[1], list):
                x2 = list(x[1])
            if x[1] == 'T' or x2.__contains__('T'):
                annotationDatRowFathmm.insert(2, 'Tolerant')

            else:
                annotationDatRowFathmm.insert(2, 'Damaging')

                self.damagingCount = self.damagingCount + 1
        if x[0] == 'converted_rankscore':
            annotationDatRowFathmm.insert(0, x[1])
        if x[0] == 'score':
            annotationDatRowFathmm.insert(1, x[1])
        if len(annotationDatRowFathmm) == 3:
            annotationDatRowFathmm.insert(0, id)
            annotationDatRowFathmm.insert(1, hgvs37Id)
            finalAnnotationDataFathmm.append(annotationDatRowFathmm.copy())
            annotationDatRowFathmm.clear()
    def annotateByFXF(self, x, annotationDatRowFXF, finalAnnotationDataFXF, id,hgvs37Id):
        x2 = ''
        if x[0] == 'coding_pred':
            if isinstance(x[1], list):
                x2 = list(x[1])
            if x[1] == 'T' or x2.__contains__('T'):
                annotationDatRowFXF.insert(2, 'Tolerant')

            else:
                annotationDatRowFXF.insert(2, 'Damaging')

                self.damagingCount = self.damagingCount + 1
        if x[0] == 'coding_rankscore':
            annotationDatRowFXF.insert(0, x[1])
        if x[0] == 'coding_score':
            annotationDatRowFXF.insert(1, x[1])
        if len(annotationDatRowFXF) == 3:
            annotationDatRowFXF.insert(0, id)
            annotationDatRowFXF.insert(1, hgvs37Id)
            finalAnnotationDataFXF.append(annotationDatRowFXF.copy())
            annotationDatRowFXF.clear()
    def annotateByFMKL(self, x, annotationDatRowFMKL, finalAnnotationDataFMKL, id,hgvs37Id):
        x2 = ''
        if x[0] == 'coding_pred':
            if isinstance(x[1], list):
                x2 = list(x[1])
            if x[1] == 'T' or x2.__contains__('T'):
                annotationDatRowFMKL.insert(2, 'Tolerant')

            else:
                annotationDatRowFMKL.insert(2, 'Damaging')

                self.damagingCount = self.damagingCount + 1
        if x[0] == 'coding_rankscore':
            annotationDatRowFMKL.insert(0, x[1])
        if x[0] == 'coding_score':
            annotationDatRowFMKL.insert(1, x[1])
        if len(annotationDatRowFMKL) == 3:
            annotationDatRowFMKL.insert(0, id)
            annotationDatRowFMKL.insert(1, hgvs37Id)
            finalAnnotationDataFMKL.append(annotationDatRowFMKL.copy())
            annotationDatRowFMKL.clear()
    def annotateByBdeladdaf(self, x, annotationDataRowBdeladdf, finalAnnotationDataBdeladdf, id,hgvs37Id):
        #if x[0] == 'add_af':
            x2 = ''
            if x[0] == 'pred':
                if isinstance(x[1], list):
                    x2 = list(x[1])
                if x[1] == 'T' or x2.__contains__('T'):
                    annotationDataRowBdeladdf.insert(2, 'Tolerant')

                else:
                    annotationDataRowBdeladdf.insert(2, 'Damaging')

                    self.damagingCount = self.damagingCount + 1
            if x[0] == 'rankscore':
                annotationDataRowBdeladdf.insert(0, x[1])
            if x[0] == 'score':
                annotationDataRowBdeladdf.insert(1, x[1])
            if len(annotationDataRowBdeladdf) == 3:
                annotationDataRowBdeladdf.insert(0, id)
                annotationDataRowBdeladdf.insert(1, hgvs37Id)
                finalAnnotationDataBdeladdf.append(annotationDataRowBdeladdf.copy())
                annotationDataRowBdeladdf.clear()
    def annotateByBdelnoaf(self, x, annotationDataRowBdelnoaf, finalAnnotationDataBdelnoaf, id,hgvs37Id):
        #if x[0] == 'no_af':
            x2 = ''
            if x[0] == 'pred':
                if isinstance(x[1], list):
                    x2 = list(x[1])
                if x[1] == 'T' or x2.__contains__('T'):
                    annotationDataRowBdelnoaf.insert(2, 'Tolerant')

                else:
                    annotationDataRowBdelnoaf.insert(2, 'Damaging')

                    self.damagingCount = self.damagingCount + 1
            if x[0] == 'rankscore':
                annotationDataRowBdelnoaf.insert(0, x[1])
            if x[0] == 'score':
                annotationDataRowBdelnoaf.insert(1, x[1])
            if len(annotationDataRowBdelnoaf) == 3:
                annotationDataRowBdelnoaf.insert(0, id)
                annotationDataRowBdelnoaf.insert(1, hgvs37Id)
                finalAnnotationDataBdelnoaf.append(annotationDataRowBdelnoaf.copy())
                annotationDataRowBdelnoaf.clear()
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
        if len(annotationDataRowVest4) == 3:
            annotationDataRowVest4.insert(0, id)
            annotationDataRowVest4.insert(1, hgvs37Id)
            finalAnnotationDataVest4.append(annotationDataRowVest4.copy())
            annotationDataRowVest4.clear()
    def annotateByDann(self, x, annotationDataRowDann, finalAnnotationDataDann, id,hgvs37Id):
        if x[0] == 'rankscore':
            annotationDataRowDann.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowDann.insert(1, x[1])
            if isinstance(x[1], list):
                tempX = x[1][0]
            else:
                tempX = x[1]
            if tempX<0.5:
                annotationDataRowDann.insert(2, 'Tolerant')

            else:
                annotationDataRowDann.insert(2, 'Damaging')

                self.damagingCount = self.damagingCount + 1
        if len(annotationDataRowDann) == 3:
            annotationDataRowDann.insert(0, id)
            annotationDataRowDann.insert(1, hgvs37Id)
            finalAnnotationDataDann.append(annotationDataRowDann.copy())
            annotationDataRowDann.clear()

    def annotateByEigen(self, x, annotationDataRowEigen, finalAnnotationDataEigen, id,hgvs37Id):
        if x[0] == 'phred_coding':
            annotationDataRowEigen.insert(0, x[1])
            if isinstance(x[1], list):
                tempX = x[1][0]
            else:
                tempX = x[1]
            if tempX<0.5:
                annotationDataRowEigen.insert(3, 'Tolerant')

            else:
                annotationDataRowEigen.insert(3, 'Damaging')

                self.damagingCount = self.damagingCount + 1
        if x[0] == 'raw_coding':
            annotationDataRowEigen.insert(1, x[1])
        if x[0] == 'raw_coding_rankscore':
            annotationDataRowEigen.insert(2, x[1])
        if len(annotationDataRowEigen) == 4:
            annotationDataRowEigen.insert(0, id)
            annotationDataRowEigen.insert(1, hgvs37Id)
            finalAnnotationDataEigen.append(annotationDataRowEigen.copy())
            annotationDataRowEigen.clear()

    def annotateByEigenpc(self, x, annotationDataRowEigenpc, finalAnnotationDataEigenpc, id,hgvs37Id):
        if x[0] == 'phred_coding':
            annotationDataRowEigenpc.insert(0, x[1])
            if isinstance(x[1], list):
                tempX = x[1][0]
            else:
                tempX = x[1]
            if tempX<0.5:
                annotationDataRowEigenpc.insert(3, 'Tolerant')

            else:
                annotationDataRowEigenpc.insert(3, 'Damaging')

                self.damagingCount = self.damagingCount + 1
        if x[0] == 'raw_coding':
            annotationDataRowEigenpc.insert(1, x[1])
        if x[0] == 'raw_coding_rankscore':
            annotationDataRowEigenpc.insert(2, x[1])
        if len(annotationDataRowEigenpc) == 4:
            annotationDataRowEigenpc.insert(0, id)
            annotationDataRowEigenpc.insert(1, hgvs37Id)
            finalAnnotationDataEigenpc.append(annotationDataRowEigenpc.copy())
            annotationDataRowEigenpc.clear()

    def annotateByDoegen2(self, x, annotationDataRowDoegen2, finalAnnotationDataDoegen2, id,hgvs37Id):
        x2 = ''
        if x[0] == 'pred':
            if isinstance(x[1], list):
                x2 = list(x[1])
            if x[1] == 'T' or x2.__contains__('T'):
                annotationDataRowDoegen2.insert(2, 'Tolerant')

            else:
                annotationDataRowDoegen2.insert(2, 'Damaging')

                self.damagingCount = self.damagingCount + 1
        if x[0] == 'rankscore':
            annotationDataRowDoegen2.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowDoegen2.insert(1, x[1])
        if len(annotationDataRowDoegen2) == 3:
            annotationDataRowDoegen2.insert(0, id)
            annotationDataRowDoegen2.insert(1, hgvs37Id)
            finalAnnotationDataDoegen2.append(annotationDataRowDoegen2.copy())
            annotationDataRowDoegen2.clear()

    def annotateByGenocanyon(self, x, annotationDataRowGenocanyon, finalAnnotationDataGenocanyon, id,hgvs37Id):
        if x[0] == 'rankscore':
            annotationDataRowGenocanyon.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowGenocanyon.insert(1, x[1])
            if isinstance(x[1], list):
                tempX = x[1][0]
            else:
                tempX = x[1]
            if tempX    <0.5:
                annotationDataRowGenocanyon.insert(2, 'Tolerant')

            else:
                annotationDataRowGenocanyon.insert(2, 'Damaging')

                self.damagingCount = self.damagingCount + 1
        if len(annotationDataRowGenocanyon) == 3:
            annotationDataRowGenocanyon.insert(0, id)
            annotationDataRowGenocanyon.insert(1, hgvs37Id)
            finalAnnotationDataGenocanyon.append(annotationDataRowGenocanyon.copy())
            annotationDataRowGenocanyon.clear()
    def getHGVS37Ids(self, varids):
        dbsnpvardata = DBSnpVarientDataRetrieval()
        grch37vardatainstance = Grch37VariantData()
        dbsnpvardata_str=''
        for varid in varids:
           print("varid", varid)
           try:
            dbsnpvardata_str = dbsnpvardata.getvariantdata(varid)
           except:
               print("NCBI server error: Data not found for", varid, ' rs id')
           grch37vardatainstance.parsevardatabystring(dbsnpvardata_str, varid)
        hgvs37Ids=grch37vardatainstance.hgvs37Ids
        return hgvs37Ids