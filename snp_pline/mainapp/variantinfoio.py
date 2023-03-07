from Bio import Entrez
import myvariant
from .dictionaryIO import DictionaryIO
import uuid
from .csvwriting import CSVFileWriting
import requests, sys, json
from .dbsnpvarientdataretrieval import DBSnpVarientDataRetrieval

from biothings_client import get_client
from .grch37variantdata import Grch37VariantData
class VarintInfoIO:
    dataNotFound=[]
    toolAnnonotationNotFound={}
    varDataNotFound=[]
    isPathogenicBySelectedTools=0
    pathogenicVariants=[]
    chekboxValues=[]
    isSift = False; isSift4G=False; isProvean=False; isMVP=False; isPolyphen2HVAR=False; isPolyphen2HDIV=False; isPrimateAI=False; isRevel=False;
    isMPC=False; isMutPred=False; isMTaster=False; isMAssessor=False; isMRNN=False; isMSVM=False; isMLR=False; isMCAP=False; isLS2=False; isFATHMM=False;
    isFXF=False; isFMKL=False; isBDelAddAF=False; isBDelNoAF=False; isVest4=False; isDANN=False; isEigen=False; isEigenPC=False; isDeogen2=False;
    isGenoCanyon=False; pathogenicVarFile='';variantNotFoundFile='';
    siftFile = '';  sift4gFile = ''; proveanFile = ''; mvpFile = '';polyphen2hvarFile = '';polyphen2hdivFile=''; primateaiFile = ''; revelFile = ''; mpcFile = '';
    mutpredFile = ''; mtasterFile = ''; massessorFile = ''; mrnnFile = ''; msvmFile = ''; mlrFile = ''; mcapFile = ''; ls2File = ''; fathmmFile = '';
    fxfFile = ''; fmklFile = ''; bdeladdafFile = ''; bdelnoafFile = ''; vest4File = ''; dannFile = ''; eigenFile = ''; eigenpcFile = '';  doegen2File = ''; genocanyonFile = '';
    siftPVars=0;sift4gPVars=0;proveanPVars=0;mvpPVars=0;polyphen2hvarPVars=0;polyphen2hdivPVars=0;primateaiPVars=0; revelPVars=0;mpcPVars=0;
    mutpredPVars = 0;mtasterPVars=0;massessorPVars=0;mrnnPVars=0;msvmPVars=0;mlrPVars=0;mcapPVars=0;ls2PVars=0;fathmmPVars=0;
    fxfPVars = 0;fmklPVars=0;bdeladdafPVars=0;bdelnoafPVars=0;vest4PVars=0;dannPVars=0;eigenPVars=0;eigenpcPVars=0;doegen2PVars=0;genocanyonPVars=0;

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
        self.isSift = False; self.isSift4G = False; self.isProvean = False; self.isMVP = False; self.isPolyphen2HVAR = False; self.isPolyphen2HDIV = False;
        self.isPrimateAI = False; self.isRevel = False; self.isMPC = False; self.isMutPred = False; self.isMTaster = False; self.isMAssessor = False;
        self.isMRNN = False; self.isMSVM = False; self.isMLR = False; self.isMCAP = False; self.isLS2 = False; self.isFATHMM = False; self.isFXF = False;
        self.isFMKL = False; self.isBDelAddAF = False; self.isBDelNoAF = False; self.isVest4 = False; self.isDANN = False; self.isEigen = False; self.isEigenPC = False;
        self.isDeogen2 = False;  self.isGenoCanyon = False;
        self.chekboxValues = str(request.POST.getlist('chekboxValues', ''))

        self.siftFile = ''; self.sift4gFile = ''; self.proveanFile = ''; self.mvpFile = ''; self.polyphen2hvarFile = ''; self.polyphen2hdivFile='';self.primateaiFile = ''; self.revelFile = ''; self.mpcFile = '';
        self.mutpredFile = ''; self.mtasterFile = ''; self.massessorFile = ''; self.mrnnFile = ''; self.msvmFile = ''; self.mlrFile = ''; self.mcapFile = ''; self.ls2File = '';
        self.fathmmFile = ''; self.fxfFile = ''; self.fmklFile = ''; self.bdeladdafFile = ''; self.bdelnoafFile = ''; self.vest4File = ''; self.dannFile = ''; self.eigenFile = '';
        self.eigenpcFile = '';  self.doegen2File = '';  self.genocanyonFile = ''; self.pathogenicVarFile='';self.variantNotFoundFile='';
        self.siftPVars = 0; self.sift4gPVars = 0;  self.proveanPVars = 0; self.mvpPVars = 0; self.polyphen2hvarPVars = 0; self.polyphen2hdivPVars = 0; self.primateaiPVars = 0;  self.revelPVars = 0;
        self.mpcPVars = 0;  self.mutpredPVars = 0;  self.mtasterPVars = 0; self.massessorPVars = 0; self.mrnnPVars = 0; self.msvmPVars = 0; self.mlrPVars = 0; self.mcapPVars = 0; self.ls2PVars = 0;
        self.fathmmPVars = 0; self.fxfPVars = 0; self.fmklPVars = 0; self.bdeladdafPVars = 0; self.bdelnoafPVars = 0; self.vest4PVars = 0; self.dannPVars = 0; self.eigenPVars = 0; self.eigenpcPVars = 0;
        self.doegen2PVars = 0;self.genocanyonPVars=0;
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
        pathogenicVariants_local=[]
        for vars in hgvs37Ids: #varids: for (a, b, c) in itertools.zip_longest(num, color, value, fillvalue=-1):
            id = vars[0]
            hgvs37Id = vars[1] #'rs' + str(varid)
            self.isPathogenicBySelectedTools = 0
            varDataNotFound_local.clear()
            pathogenicVariants_local.clear()
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

                if self.chekboxValues.__contains__('sft') and  'sift' in dbnsfpVarDic:
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('sift'), 'sift' ):
                        #print(x)
                        #if len(x)>0: #if x[0] == 'sift':
                        self.annotateBySift(x, annotationDataRowSift, finalAnnotationDataSift, id,hgvs37Id)
                    if loopCount==0:
                        self.siftFile='media/' + str(uuid.uuid4()) +'_sift' + '.csv'

                else:
                    self.toolAnnonotationNotFound['sift' + str(loopCount)] = id

                if self.chekboxValues.__contains__('sift4g') and  'sift4g' in dbnsfpVarDic:
                   for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('sift4g'), 'sift4g'):
                        #print(x)
                        #if len(x)>0:#if x[0] == 'sift4g':
                        self.annotateBySift4g(x,annotationDataRowSift4g,finalAnnotationDataSift4g,id,hgvs37Id)
                   if loopCount==0:
                        self.sift4gFile='media/' + str(uuid.uuid4()) +'_sift4g' + '.csv'

                else:
                    self.toolAnnonotationNotFound['sift4g' + str(loopCount)] = id

                if self.chekboxValues.__contains__('provean') and  'provean' in dbnsfpVarDic:
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('provean'), 'provean'):
                        #print(x)
                        #if len(x)>0:#if x[0] == 'provean':
                            self.annotateByProvean(x,annotationDataRowProvean,finalAnnotationDataProvean,id,hgvs37Id)
                    if loopCount==0:
                        self.proveanFile='media/' + str(uuid.uuid4()) +'_provean' + '.csv'

                else:
                    self.toolAnnonotationNotFound['provean' + str(loopCount)] = id

                if self.chekboxValues.__contains__('mvp') and  'mvp' in dbnsfpVarDic:
                   for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('mvp'), 'mvp'):
                        #print(x)
                        #if len(x)>0:#if x[0] == 'mvp':
                        self.annotateByMVP(x,annotationDataRowMVP,finalAnnotationDataMVP,id,hgvs37Id)
                   if loopCount==0:
                        self.mvpFile='media/' + str(uuid.uuid4()) +'_mvp' + '.csv'

                else:
                    self.toolAnnonotationNotFound['mvp' + str(loopCount)] = id

                if self.chekboxValues.__contains__('polyphen2hvar') and dbnsfpVarDic.__contains__('polyphen2'):
                    tempDic=dbnsfpVarDic.get('polyphen2')
                    #print('dbnsfpVarDic.get(hvar)==', dbnsfpVarDic.get('hvar'))
                    for x in dictIO.nested_dict_pairs_iterator(tempDic.get('hvar'), 'polyphen2hvar'): #Polyphen2hdiv
                        #print(x)
                        #if len(x)>0:# 'hvar':
                        self.annotateByPolyphen2hvar(x,annotationDataRowPolyphen2hvar,finalAnnotationDataPolyphen2hvar,id,hgvs37Id)

                    if loopCount==0:
                        self.polyphen2hvarFile='media/' + str(uuid.uuid4()) +'_polyphen2hvar' + '.csv'

                else:
                    self.toolAnnonotationNotFound['polyphen2hvar' + str(loopCount)] = id

                if self.chekboxValues.__contains__('polyphen2hdiv') and dbnsfpVarDic.__contains__('polyphen2'):
                    tempDic = dbnsfpVarDic.get('polyphen2')
                    #print('dbnsfpVarDic.get(hvar)==', tempDic.get('hdiv'))
                    for x in dictIO.nested_dict_pairs_iterator(tempDic.get('hdiv'),'polyphen2hdiv'): #Polyphen2hdiv
                       # print(x)
                        #if len(x)>0:#'hdiv':
                        self.annotateByPolyphen2hdiv(x,annotationDataRowPolyphen2hdiv,finalAnnotationDataPolyphen2hdiv,id,hgvs37Id)
                    if loopCount==0:
                        self.polyphen2hdivFile='media/' + str(uuid.uuid4()) +'_polyphen2hdiv' + '.csv'

                else:
                     self.toolAnnonotationNotFound['polyphen2hdiv' + str(loopCount)] = id

                if self.chekboxValues.__contains__('primateai') and  'primateai' in dbnsfpVarDic:
                   for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('primateai'),'primateai'):
                        #print(x)
                        #if len(x)>0:#if x[0] == 'primateai':
                        self.annotateByPrimateai(x,annotationDataRowPrimateai,finalAnnotationDataPrimateai,id,hgvs37Id)
                   if loopCount==0:
                        self.primateaiFile='media/' + str(uuid.uuid4()) +'_primateai' + '.csv'

                else:
                    self.toolAnnonotationNotFound['primateai' + str(loopCount)] = id

                if self.chekboxValues.__contains__('revel') and  'revel' in dbnsfpVarDic:
                   for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('revel'),'revel'):
                        #print(x)
                        #if len(x)>0:#if x[0] == 'primateai':
                        self.annotateByRevel(x,annotationDataRowRevel,finalAnnotationDataRevel,id,hgvs37Id)
                   if loopCount==0:
                        self.revelFile='media/' + str(uuid.uuid4()) +'_revel' + '.csv'

                else:
                    self.toolAnnonotationNotFound['revel' + str(loopCount)] = id

                if self.chekboxValues.__contains__('mpc') and  'mpc' in dbnsfpVarDic:
                   for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('mpc'),'mpc'):
                        #print(x)
                        #if len(x)>0:#if x[0] == 'primateai':
                        self.annotateByMPC(x,annotationDataRowMPC,finalAnnotationDataMPC,id,hgvs37Id)
                   if loopCount==0:
                        self.mpcFile='media/' + str(uuid.uuid4()) +'_mpc' + '.csv'

                else:
                    self.toolAnnonotationNotFound['mpc' + str(loopCount)] = id

                if self.chekboxValues.__contains__('mutpred') and  'mutpred' in dbnsfpVarDic:
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('mutpred'),'mutpred'):
                        #print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByMutpred(x, annotationDataRowMutpred, finalAnnotationDataMutpred, id,hgvs37Id)
                    if loopCount==0:
                        self.mutpredFile='media/' + str(uuid.uuid4()) +'_mutpred' + '.csv'

                else:
                    self.toolAnnonotationNotFound['mutpred' + str(loopCount)] = id

                if self.chekboxValues.__contains__('mtaster') and  'mutationtaster' in dbnsfpVarDic:
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('mutationtaster'),'mtaster'):
                        #print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                            self.annotateByMtaster(x, annotationDataRowMtaster, finalAnnotationDataMtaster, id,hgvs37Id)
                    if loopCount==0:
                        self.mtasterFile='media/' + str(uuid.uuid4()) +'_mtaster' + '.csv'

                else:
                    self.toolAnnonotationNotFound['mtaster' + str(loopCount)] = id

                if self.chekboxValues.__contains__('massessor') and  'mutationassessor' in dbnsfpVarDic:
                   for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('mutationassessor'),'massessor'):
                        #print(x)
                        if len(x) > 0:  # if x[0] == 'primateai':
                            self.annotateByMassessor(x, annotationDataRowMassessor, finalAnnotationDataMassessor, id,hgvs37Id)
                   if loopCount==0:
                        self.massessorFile='media/' + str(uuid.uuid4()) +'_massessor' + '.csv'

                else:
                    self.toolAnnonotationNotFound['massessor' + str(loopCount)] = id

                if self.chekboxValues.__contains__('mrnn') and  'metarnn' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.metarnn.pred, dbnsfp.metarnn.rankscore, dbnsfp.metarnn.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('metarnn'),'mrnn'):
                        #print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                            self.annotateByMRNN(x, annotationDataRowMRNN, finalAnnotationDataMRNN, id,hgvs37Id)

                    if loopCount==0:
                        self.mrnnFile='media/' + str(uuid.uuid4()) +'_mrnn' + '.csv'

                else:
                    self.toolAnnonotationNotFound['mrnn' + str(loopCount)] = id

                if self.chekboxValues.__contains__('msvm') and  'metasvm' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.metasvm.pred, dbnsfp.metasvm.rankscore, dbnsfp.metasvm.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('metasvm'),'msvm'):
                        #print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByMRNN(x, annotationDataRowMSVM, finalAnnotationDataMSVM, id,hgvs37Id)
                    if loopCount==0:
                        self.msvmFile='media/' + str(uuid.uuid4()) +'_msvm' + '.csv'

                else:
                    self.toolAnnonotationNotFound['msvm' + str(loopCount)] = id

                if self.chekboxValues.__contains__('mlr') and  'metalr' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.metalr.pred, dbnsfp.metalr.rankscore, dbnsfp.metalr.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('metalr'),'mlr'):
                        #print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByMLR(x, annotationDataRowMLR, finalAnnotationDataMLR, id,hgvs37Id)
                    if loopCount==0:
                        self.mlrFile='media/' + str(uuid.uuid4()) +'_mlr' + '.csv'

                else:
                    self.toolAnnonotationNotFound['mlr' + str(loopCount)] = id

                if self.chekboxValues.__contains__('mcap') and  'm-cap' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.m-cap.pred, dbnsfp.m-cap.rankscore, dbnsfp.m-cap.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('m-cap'),'mcap'):
                        #print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByMCAP(x, annotationDataRowMCAP, finalAnnotationDataMCAP, id,hgvs37Id)
                    if loopCount==0:
                            self.mcapFile='media/' + str(uuid.uuid4()) +'_mcap' + '.csv'

                else:
                    self.toolAnnonotationNotFound['mcap' + str(loopCount)] = id

                if self.chekboxValues.__contains__('ls2') and  'list-s2' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.list-s2.pred, dbnsfp.list-s2.rankscore, dbnsfp.list-s2.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('list-s2'),'ls2'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByLS2(x, annotationDataRowLS2, finalAnnotationDataLS2, id,hgvs37Id)
                    if loopCount==0:
                        self.ls2File='media/' + str(uuid.uuid4()) +'_ls2' + '.csv'

                else:
                    self.toolAnnonotationNotFound['ls2' + str(loopCount)] = id

                if self.chekboxValues.__contains__('fathmm') and  'fathmm' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.fathmm.pred, dbnsfp.fathmm.rankscore, dbnsfp.fathmm.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('fathmm'),'fathmm'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByFathmm(x, annotationDataRowFathmm, finalAnnotationDataFathmm, id,hgvs37Id)
                    if loopCount==0:
                        self.fathmmFile='media/' + str(uuid.uuid4()) +'_fathmm' + '.csv'

                else:
                    self.toolAnnonotationNotFound['fathmm' + str(loopCount)] = id
                if self.chekboxValues.__contains__('fxf') and  'fathmm-xf' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.fathmm-xf.coding_pred,dbnsfp.fathmm-xf.coding_rankscore, dbnsfp.fathmm-xf.coding_score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('fathmm-xf'),'fathmm-xf'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByFXF(x, annotationDataRowFXF, finalAnnotationDataFXF, id,hgvs37Id)
                    if loopCount==0:
                        self.fxfFile='media/' + str(uuid.uuid4()) +'_fxf' + '.csv'

                else:
                    self.toolAnnonotationNotFound['fxf' + str(loopCount)] = id

                if self.chekboxValues.__contains__('fmkl') and  'fathmm-mkl' in dbnsfpVarDic:
                    #fabricatedFields.append('fathmm-mkl.coding_pred,dbnsfp.fathmm-mkl.coding_rankscore,dbnsfp.fathmm-mkl.coding_score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('fathmm-mkl'),'fathmm-mkl'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByFMKL(x, annotationDataRowFMKL, finalAnnotationDataFMKL, id,hgvs37Id)
                    if loopCount==0:
                        self.fmklFile='media/' + str(uuid.uuid4()) +'_fmkl' + '.csv'

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
                    if loopCount==0:
                        self.bdeladdafFile='media/' + str(uuid.uuid4()) +'_bdeladdaf' + '.csv'

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
                    if loopCount==0:
                        self.bdelnoafFile='media/' + str(uuid.uuid4()) +'_bdelnoaf' + '.csv'

                else:
                    self.toolAnnonotationNotFound['bdelnoaf' + str(loopCount)] = id

                if self.chekboxValues.__contains__('vest4') and  'vest4' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.vest4.rankscore,dbnsfp.vest4.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('vest4'),'vest4'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByVest4(x, annotationDataRowVest4, finalAnnotationDataVest4, id,hgvs37Id)
                    if loopCount==0:
                        self.vest4File='media/' + str(uuid.uuid4()) +'_vest4' + '.csv'

                else:
                    self.toolAnnonotationNotFound['vest4' + str(loopCount)] = id

                if self.chekboxValues.__contains__('dann') and  'dann' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.dann.rankscore, dbnsfp.dann.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('dann'),'dann'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByDann(x, annotationDataRowDann, finalAnnotationDataDann, id,hgvs37Id)
                    if loopCount==0:
                        self.dannFile='media/' + str(uuid.uuid4()) +'_dann' + '.csv'

                else:
                    self.toolAnnonotationNotFound['dann' + str(loopCount)] = id

                if self.chekboxValues.__contains__('eign') and  'eigen' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.eigen.raw_coding, dbnsfp.eigen.raw_coding_rankscore')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('eigen'),'eigen'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByEigen(x, annotationDataRowEigen, finalAnnotationDataEigen, id,hgvs37Id)
                    if loopCount==0:
                        self.eigenFile='media/' + str(uuid.uuid4()) +'_eigen' + '.csv'

                else:
                    self.toolAnnonotationNotFound['eigen' + str(loopCount)] = id

                if self.chekboxValues.__contains__('eigenpc') and  'eigen-pc' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.eigen-pc.raw_coding, dbnsfp.eigen-pc.raw_coding_rankscore')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('eigen-pc'),'eigenpc'):
                        # print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByEigenpc(x, annotationDataRowEigenpc, finalAnnotationDataEigenpc, id,hgvs37Id)
                    if loopCount==0:
                        self.eigenpcFile='media/' + str(uuid.uuid4()) +'_eigenpc' + '.csv'

                else:
                    self.toolAnnonotationNotFound['eigenpc' + str(loopCount)] = id

                if self.chekboxValues.__contains__('deogen2') and  'deogen2' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.deogen2.pred,dbnsfp.deogen2.rankscore, dbnsfp.deogen2.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('deogen2'),'doegen2'):
                        # print(x)
                        if len(x) > 0:  # if x[0] == 'primateai':
                            self.annotateByDoegen2(x, annotationDataRowDoegen2, finalAnnotationDataDoegen2, id,hgvs37Id)
                    if loopCount==0:
                        self.doegen2File='media/' + str(uuid.uuid4()) +'_deogen2' + '.csv'

                else:
                    self.toolAnnonotationNotFound['deogen2' + str(loopCount)] = id

                if self.chekboxValues.__contains__('genocanyon') and  'genocanyon' in dbnsfpVarDic:
                    #fabricatedFields.append('dbnsfp.genocanyon.rankscore,dbnsfp.genocanyon.score')
                    for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('genocanyon'),'genocanyon'):
                        #print(x)
                        #if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByGenocanyon(x, annotationDataRowGenocanyon, finalAnnotationDataGenocanyon, id,hgvs37Id)

                    if loopCount==0:
                        self.genocanyonFile='media/' + str(uuid.uuid4()) +'_genocanyon' + '.csv'


                else:
                    self.toolAnnonotationNotFound['genocanyon' + str(loopCount)] = id
                loopCount +=1
            else:
                varDataNotFound_local.append(id)
                varDataNotFound_local.append(hgvs37Id)
                self.varDataNotFound.append(varDataNotFound_local.copy())
                continue
            if self.isPathogenicBySelectedTools==2:
                pathogenicVariants_local.append(id)
                pathogenicVariants_local.append(hgvs37Id)
                self.pathogenicVariants.append(pathogenicVariants_local.copy())
        ###ends outmost for loop
        if len(finalAnnotationDataSift) > 0:
            fieldnames = ['rsid','HGVSId', 'Converted Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataSift, self.siftFile,fieldnames)
            self.isSift = True
        else:
            self.dataNotFound.append('sift')

        if len(finalAnnotationDataSift4g) > 0:
            fieldnames = ['rsid','HGVSId', 'Converted Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataSift4g, self.sift4gFile,fieldnames)
            self.isSift4G=True
        else:
            self.dataNotFound.append('sift4g')

        if len(finalAnnotationDataProvean) > 0:
            fieldnames = ['rsid','HGVSId', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataProvean, self.proveanFile, fieldnames)
            self.isProvean=True
        else:
            self.dataNotFound.append('provean')

        if len(finalAnnotationDataMVP) > 0:
            fieldnames = ['rsid','HGVSId', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMVP, self.mvpFile, fieldnames)
            self.isMVP=True
        else:
            self.dataNotFound.append('mvp')


        if len(finalAnnotationDataPolyphen2hvar) > 0:
            fieldnames = ['rsid','HGVSId', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataPolyphen2hvar, self.polyphen2hvarFile,fieldnames)
            self.isPolyphen2HVAR=True
        else:
            self.dataNotFound.append('polyphen2hvar')

        if len(finalAnnotationDataPolyphen2hdiv) > 0:
            fieldnames = ['rsid','HGVSId', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataPolyphen2hdiv, self.polyphen2hdivFile,fieldnames)
            self.isPolyphen2HDIV=True
        else:
            self.dataNotFound.append('polyphen2hdiv')

        if len(finalAnnotationDataPrimateai) > 0:
            fieldnames = ['rsid','HGVSId', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataPrimateai, self.primateaiFile,fieldnames)
            self.isPrimateAI=True
        else:
            self.dataNotFound.append('primateai') #

        if len(finalAnnotationDataRevel) > 0:
            fieldnames = ['rsid','HGVSId','Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataRevel, self.revelFile,fieldnames)
            self.isRevel=True
        else:
            self.dataNotFound.append('revel')

        if len(finalAnnotationDataMPC) > 0:
            fieldnames = ['rsid','HGVSId', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMPC, self.mpcFile, fieldnames)
            self.isMPC=True
        else:
            self.dataNotFound.append('mpc')  #

        if len(finalAnnotationDataMutpred) > 0:
            fieldnames = ['rsid','HGVSId', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMutpred, self.mutpredFile, fieldnames)
            self.isMutPred=True
        else:
            self.dataNotFound.append('mutpred')

        if len(finalAnnotationDataMtaster) > 0:
            fieldnames = ['rsid','HGVSId', 'Converted Rankscore', 'Score', 'Prediction']
            #print('finalAnnotationDataMtaster== ',finalAnnotationDataMtaster)
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMtaster, self.mtasterFile, fieldnames)
            self.isMTaster=True
        else:
            self.dataNotFound.append('mtaster')  # 'finalAnnotationDataMtaster'

        if len(finalAnnotationDataMassessor) > 0:
            fieldnames = ['rsid','HGVSId', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMassessor, self.massessorFile, fieldnames)
            self.isMAssessor=True
        else:
            self.dataNotFound.append('massessor')

        if len(finalAnnotationDataMRNN) > 0:
            fieldnames = ['rsid','HGVSId', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMRNN, self.mrnnFile, fieldnames)
            self.isMRNN=True
        else:
            self.dataNotFound.append('mrnn')

        if len(finalAnnotationDataMSVM) > 0:
            fieldnames = ['rsid','HGVSId', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMSVM, self.msvmFile, fieldnames)
            self.isMSVM=True
        else:
            self.dataNotFound.append('msvm')

        if len(finalAnnotationDataMLR) > 0:
            fieldnames = ['rsid','HGVSId', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMLR, self.mlrFile, fieldnames)
            self.isMLR=True
        else:
            self.dataNotFound.append('mlr')

        if len(finalAnnotationDataMCAP) > 0:
            fieldnames = ['rsid','HGVSId', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMCAP, self.mcapFile, fieldnames)
            self.isMCAP=True
        else:
            self.dataNotFound.append('mcap')

        if len(finalAnnotationDataLS2) > 0:
            fieldnames = ['rsid','HGVSId', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataLS2, self.ls2File, fieldnames)
            self.isLS2-=True
        else:
            self.dataNotFound.append('ls2')

        if len(finalAnnotationDataFathmm) > 0:
            fieldnames = ['rsid','HGVSId', 'Converted Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataFathmm, self.fathmmFile, fieldnames)
            self.isFATHMM=True
        else:
            self.dataNotFound.append('fathmm')

        if len(finalAnnotationDataFXF) > 0:
            fieldnames = ['rsid','HGVSId', 'Coding Rankscore', 'CodingScore', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataFXF, self.fxfFile, fieldnames)
            self.isFXF=True
        else:
            self.dataNotFound.append('fxf')
        if len(finalAnnotationDataFMKL) > 0:
            fieldnames = ['rsid','HGVSId', 'Coding Rankscore', 'Coding Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataFMKL, self.fmklFile, fieldnames)
            self.isFMKL=True
        else:
            self.dataNotFound.append('fmkl') #finalAnnotationDataBdeladdaf

        if len(finalAnnotationDataBdeladdaf) > 0:
            fieldnames = ['rsid', 'HGVSId','Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataBdeladdaf, self.bdeladdafFile, fieldnames)
            self.isBDelAddAF=True
        else:
            self.dataNotFound.append('bdeladdaf')

        if len(finalAnnotationDataBdelnoaf) > 0:
            fieldnames = ['rsid','HGVSId','Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataBdelnoaf, self.bdelnoafFile, fieldnames)
            self.isBDelNoAF=True
        else:
            self.dataNotFound.append('bdelnoaf')

        if len(finalAnnotationDataVest4) > 0:
            fieldnames = ['rsid','HGVSId', 'Rankscore', 'Score']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataVest4, self.vest4File, fieldnames)
            self.isVest4=True
        else:
            self.dataNotFound.append('vest4')

        if len(finalAnnotationDataDann) > 0:
            fieldnames = ['rsid','HGVSId', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataDann, self.dannFile, fieldnames)
            self.isDANN=True
        else:
            self.dataNotFound.append('dann')

        if len(finalAnnotationDataDoegen2) > 0:
            fieldnames = ['rsid','HGVSId', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataDoegen2, self.doegen2File, fieldnames)
            self.isDeogen2=True
        else:
            self.dataNotFound.append('doegen2')

        if len(finalAnnotationDataGenocanyon) > 0:
            fieldnames = ['rsid','HGVSId', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataGenocanyon, self.genocanyonFile, fieldnames)
            self.isGenoCanyon=True
        else:
            self.dataNotFound.append('genocanyon')

        if len(finalAnnotationDataEigen) > 0:
            fieldnames = ['rsid','HGVSId', 'Phred Coding', 'Raw Coding','Raw Coding Rankscore' , 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataEigen, self.eigenFile, fieldnames)
            self.isEigen=True
        else:
            self.dataNotFound.append('eigen')

        if len(finalAnnotationDataEigenpc) > 0:
            fieldnames = ['rsid','HGVSId', 'Phred Coding', 'Raw Coding','Raw Coding Rankscore', 'Prediction' ]
            csvfw.writeAnnotationDateCSV(finalAnnotationDataEigenpc, self.eigenpcFile, fieldnames)
            self.isEigenPC=True
        else:
            self.dataNotFound.append('eigenpc')
        if len(self.pathogenicVariants)>0:
            fieldnames = ['rsid', 'HGVSId']
            self.pathogenicVarFile = 'media/' + str(uuid.uuid4()) + '_pathogenic_variants' + '.csv'
            csvfw.writeAnnotationDateCSV(self.pathogenicVariants, self.pathogenicVarFile, fieldnames)

        if len(self.varDataNotFound)>0:
            fieldnames = ['rsid', 'HGVSId']
            self.variantNotFoundFile = 'media/' + str(uuid.uuid4()) + '_variantnotfound' + '.csv'
            csvfw.writeAnnotationDateCSV(self.varDataNotFound, self.variantNotFoundFile, fieldnames)
    ###ends startParsing method

    def annotateBySift(self, x, annotationDataRowSift,finalAnnotationDataSift,id,hgvs37Id):
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDataRowSift.insert(2, 'Tolerant')
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowSift.insert(2, 'Damaging')
                self.siftPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools=  2
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
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDataRowSift4g.insert(2, 'Tolerant')
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowSift4g.insert(2, 'Damaging')
                self.sift4gPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
        if x[0] == 'pred':
            if x[1] == 'N':
                annotationDataRowProvean.insert(2, 'Tolerant')
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowProvean.insert(2, 'Damaging')
                self.proveanPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowMVP.insert(2, 'Damaging')
                self.mvpPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
        if len(annotationDataRowMVP) == 3:
            annotationDataRowMVP.insert(0, id)
            annotationDataRowMVP.insert(1, hgvs37Id)
            finalAnnotationDataMVP.append(annotationDataRowMVP.copy())
            annotationDataRowMVP.clear()
    #Polyphen2hdiv
    def annotateByPolyphen2hvar(self, x, annotationDataRowPolyphen2hvar, finalAnnotationDataPolyphen2hvar, id,hgvs37Id):
        if x[0] == 'pred':
            if x[1] == 'B':
                annotationDataRowPolyphen2hvar.insert(2, 'Tolerant')
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowPolyphen2hvar.insert(2, 'Damaging')
                self.polyphen2hvarPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowPolyphen2hdiv.insert(2, 'Damaging')
                self.polyphen2hdivPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowPrimateai.insert(2, 'Damaging')
                self.primateaiPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowRevel.insert(2, 'Damaging')
                self.revelPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowMPC.insert(2, 'Damaging')
                self.mpcPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowMutpred.insert(2, 'Damaging')
                self.mutpredPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
        if len(annotationDataRowMutpred) == 3:
            annotationDataRowMutpred.insert(0, id)
            annotationDataRowMutpred.insert(1, hgvs37Id)
            finalAnnotationDataMutpred.append(annotationDataRowMutpred.copy())
            annotationDataRowMutpred.clear()

    def annotateByMtaster(self, x, annotationDataRowMtaster, finalAnnotationDataMtaster, id,hgvs37Id):
        if x[0] == 'pred':
            if x[1] == 'N':
                annotationDataRowMtaster.insert(2, 'Tolerant')
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowMtaster.insert(2, 'Damaging')
                self.mtasterPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
        if x[0] == 'pred':
            if x[1] == 'N':
                annotationDataRowMassessor.insert(2, 'Tolerant')
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowMassessor.insert(2, 'Damaging')
                self.massessorPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDatRowMRNN.insert(2, 'Tolerant')
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDatRowMRNN.insert(2, 'Damaging')
                self.mrnnPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDatRowMSVM.insert(2, 'Tolerant')
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDatRowMSVM.insert(2, 'Damaging')
                self.msvmPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDatRowMLR.insert(2, 'Tolerant')
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDatRowMLR.insert(2, 'Damaging')
                self.mlrPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDatRowMCAP.insert(2, 'Tolerant')
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDatRowMCAP.insert(2, 'Damaging')
                self.mcapPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDatRowLS2.insert(2, 'Tolerant')
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDatRowLS2.insert(2, 'Damaging')
                self.ls2PVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDatRowFathmm.insert(2, 'Tolerant')
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDatRowFathmm.insert(2, 'Damaging')
                self.fathmmPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
        if x[0] == 'coding_pred':
            if x[1] == 'N':
                annotationDatRowFXF.insert(2, 'Tolerant')
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDatRowFXF.insert(2, 'Damaging')
                self.fxfPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
        if x[0] == 'coding_pred':
            if x[1] == 'N':
                annotationDatRowFMKL.insert(2, 'Tolerant')
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDatRowFMKL.insert(2, 'Damaging')
                self.fmklPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
            if x[0] == 'pred':
                if x[1] == 'T':
                    annotationDataRowBdeladdf.insert(2, 'Tolerant')
                    self.isPathogenicBySelectedTools = 1
                else:
                    annotationDataRowBdeladdf.insert(2, 'Damaging')
                    self.bdeladdafPVars+=1
                    if self.isPathogenicBySelectedTools != 1:
                        self.isPathogenicBySelectedTools = 2
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
            if x[0] == 'pred':
                if x[1] == 'T':
                    annotationDataRowBdelnoaf.insert(2, 'Tolerant')
                    self.isPathogenicBySelectedTools = 1
                else:
                    annotationDataRowBdelnoaf.insert(2, 'Damaging')
                    self.bdelnoafPVars+=1
                    if self.isPathogenicBySelectedTools != 1:
                        self.isPathogenicBySelectedTools = 2
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
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowVest4.insert(2, 'Damaging')
                self.vest4PVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowDann.insert(2, 'Damaging')
                self.dannPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowEigen.insert(3, 'Damaging')
                self.eigenPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowEigenpc.insert(3, 'Damaging')
                self.eigenpcPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDataRowDoegen2.insert(2, 'Tolerant')
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowDoegen2.insert(2, 'Damaging')
                self.doegen2PVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
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
                self.isPathogenicBySelectedTools = 1
            else:
                annotationDataRowGenocanyon.insert(2, 'Damaging')
                self.genocanyonPVars+=1
                if self.isPathogenicBySelectedTools != 1:
                    self.isPathogenicBySelectedTools = 2
        if len(annotationDataRowGenocanyon) == 3:
            annotationDataRowGenocanyon.insert(0, id)
            annotationDataRowGenocanyon.insert(1, hgvs37Id)
            finalAnnotationDataGenocanyon.append(annotationDataRowGenocanyon.copy())
            annotationDataRowGenocanyon.clear()
    def getHGVS37Ids(self, varids):
        dbsnpvardata = DBSnpVarientDataRetrieval()
        grch37vardatainstance = Grch37VariantData()
        for varid in varids:
           dbsnpvardata_str = dbsnpvardata.getvariantdata(varid)
           grch37vardatainstance.parsevardatabystring(dbsnpvardata_str, varid)
        hgvs37Ids=grch37vardatainstance.hgvs37Ids
        return hgvs37Ids