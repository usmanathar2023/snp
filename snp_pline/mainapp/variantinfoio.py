from Bio import Entrez
import myvariant
from .dictionaryIO import DictionaryIO
import uuid
from .csvwriting import CSVFileWriting
import requests, sys, json
from biothings_client import get_client
class VarintInfoIO:
    dataNotFound=[]
    variantNotFound={}
    def parseMyVarinats(self, request,fabricatedTerm):
        dictIO = DictionaryIO()
        csvfw = CSVFileWriting()
        mv = myvariant.MyVariantInfo()
        mvVar = ''
        Entrez.email = "usman.athar@gmail.com"
        handle = Entrez.esearch(db="snp", term=fabricatedTerm, retmax=2)
        variantData = Entrez.read(handle)
        totalSNPs = variantData['Count']
        varids = variantData["IdList"]
        chekboxValues = str(request.POST.getlist('chekboxValues', ''))
        siftFile='';sift4gFile=''; proveanFile='';  caddFile='';mvpFile='';polyphen2hvarFile=''; polyphen2hdiv='';
        primateaiFile='';revelFile=''; mpcFile=''; mutpredFile='';mtasterFile=''; massessorFile=''; mrnnFile=''; msvmFile='';
        mlrFile=''; mcapFile=''; lrtFile=''; ls2File=''; fathmmFile=''; fxfFile=''; fmklFile=''; bdeladdafFile=''; bdelnoafFile='';
        aloftFile=''; vest4File=''; dannFile=''; eigenFile=''; eigenpcFile=''; doegen2File=''; genocanyonFile='';
        rsIds = []
        var = ''
        loopCount=0;
        print('chekboxValues=== ', chekboxValues)
        fabricatedFields=[]
        listOfTuplesSift = []
        annotationDataRowSift = []
        finalAnnotationDataSift = []
        annotationDataRowSift4g = []
        finalAnnotationDataSift4g = []
        annotationDataRowProvean = []
        finalAnnotationDataProvean = []
        annotationDataRowMVP = []
        finalAnnotationDataMVP = []
        annotationDataRowPolyphen2hvar = []
        finalAnnotationDataPolyphen2hvar = []
        annotationDataRowPolyphen2hdiv = []
        finalAnnotationDataPolyphen2hdiv = []
        annotationDataRowPrimateai = []
        finalAnnotationDataPrimateai = []
        annotationDataRowRevel = []
        finalAnnotationDataRevel = []
        annotationDataRowMPC = []
        finalAnnotationDataMPC = []
        annotationDataRowMutpred = []
        finalAnnotationDataMutpred = []
        annotationDataRowMtaster = []
        finalAnnotationDataMtaster = []
        annotationDataRowMassessor = []
        finalAnnotationDataMassessor = []
        annotationDataRowMRNN = []
        finalAnnotationDataMRNN = []
        annotationDataRowMSVM = []
        finalAnnotationDataMSVM = []
        annotationDataRowMLR = []
        finalAnnotationDataMLR = []
        annotationDataRowMCAP = []
        finalAnnotationDataMCAP = []
        annotationDataRowLS2 = []
        finalAnnotationDataLS2 = []
        annotationDataRowFathmm = []
        finalAnnotationDataFathmm = []
        annotationDataRowFXF = []
        finalAnnotationDataFXF = []
        annotationDataRowFMKL = []
        finalAnnotationDataFMKL = []
        annotationDataRowBdeladdaf = []
        finalAnnotationDataBdeladdaf = []
        annotationDataRowBdelnoaf = []
        finalAnnotationDataBdelnoaf = []
        annotationDataRowAloft = []
        finalAnnotationDataAloft = []
        annotationDataRowVest4 = []
        finalAnnotationDataVest4 = []
        annotationDataRowDann = []
        finalAnnotationDataDann = []
        annotationDataRowEigen = []
        finalAnnotationDataEigen = []
        annotationDataRowEigenpc = []
        finalAnnotationDataEigenpc = []
        annotationDataRowDoegen2 = []
        finalAnnotationDataDoegen2 = []
        annotationDataRowGenocanyon = []
        finalAnnotationDataGenocanyon = []

        for varid in varids:
            id = 'rs' + str(varid)
            r = requests.get('http://myvariant.info/v1/variant/'+id)  # , headers={"Content-Type": "application/json"})
            varDic = r.json()
            dbnsfpVarDic=varDic.get('dbnsfp')
            print('dbnsfpVarDic= ',dbnsfpVarDic)


            fabricatedFields.clear()
            '''
            if chekboxValues.__contains__('All'):
                print('fld=== ', mv.getvariant(id, fields=['dbnsfp.cadd','dbnsfp.cadd.pred', 'dbnsfp.cadd.raw_rankscore',
                                                           'dbnsfp.cadd.raw_score',
                                                           'dbnsfp.sift','dbnsfp.sift.pred',
                                                           'dbnsfp.sift.converted_rankscore', 'dbnsfp.sift.score',
                                                           'dbnsfp.sift4g,dbnsfp.sift4g.pred',
                                                           'dbnsfp.sift4g.converted_rankscore', 'dbnsfp.sift4g.score',
                                                           'dbnsfp.provean,dbnsfp.provean.pred',
                                                           'dbnsfp.provean.rankscore', 'dbnsfp.provean.score',
                                                           'dbnsfp.primateai,dbnsfp.primateai.pred',
                                                           'dbnsfp.primateai.rankscore', 'dbnsfp.primateai.score',
                                                           'dbnsfp.polyphen2.hdiv,dbnsfp.polyphen2.hdiv.pred',
                                                           'dbnsfp.polyphen2.hdiv.rankscore',
                                                           'dbnsfp.polyphen2.hdiv.score',
                                                           'dbnsfp.polyphen2.hvar', 'dbnsfp.polyphen2.hvar.pred',
                                                           'dbnsfp.polyphen2.hvar.rankscore',
                                                           'dbnsfp.polyphen2.hvar.score',
                                                           'dbnsfp.revel', 'dbnsfp.revel.rankscore',
                                                           'dbnsfp.revel.score',
                                                           'dbnsfp.mvp', 'dbnsfp.mvp.rankscore', 'dbnsfp.mvp.score',
                                                           'dbnsfp.mpc', 'dbnsfp.mpc.rankscore', 'dbnsfp.mpc.score',
                                                           'dbnsfp.mutpred,dbnsfp.mutpred.pred',
                                                           'dbnsfp.mutpred.rankscore', 'dbnsfp.mutpred.score',
                                                           'dbnsfp.mutationtaster,dbnsfp.mutationtaster.pred',
                                                           'dbnsfp.mutationtaster.converted_rankscore',
                                                           'dbnsfp.mutationtaster.score',
                                                           'dbnsfp.mutationassessor,dbnsfp.mutationassessor.pred',
                                                           'dbnsfp.mutationassessor.rankscore',
                                                           'dbnsfp.mutationassessor.score',
                                                           'dbnsfp.metarnn,dbnsfp.metarnn.pred',
                                                           'dbnsfp.metarnn.rankscore', 'dbnsfp.metarnn.score',
                                                           'dbnsfp.metasvm,dbnsfp.metasvm.pred',
                                                           'dbnsfp.metasvm.rankscore', 'dbnsfp.metasvm.score',
                                                           'dbnsfp.metalr,dbnsfp.metalr.pred',
                                                           'dbnsfp.metalr.rankscore', 'dbnsfp.metalr.score',
                                                           'dbnsfp.m-cap,dbnsfp.m-cap.pred', 'dbnsfp.m-cap.rankscore',
                                                           'dbnsfp.m-cap.score',
                                                           'dbnsfp.lrt,dbnsfp.lrt.pred',
                                                           'dbnsfp.lrt.converted_rankscore', 'dbnsfp.lrt.score',
                                                           'dbnsfp.list-s2,dbnsfp.list-s2.pred',
                                                           'dbnsfp.list-s2.rankscore', 'dbnsfp.list-s2.score',

                                                           'dbnsfp.fathmm,dbnsfp.fathmm.pred',
                                                           'dbnsfp.fathmm.rankscore', 'dbnsfp.fathmm.score',
                                                           'dbnsfp.fathmm-mkl,dbnsfp.fathmm-mkl.coding_pred',
                                                           'dbnsfp.fathmm-mkl.coding_rankscore',
                                                           'dbnsfp.fathmm-mkl.coding_score',
                                                           'dbnsfp.fathmm-xf,dbnsfp.fathmm-xf.coding_pred',
                                                           'dbnsfp.fathmm-xf.coding_rankscore',
                                                           'dbnsfp.fathmm-xf.coding_score',
                                                           'dbnsfp.bayesdel.add_af,dbnsfp.bayesdel.add_af.pred',
                                                           'dbnsfp.bayesdel.add_af.rankscore',
                                                           'dbnsfp.bayesdel.add_af.score',
                                                           'dbnsfp.bayesdel.no_af,dbnsfp.bayesdel.no_af.pred',
                                                           'dbnsfp.bayesdel.no_af.rankscore',
                                                           'dbnsfp.bayesdel.no_af.score',
                                                           'dbnsfp.aloft,dbnsfp.aloft.pred', 'dbnsfp.aloft.confidence',
                                                           'dbnsfp.vest4', 'dbnsfp.vest4.rankscore',
                                                           'dbnsfp.vest4.score',
                                                           'dbnsfp.dann', 'dbnsfp.dann.rankscore', 'dbnsfp.dann.score',
                                                           'dbnsfp.eigen', 'dbnsfp.eigen.raw_coding',
                                                           'dbnsfp.eigen.raw_coding_rankscore',
                                                           'dbnsfp.eigen-pc', 'dbnsfp.eigen-pc.raw_coding',
                                                           'dbnsfp.eigen-pc.raw_coding_rankscore',
                                                           'dbnsfp.deogen2', 'dbnsfp.deogen2.pred',
                                                           'dbnsfp.deogen2.rankscore', 'dbnsfp.deogen2.score',
                                                           'dbnsfp.genocanyon', 'dbnsfp.genocanyon.rankscore',
                                                           'dbnsfp.genocanyon.score', ]))
             '''
            if chekboxValues.__contains__('sift'):
                #fabricatedFields.append('dbnsfp.sift.pred, dbnsfp.sift.converted_rankscore, dbnsfp.sift.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('sift')):
                    #print(x)
                    if len(x)>0: #if x[0] == 'sift':
                        self.annotateBySift(x, annotationDataRowSift, finalAnnotationDataSift, id)
                    else:
                        self.variantNotFound['sift' + str(loopCount)] = id
                if loopCount==0:
                    siftFile='media/' + str(uuid.uuid4()) +'_sift' + '.csv'

                #fw = FileWriting()
                #fw.writeProteinData(protData, prot_var_data.prot_var_data_dict, prot_var_data_filename)
            if chekboxValues.__contains__('sift4g'):
                #fabricatedFields.append('dbnsfp.sift4g.pred, dbnsfp.sift4g.converted_rankscore, dbnsfp.sift4g.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('sift4g')):
                    #print(x)
                    if len(x)>0:#if x[0] == 'sift4g':
                        self.annotateBySift4g(x,annotationDataRowSift4g,finalAnnotationDataSift4g,id)
                    else:
                        self.variantNotFound['sift4g'+str(loopCount)]=id
                if loopCount==0:
                    sift4gFile='media/' + str(uuid.uuid4()) +'_sift4g' + '.csv'
            if chekboxValues.__contains__('provean'):
                #fabricatedFields.append('dbnsfp.provean.pred,dbnsfp.provean.rankscore, dbnsfp.provean.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('provean')):
                    print(x)
                    if len(x)>0:#if x[0] == 'provean':
                        self.annotateByProvean(x,annotationDataRowProvean,finalAnnotationDataProvean,id)
                    else:
                        self.variantNotFound['provean'+str(loopCount)]=id
                if loopCount==0:
                    proveanFile='media/' + str(uuid.uuid4()) +'_provean' + '.csv'
            if chekboxValues.__contains__('mvp'):
                #fabricatedFields.append('dbnsfp.mvp.rankscore, dbnsfp.mvp.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('mvp')):
                    print(x)
                    if len(x)>0:#if x[0] == 'mvp':
                        self.annotateByMVP(x,annotationDataRowMVP,finalAnnotationDataMVP,id)
                    else:
                        self.variantNotFound['mvp'+str(loopCount)]=id
                if loopCount==0:
                    mvpFile='media/' + str(uuid.uuid4()) +'_mvp' + '.csv'
            if chekboxValues.__contains__('polyphen2hvar'):
                #fabricatedFields.append('dbnsfp.polyphen2.hvar.pred, dbnsfp.polyphen2.hvar.rankscore,dbnsfp.polyphen2.hvar.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('polyphen2')): #Polyphen2hdiv
                    print(x)
                    if x[0] == 'hvar':
                        self.annotateByPolyphen2hvar(x,annotationDataRowPolyphen2hvar,finalAnnotationDataPolyphen2hvar,id)
                    else:
                        self.variantNotFound['polyphen2hvar'+str(loopCount)]=id
                if loopCount==0:
                    polyphen2hvarFile='media/' + str(uuid.uuid4()) +'_polyphen2hvar' + '.csv'
            if chekboxValues.__contains__('polyphen2hdiv'):
                #fabricatedFields.append('dbnsfp.polyphen2.hdiv.pred,dbnsfp.polyphen2.hdiv.rankscore,dbnsfp.polyphen2.hdiv.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('polyphen2')): #Polyphen2hdiv
                    print(x)
                    if x[0] == 'hdiv':
                        self.annotateByPolyphen2hdiv(x,annotationDataRowPolyphen2hdiv,finalAnnotationDataPolyphen2hdiv,id)
                    else:
                        self.variantNotFound['polyphen2hdiv'+str(loopCount)]=id
                if loopCount==0:
                    polyphen2hdivFile='media/' + str(uuid.uuid4()) +'_polyphen2hdiv' + '.csv'
            if chekboxValues.__contains__('primateai'):
                #fabricatedFields.append('dbnsfp.primateai.pred, dbnsfp.primateai.rankscore, dbnsfp.primateai.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('primateai')):
                    print(x)
                    if len(x)>0:#if x[0] == 'primateai':
                        self.annotateByPrimateai(x,annotationDataRowPrimateai,finalAnnotationDataPrimateai,id)
                    else:
                        self.variantNotFound['primateai'+str(loopCount)]=id
                if loopCount==0:
                    primateaiFile='media/' + str(uuid.uuid4()) +'_primateai' + '.csv'
            if chekboxValues.__contains__('revel'):
                #fabricatedFields.append('dbnsfp.revel.rankscore, dbnsfp.revel.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('revel')):
                    #print(x)
                    if len(x)>0:#if x[0] == 'primateai':
                        self.annotateByRevel(x,annotationDataRowRevel,finalAnnotationDataRevel,id)
                    else:
                        self.variantNotFound['revel'+str(loopCount)]=id
                if loopCount==0:
                    revelFile='media/' + str(uuid.uuid4()) +'_revel' + '.csv'
            if chekboxValues.__contains__('mpc'):
                #fabricatedFields.append('dbnsfp.mpc.rankscore, dbnsfp.mpc.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('mpc')):
                    #print(x)
                    if len(x)>0:#if x[0] == 'primateai':
                        self.annotateByMPC(x,annotationDataRowMPC,finalAnnotationDataMPC,id)
                    else:
                        self.variantNotFound['mpc'+str(loopCount)]=id
                if loopCount==0:
                    mpcFile='media/' + str(uuid.uuid4()) +'_mpc' + '.csv'
            if chekboxValues.__contains__('mutpred'):
                #fabricatedFields.append('dbnsfp.mutpred.pred, dbnsfp.mutpred.rankscore, dbnsfp.mutpred.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('mutpred')):
                    #print(x)
                    if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByMutpred(x, annotationDataRowMutpred, finalAnnotationDataMutpred, id)
                    else:
                        self.variantNotFound['mutpred' + str(loopCount)] = id
                if loopCount==0:
                    mutpredFile='media/' + str(uuid.uuid4()) +'_mutpred' + '.csv'
            if chekboxValues.__contains__('mtaster'):
                #fabricatedFields.append('dbnsfp.mutationtaster.pred,dbnsfp.mutationtaster.converted_rankscore,dbnsfp.mutationtaster.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('mutationtaster')):
                    #print(x)
                    if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByMtaster(x, annotationDataRowMtaster, finalAnnotationDataMtaster, id)
                    else:
                        self.variantNotFound['mtaster' + str(loopCount)] = id
                if loopCount==0:
                    mtasterFile='media/' + str(uuid.uuid4()) +'_mtaster' + '.csv'
            if chekboxValues.__contains__('massessor'):
                #fabricatedFields.append('dbnsfp.mutationassessor.pred,dbnsfp.mutationassessor.rankscore,dbnsfp.mutationassessor.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('mutationassessor')):
                    #print(x)
                    if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByMassessor(x, annotationDataRowMassessor, finalAnnotationDataMassessor, id)
                    else:
                        self.variantNotFound['massessor' + str(loopCount)] = id
                if loopCount==0:
                    massessorFile='media/' + str(uuid.uuid4()) +'_massessor' + '.csv'
            if chekboxValues.__contains__('mrnn'):
                #fabricatedFields.append('dbnsfp.metarnn.pred, dbnsfp.metarnn.rankscore, dbnsfp.metarnn.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('metarnn')):
                    #print(x)
                    if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByMRNN(x, annotationDataRowMRNN, finalAnnotationDataMRNN, id)
                    else:
                        self.variantNotFound['mrnn' + str(loopCount)] = id
                if loopCount==0:
                    mrnnFile='media/' + str(uuid.uuid4()) +'_mrnn' + '.csv'
            if chekboxValues.__contains__('msvm'):
                #fabricatedFields.append('dbnsfp.metasvm.pred, dbnsfp.metasvm.rankscore, dbnsfp.metasvm.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('metasvm')):
                    #print(x)
                    if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByMRNN(x, annotationDataRowMSVM, finalAnnotationDataMSVM, id)
                    else:
                        self.variantNotFound['msvm' + str(loopCount)] = id
                if loopCount==0:
                    msvmFile='media/' + str(uuid.uuid4()) +'_msvm' + '.csv'
            if chekboxValues.__contains__('mlr'):
                #fabricatedFields.append('dbnsfp.metalr.pred, dbnsfp.metalr.rankscore, dbnsfp.metalr.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('metalr')):
                    #print(x)
                    if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByMLR(x, annotationDataRowMLR, finalAnnotationDataMLR, id)
                    else:
                        self.variantNotFound['mlr' + str(loopCount)] = id
                if loopCount==0:
                    mlrFile='media/' + str(uuid.uuid4()) +'_mlr' + '.csv'
            if chekboxValues.__contains__('mcap'):
                #fabricatedFields.append('dbnsfp.m-cap.pred, dbnsfp.m-cap.rankscore, dbnsfp.m-cap.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('m-cap')):
                    #print(x)
                    if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByMCAP(x, annotationDataRowMCAP, finalAnnotationDataMCAP, id)
                    else:
                        self.variantNotFound['mcap' + str(loopCount)] = id
                if loopCount==0:
                    mcapFile='media/' + str(uuid.uuid4()) +'_mcap' + '.csv'

            if chekboxValues.__contains__('ls2'):
                #fabricatedFields.append('dbnsfp.list-s2.pred, dbnsfp.list-s2.rankscore, dbnsfp.list-s2.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('list-s2')):
                    # print(x)
                    if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByLS2(x, annotationDataRowLS2, finalAnnotationDataLS2, id)
                    else:
                        self.variantNotFound['ls2' + str(loopCount)] = id
                if loopCount==0:
                    ls2File='media/' + str(uuid.uuid4()) +'_ls2' + '.csv'
            if chekboxValues.__contains__('fathmm'):
                #fabricatedFields.append('dbnsfp.fathmm.pred, dbnsfp.fathmm.rankscore, dbnsfp.fathmm.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('fathmm')):
                    # print(x)
                    if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByFathmm(x, annotationDataRowFathmm, finalAnnotationDataFathmm, id)
                    else:
                        self.variantNotFound['fathmm' + str(loopCount)] = id
                if loopCount==0:
                    fathmmFile='media/' + str(uuid.uuid4()) +'_fathmm' + '.csv'
            if chekboxValues.__contains__('fxf'):
                #fabricatedFields.append('dbnsfp.fathmm-xf.coding_pred,dbnsfp.fathmm-xf.coding_rankscore, dbnsfp.fathmm-xf.coding_score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('fathmm-xf')):
                    # print(x)
                    if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByFXF(x, annotationDataRowFXF, finalAnnotationDataFXF, id)
                    else:
                        self.variantNotFound['fxf' + str(loopCount)] = id
                if loopCount==0:
                    fxfFile='media/' + str(uuid.uuid4()) +'_fxf' + '.csv'
            if chekboxValues.__contains__('fmkl'):
                #fabricatedFields.append('fathmm-mkl.coding_pred,dbnsfp.fathmm-mkl.coding_rankscore,dbnsfp.fathmm-mkl.coding_score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('fathmm-mkl')):
                    # print(x)
                    if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByFMKL(x, annotationDataRowFMKL, finalAnnotationDataFMKL, id)
                    else:
                        self.variantNotFound['fmkl' + str(loopCount)] = id
                if loopCount==0:
                    fmklFile='media/' + str(uuid.uuid4()) +'_fmkl' + '.csv'
            if chekboxValues.__contains__('bdeladdaf'):
                #fabricatedFields.append('dbnsfp.bayesdel.add_af.pred, dbnsfp.bayesdel.add_af.rankscore,dbnsfp.bayesdel.add_af.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('bayesdel')):
                    print(x)
                    if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByBdeladdaf(x, annotationDataRowBdeladdaf, finalAnnotationDataBdeladdaf, id)
                    else:
                        self.variantNotFound['bdeladdaf' + str(loopCount)] = id
                if loopCount==0:
                    bdeladdafFile='media/' + str(uuid.uuid4()) +'_bdeladdaf' + '.csv'
            if chekboxValues.__contains__('bdelnoaf'):
                #fabricatedFields.append('dbnsfp.bayesdel.no_af.pred, dbnsfp.bayesdel.no_af.rankscore,dbnsfp.bayesdel.no_af.score')
                for x in dictIO.nested_dict_pairs_iterator(dbnsfpVarDic.get('bayesdel')):
                    # print(x)
                    if len(x) > 0:  # if x[0] == 'primateai':
                        self.annotateByBdelnoaf(x, annotationDataRowBdelnoaf, finalAnnotationDataBdelnoaf, id)
                    else:
                        self.variantNotFound['bdelnoaf' + str(loopCount)] = id
                if loopCount==0:
                    bdelnoafFile='media/' + str(uuid.uuid4()) +'_bdelnoaf' + '.csv'

            if chekboxValues.__contains__('vest4'):
                fabricatedFields.append('dbnsfp.vest4.rankscore,dbnsfp.vest4.score')
                if loopCount==0:
                    vest4File='media/' + str(uuid.uuid4()) +'_vest4' + '.csv'
            if chekboxValues.__contains__('dann'):
                fabricatedFields.append('dbnsfp.dann.rankscore, dbnsfp.dann.score')
                if loopCount==0:
                    dannFile='media/' + str(uuid.uuid4()) +'_dann' + '.csv'
            if chekboxValues.__contains__('eigen'):
                fabricatedFields.append('dbnsfp.eigen.raw_coding, dbnsfp.eigen.raw_coding_rankscore')
                if loopCount==0:
                    eigenFile='media/' + str(uuid.uuid4()) +'_eigen' + '.csv'
            if chekboxValues.__contains__('eigenpc'):
                fabricatedFields.append('dbnsfp.eigen-pc.raw_coding, dbnsfp.eigen-pc.raw_coding_rankscore')
                if loopCount==0:
                    eigenpcFile='media/' + str(uuid.uuid4()) +'_eigenpc' + '.csv'
            if chekboxValues.__contains__('doegen2'):
                fabricatedFields.append('dbnsfp.deogen2.pred,dbnsfp.deogen2.rankscore, dbnsfp.deogen2.score')
                if loopCount==0:
                    doegen2File='media/' + str(uuid.uuid4()) +'_doegen2' + '.csv'
            if chekboxValues.__contains__('genocanyon'):
                fabricatedFields.append('dbnsfp.genocanyon.rankscore,dbnsfp.genocanyon.score')
                if loopCount==0:
                    genocanyonFile='media/' + str(uuid.uuid4()) +'_genocanyon' + '.csv'
            annotationData=mv.getvariant(id, fields=fabricatedFields)
            #print('sift--- ', annotationData.get('dbnsfp'))

            loopCount +=1
        ###ends outmost for loop
        if len(finalAnnotationDataSift) > 0:
            fieldnames = ['rsid', 'Converted Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataSift, siftFile,fieldnames)
        else:
            self.dataNotFound.append('sift')
        if len(finalAnnotationDataSift4g) > 0:
            fieldnames = ['rsid', 'Converted Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataSift4g, sift4gFile,fieldnames)
        else:
            self.dataNotFound.append('sift4g')
        if len(finalAnnotationDataProvean) > 0:
            fieldnames = ['rsid', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataProvean, proveanFile, fieldnames)
        else:
            self.dataNotFound.append('provean')
        if len(finalAnnotationDataMVP) > 0:
            fieldnames = ['rsid', 'Rankscore', 'Score']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMVP, mvpFile, fieldnames)
        else:
            self.dataNotFound.append('mvp')

        if len(finalAnnotationDataPolyphen2hvar) > 0:
            fieldnames = ['rsid', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataPolyphen2hvar, polyphen2hvarFile,fieldnames)
        else:
            self.dataNotFound.append('polyphen2hvar')
        if len(finalAnnotationDataPolyphen2hdiv) > 0:
            fieldnames = ['rsid', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataPolyphen2hdiv, polyphen2hdivFile,fieldnames)
        else:
            self.dataNotFound.append('polyphen2hdiv')
        if len(finalAnnotationDataPrimateai) > 0:
            fieldnames = ['rsid', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataPrimateai, primateaiFile,fieldnames)
        else:
            self.dataNotFound.append('polyphen2hdiv') #
        if len(finalAnnotationDataRevel) > 0:
            fieldnames = ['rsid','Rankscore', 'Score']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataRevel, revelFile,fieldnames)
        else:
            self.dataNotFound.append('revel')
        if len(finalAnnotationDataMPC) > 0:
            fieldnames = ['rsid', 'Rankscore', 'Score']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMPC, mpcFile, fieldnames)
        else:
            self.dataNotFound.append('mpc')  #
        if len(finalAnnotationDataMutpred) > 0:
            fieldnames = ['rsid', 'Rankscore', 'Score']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMutpred, mutpredFile, fieldnames)
        else:
            self.dataNotFound.append('mutpred')
        if len(finalAnnotationDataMtaster) > 0:
            fieldnames = ['rsid', 'Converted Rankscore', 'Score', 'Prediction']
            print('finalAnnotationDataMtaster== ',finalAnnotationDataMtaster)
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMtaster, mtasterFile, fieldnames)
        else:
            self.dataNotFound.append('mtaster')  # 'finalAnnotationDataMtaster'
        if len(finalAnnotationDataMassessor) > 0:
            fieldnames = ['rsid', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMassessor, massessorFile, fieldnames)
        else:
            self.dataNotFound.append('massessor')
        if len(finalAnnotationDataMRNN) > 0:
            fieldnames = ['rsid', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMRNN, mrnnFile, fieldnames)
        else:
            self.dataNotFound.append('mrnn')
        if len(finalAnnotationDataMSVM) > 0:
            fieldnames = ['rsid', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMSVM, msvmFile, fieldnames)
        else:
            self.dataNotFound.append('msvm')
        if len(finalAnnotationDataMLR) > 0:
            fieldnames = ['rsid', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMLR, mlrFile, fieldnames)
        else:
            self.dataNotFound.append('mlr')
        if len(finalAnnotationDataMCAP) > 0:
            fieldnames = ['rsid', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataMCAP, mcapFile, fieldnames)
        else:
            self.dataNotFound.append('mcap')
        if len(finalAnnotationDataLS2) > 0:
            fieldnames = ['rsid', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataLS2, ls2File, fieldnames)
        else:
            self.dataNotFound.append('ls2')
        if len(finalAnnotationDataFathmm) > 0:
            fieldnames = ['rsid', 'Converted Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataFathmm, fathmmFile, fieldnames)
        else:
            self.dataNotFound.append('fathmm')
        if len(finalAnnotationDataFXF) > 0:
            fieldnames = ['rsid', 'Coding Rankscore', 'CodingScore', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataFXF, fxfFile, fieldnames)
        else:
            self.dataNotFound.append('fxf')
        if len(finalAnnotationDataFMKL) > 0:
            fieldnames = ['rsid', 'Coding Rankscore', 'Coding Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataFMKL, fmklFile, fieldnames)
        else:
            self.dataNotFound.append('fmkl') #finalAnnotationDataBdeladdaf
        if len(finalAnnotationDataBdeladdaf) > 0:
            fieldnames = ['rsid', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataBdeladdaf, bdeladdafFile, fieldnames)
        else:
            self.dataNotFound.append('bdeladdaf')
        if len(finalAnnotationDataBdelnoaf) > 0:
            fieldnames = ['rsid', 'Rankscore', 'Score', 'Prediction']
            csvfw.writeAnnotationDateCSV(finalAnnotationDataBdelnoaf, bdelnoafFile, fieldnames)
        else:
            self.dataNotFound.append('bdelnoaf')
        ###ends parseMyVarinats
    def annotateBySift(self, x, annotationDataRowSift,finalAnnotationDataSift,id):
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDataRowSift.insert(2, 'Tolerant')
            else:
                annotationDataRowSift.insert(2, 'Damaging')
        if x[0] == 'converted_rankscore':
            annotationDataRowSift.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowSift.insert(1, x[1])
        if len(annotationDataRowSift) == 3:
            annotationDataRowSift.insert(0, id)
            finalAnnotationDataSift.append(annotationDataRowSift.copy())
            annotationDataRowSift.clear()

    def annotateBySift4g(self, x,annotationDataRowSift4g,finalAnnotationDataSift4g, id):
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDataRowSift4g.insert(2, 'Tolerant')
            else:
                annotationDataRowSift4g.insert(2, 'Damaging')
        if x[0] == 'converted_rankscore':
            annotationDataRowSift4g.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowSift4g.insert(1, x[1])
        if len(annotationDataRowSift4g) == 3:
            annotationDataRowSift4g.insert(0, id)
            finalAnnotationDataSift4g.append(annotationDataRowSift4g.copy())
            annotationDataRowSift4g.clear()

    def annotateByProvean(self, x, annotationDataRowProvean, finalAnnotationDataProvean, id):
        if x[0] == 'pred':
            if x[1] == 'N':
                annotationDataRowProvean.insert(2, 'Tolerant')
            else:
                annotationDataRowProvean.insert(2, 'Damaging')
        if x[0] == 'converted_rankscore':
            annotationDataRowProvean.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowProvean.insert(1, x[1])
        if len(annotationDataRowProvean) == 3:
            annotationDataRowProvean.insert(0, id)
            finalAnnotationDataProvean.append(annotationDataRowProvean.copy())
            annotationDataRowProvean.clear()

    def annotateByMVP(self, x, annotationDataRowMVP, finalAnnotationDataMVP, id):
        if x[0] == 'rankscore':
            annotationDataRowMVP.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowMVP.insert(1, x[1])
        if len(annotationDataRowMVP) == 2:
            annotationDataRowMVP.insert(0, id)
            finalAnnotationDataMVP.append(annotationDataRowMVP.copy())
            annotationDataRowMVP.clear()
    #Polyphen2hdiv
    def annotateByPolyphen2hvar(self, x, annotationDataRowPolyphen2hvar, finalAnnotationDataPolyphen2hvar, id):
        if x[1] == 'pred':
            if x[2] == 'B':
                annotationDataRowPolyphen2hvar.insert(2, 'Tolerant')
            else:
                annotationDataRowPolyphen2hvar.insert(2, 'Damaging')
        if x[1] == 'rankscore':
            annotationDataRowPolyphen2hvar.insert(0, x[2])
        if x[1] == 'score':
            annotationDataRowPolyphen2hvar.insert(1, x[2])
        if len(annotationDataRowPolyphen2hvar) == 3:
            annotationDataRowPolyphen2hvar.insert(0, id)
            finalAnnotationDataPolyphen2hvar.append(annotationDataRowPolyphen2hvar.copy())
            annotationDataRowPolyphen2hvar.clear()

    def annotateByPolyphen2hdiv(self, x, annotationDataRowPolyphen2hdiv, finalAnnotationDataPolyphen2hdiv, id):
        if x[1] == 'pred':
            if x[2] == 'B':
                annotationDataRowPolyphen2hdiv.insert(2, 'Tolerant')
            else:
                annotationDataRowPolyphen2hdiv.insert(2, 'Damaging')
        if x[1] == 'rankscore':
            annotationDataRowPolyphen2hdiv.insert(0, x[2])
        if x[1] == 'score':
            annotationDataRowPolyphen2hdiv.insert(1, x[2])
        if len(annotationDataRowPolyphen2hdiv) == 3:
            annotationDataRowPolyphen2hdiv.insert(0, id)
            finalAnnotationDataPolyphen2hdiv.append(annotationDataRowPolyphen2hdiv.copy())
            annotationDataRowPolyphen2hdiv.clear()

    def annotateByPrimateai(self, x, annotationDataRowPrimateai, finalAnnotationDataPrimateai, id):
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDataRowPrimateai.insert(2, 'Tolerant')
            else:
                annotationDataRowPrimateai.insert(2, 'Damaging')
        if x[0] == 'rankscore':
            annotationDataRowPrimateai.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowPrimateai.insert(1, x[1])
        if len(annotationDataRowPrimateai) == 3:
            annotationDataRowPrimateai.insert(0, id)
            finalAnnotationDataPrimateai.append(annotationDataRowPrimateai.copy())
            annotationDataRowPrimateai.clear()
#annotationDataRowRevel
    def annotateByRevel(self, x, annotationDataRowRevel, finalAnnotationDataRevel, id):
        if x[0] == 'rankscore':
            annotationDataRowRevel.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowRevel.insert(1, x[1])
        if len(annotationDataRowRevel) == 2:
            annotationDataRowRevel.insert(0, id)
            finalAnnotationDataRevel.append(annotationDataRowRevel.copy())
            annotationDataRowRevel.clear()

    def annotateByMPC(self, x, annotationDataRowMPC, finalAnnotationDataMPC, id):
        if x[0] == 'rankscore':
            annotationDataRowMPC.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowMPC.insert(1, x[1])
        if len(annotationDataRowMPC) == 2:
            annotationDataRowMPC.insert(0, id)
            finalAnnotationDataMPC.append(annotationDataRowMPC.copy())
            annotationDataRowMPC.clear()

    def annotateByMutpred(self, x, annotationDataRowMutpred, finalAnnotationDataMutpred, id):
        if x[0] == 'rankscore':
            annotationDataRowMutpred.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowMutpred.insert(1, x[1])
        if len(annotationDataRowMutpred) == 2:
            annotationDataRowMutpred.insert(0, id)
            finalAnnotationDataMutpred.append(annotationDataRowMutpred.copy())
            annotationDataRowMutpred.clear()

    def annotateByMtaster(self, x, annotationDataRowMtaster, finalAnnotationDataMtaster, id):
        if x[0] == 'pred':
            if x[1] == 'N':
                annotationDataRowMtaster.insert(2, 'Tolerant')
            else:
                annotationDataRowMtaster.insert(2, 'Damaging')
        if x[0] == 'converted_rankscore':
            annotationDataRowMtaster.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowMtaster.insert(1, x[1])
        if len(annotationDataRowMtaster) == 3:
            annotationDataRowMtaster.insert(0, id)
            finalAnnotationDataMtaster.append(annotationDataRowMtaster.copy())
            annotationDataRowMtaster.clear()  # annotationDataRowMtaster

    def annotateByMassessor(self, x, annotationDataRowMassessor, finalAnnotationDataMassessor, id):
        if x[0] == 'pred':
            if x[1] == 'N':
                annotationDataRowMassessor.insert(2, 'Tolerant')
            else:
                annotationDataRowMassessor.insert(2, 'Damaging')
        if x[0] == 'rankscore':
            annotationDataRowMassessor.insert(0, x[1])
        if x[0] == 'score':
            annotationDataRowMassessor.insert(1, x[1])
        if len(annotationDataRowMassessor) == 3:
            annotationDataRowMassessor.insert(0, id)
            finalAnnotationDataMassessor.append(annotationDataRowMassessor.copy())
            annotationDataRowMassessor.clear()

    def annotateByMRNN(self, x, annotationDatRowMRNN, finalAnnotationDataMRNN, id):
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDatRowMRNN.insert(2, 'Tolerant')
            else:
                annotationDatRowMRNN.insert(2, 'Damaging')
        if x[0] == 'rankscore':
            annotationDatRowMRNN.insert(0, x[1])
        if x[0] == 'score':
            annotationDatRowMRNN.insert(1, x[1])
        if len(annotationDatRowMRNN) == 3:
            annotationDatRowMRNN.insert(0, id)
            finalAnnotationDataMRNN.append(annotationDatRowMRNN.copy())
            annotationDatRowMRNN.clear()

    def annotateByMSVM(self, x, annotationDatRowMSVM, finalAnnotationDataMSVM, id):
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDatRowMSVM.insert(2, 'Tolerant')
            else:
                annotationDatRowMSVM.insert(2, 'Damaging')
        if x[0] == 'rankscore':
            annotationDatRowMSVM.insert(0, x[1])
        if x[0] == 'score':
            annotationDatRowMSVM.insert(1, x[1])
        if len(annotationDatRowMSVM) == 3:
            annotationDatRowMSVM.insert(0, id)
            finalAnnotationDataMSVM.append(annotationDatRowMSVM.copy())
            annotationDatRowMSVM.clear() #annotationDataRowMLR
    def annotateByMLR(self, x, annotationDatRowMLR, finalAnnotationDataMLR, id):
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDatRowMLR.insert(2, 'Tolerant')
            else:
                annotationDatRowMLR.insert(2, 'Damaging')
        if x[0] == 'rankscore':
            annotationDatRowMLR.insert(0, x[1])
        if x[0] == 'score':
            annotationDatRowMLR.insert(1, x[1])
        if len(annotationDatRowMLR) == 3:
            annotationDatRowMLR.insert(0, id)
            finalAnnotationDataMLR.append(annotationDatRowMLR.copy())
            annotationDatRowMLR.clear()

    def annotateByMCAP(self, x, annotationDatRowMCAP, finalAnnotationDataMCAP, id):
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDatRowMCAP.insert(2, 'Tolerant')
            else:
                annotationDatRowMCAP.insert(2, 'Damaging')
        if x[0] == 'rankscore':
            annotationDatRowMCAP.insert(0, x[1])
        if x[0] == 'score':
            annotationDatRowMCAP.insert(1, x[1])
        if len(annotationDatRowMCAP) == 3:
            annotationDatRowMCAP.insert(0, id)
            finalAnnotationDataMCAP.append(annotationDatRowMCAP.copy())
            annotationDatRowMCAP.clear()
    def annotateByLS2(self, x, annotationDatRowLS2, finalAnnotationDataLS2, id):
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDatRowLS2.insert(2, 'Tolerant')
            else:
                annotationDatRowLS2.insert(2, 'Damaging')
        if x[0] == 'rankscore':
            annotationDatRowLS2.insert(0, x[1])
        if x[0] == 'score':
            annotationDatRowLS2.insert(1, x[1])
        if len(annotationDatRowLS2) == 3:
            annotationDatRowLS2.insert(0, id)
            finalAnnotationDataLS2.append(annotationDatRowLS2.copy())
            annotationDatRowLS2.clear()
    def annotateByFathmm(self, x, annotationDatRowFathmm, finalAnnotationDataFathmm, id):
        if x[0] == 'pred':
            if x[1] == 'T':
                annotationDatRowFathmm.insert(2, 'Tolerant')
            else:
                annotationDatRowFathmm.insert(2, 'Damaging')
        if x[0] == 'converted_rankscore':
            annotationDatRowFathmm.insert(0, x[1])
        if x[0] == 'score':
            annotationDatRowFathmm.insert(1, x[1])
        if len(annotationDatRowFathmm) == 3:
            annotationDatRowFathmm.insert(0, id)
            finalAnnotationDataFathmm.append(annotationDatRowFathmm.copy())
            annotationDatRowFathmm.clear()
    def annotateByFXF(self, x, annotationDatRowFXF, finalAnnotationDataFXF, id):
        if x[0] == 'coding_pred':
            if x[1] == 'N':
                annotationDatRowFXF.insert(2, 'Tolerant')
            else:
                annotationDatRowFXF.insert(2, 'Damaging')
        if x[0] == 'coding_rankscore':
            annotationDatRowFXF.insert(0, x[1])
        if x[0] == 'coding_score':
            annotationDatRowFXF.insert(1, x[1])
        if len(annotationDatRowFXF) == 3:
            annotationDatRowFXF.insert(0, id)
            finalAnnotationDataFXF.append(annotationDatRowFXF.copy())
            annotationDatRowFXF.clear()
    def annotateByFMKL(self, x, annotationDatRowFMKL, finalAnnotationDataFMKL, id):
        if x[0] == 'coding_pred':
            if x[1] == 'N':
                annotationDatRowFMKL.insert(2, 'Tolerant')
            else:
                annotationDatRowFMKL.insert(2, 'Damaging')
        if x[0] == 'coding_rankscore':
            annotationDatRowFMKL.insert(0, x[1])
        if x[0] == 'coding_score':
            annotationDatRowFMKL.insert(1, x[1])
        if len(annotationDatRowFMKL) == 3:
            annotationDatRowFMKL.insert(0, id)
            finalAnnotationDataFMKL.append(annotationDatRowFMKL.copy())
            annotationDatRowFMKL.clear()
    def annotateByBdeladdaf(self, x, annotationDataRowBdeladdf, finalAnnotationDataBdeladdf, id):
        if x[0] == 'add_af':
            if x[1] == 'pred':
                if x[2] == 'T':
                    annotationDataRowBdeladdf.insert(2, 'Tolerant')
                else:
                    annotationDataRowBdeladdf.insert(2, 'Damaging')
            if x[1] == 'rankscore':
                annotationDataRowBdeladdf.insert(0, x[2])
            if x[1] == 'score':
                annotationDataRowBdeladdf.insert(1, x[2])
            if len(annotationDataRowBdeladdf) == 3:
                annotationDataRowBdeladdf.insert(0, id)
                finalAnnotationDataBdeladdf.append(annotationDataRowBdeladdf.copy())
                annotationDataRowBdeladdf.clear()
    def annotateByBdelnoaf(self, x, annotationDataRowBdelnoaf, finalAnnotationDataBdelnoaf, id):
        if x[0] == 'no_af':
            if x[1] == 'pred':
                if x[2] == 'T':
                    annotationDataRowBdelnoaf.insert(2, 'Tolerant')
                else:
                    annotationDataRowBdelnoaf.insert(2, 'Damaging')
            if x[1] == 'rankscore':
                annotationDataRowBdelnoaf.insert(0, x[2])
            if x[1] == 'score':
                annotationDataRowBdelnoaf.insert(1, x[2])
            if len(annotationDataRowBdelnoaf) == 3:
                annotationDataRowBdelnoaf.insert(0, id)
                finalAnnotationDataBdelnoaf.append(annotationDataRowBdelnoaf.copy())
                annotationDataRowBdelnoaf.clear()