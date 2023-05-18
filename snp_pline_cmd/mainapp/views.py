

from Bio import Entrez
import requests, sys, json
import grch38variantdata
from grch37variantdata import Grch37VariantData
from dbsnpvarientdataretrieval import DBSnpVarientDataRetrieval
from filewriting import FileWriting
from filedownloading import FileDownload
from excelrw import ExcelRW
from csvwriting import CSVFileWriting
from proteinvariantdata import ProteinVariantData
from refseqid_to_uniprotid import RefSeqId_to_UniProtId
#from variantinfoio import VarintInfoIO
import variantinfoio
import uuid
import myvariant
import numpy as np
# Create your views here.
noOk=0
rsIds_filename=''
grch37_filename=''
grch38_filename=''
prot_var_data_filename=''

varInfoIO=variantinfoio.VarintInfoIO()


def index(request):
    return render(request, 'index.html')
def runTools(request):

    return render(request, 'runtools.html')
def varDataRetrieval(request):
    return render(request, 'vardataretrieval.html')
results =''
def vardataretrievalprocessing():
    grch37B = False;
    grch38B = False;
    protDataB = False;
    rsIdsB = False
    varids=[]; givenTerm=''; totalSNPs='';
    csvfw = CSVFileWriting()
    global fName;
    fName = varInfoIO.variantNotFoundFile
    print('fName...', fName)

    if fName == '':
        givenTerm = input('Enter gene id')
        fabricatedTerm = givenTerm + ' AND missense variant[Function_Class]'
        Entrez.email = "usman.athar@gmail.com"
        handle = Entrez.esearch(db="snp", term=fabricatedTerm, retmax=5)
        variantData = Entrez.read(handle)
        totalSNPs = variantData['Count']
        varids = variantData["IdList"]
        chekboxValues = list(input(
            'Enter 1 for GRCh37 data , 2 for GRCh38 data, 3 for protein data,4 for variant ids, 5 for all options or enter multiple options separated by space\n').split())

        if chekboxValues.__contains__('5'):
            for x in range(4):
                chekboxValues.append(str(x + 1))
    else:
        #print('csvfw...', csvfw)
        givenTerm = varInfoIO.genId
        totalSNPs = varInfoIO.totalSNPs
        #newfName = 'media/' + fName
        newfName =  fName
        newfName = newfName.removesuffix('//')
        print('inside else givenTerm===', givenTerm)
        print('inside else totalSNPs===', totalSNPs)
        givenTerm=givenTerm.removesuffix('//')
        totalSNPs=totalSNPs.removesuffix('//')
        varData=csvfw.readCSVFile(newfName)
        for rsIds in varData:
            rsid=str(rsIds[0])
            rsid=rsid.removeprefix('rs')
            #print('rsid==',rsid)
            varids.append(rsid)
        #print('inside else...',varids)
        chekboxValues = list(input(
            'Enter 1 for GRCh37 data , 2 for GRCh38 data, 3 for protein data,4 for all options or enter multiple options separated by space\n').split())

        if chekboxValues.__contains__('4'):
            for x in range(3):
                chekboxValues.append(str(x + 1))

    print('chekboxValues=== ', chekboxValues)
    dbsnpvardata = DBSnpVarientDataRetrieval()
    rsIds = []
    grch37vardatainstance = Grch37VariantData()
    grch38vardatainstance = grch38variantdata.Grch38VariantData()

    csvfw=CSVFileWriting()
    prot_var_data = ProteinVariantData()
    refseqid_to_uniprotid = RefSeqId_to_UniProtId()
    count=0; protData=''
    for varid in varids:
        print('varid==', varid)
        dbsnpvardata_str = dbsnpvardata.getvariantdata(varid)
        #print('dbsnpvardata_str==', dbsnpvardata_str)
        if chekboxValues.__contains__('1'):  #writing GRCh37 chromosome coordinates
            grch37vardatainstance.parsevardatabystring(dbsnpvardata_str, varid)
        if chekboxValues.__contains__('2'): # writing GRCh38 chromosome coordinates in csv format
            grch38vardatainstance.parsevardatabystring(dbsnpvardata_str, varid)
        if chekboxValues.__contains__('3'):
            prot_var_data.parsevardatabystring(dbsnpvardata_str)
            if count==0:
                protData=refseqid_to_uniprotid.get_uniprotid_from_refseqid(givenTerm,prot_var_data.prot_var_data_dict['refSeqProtId'])

        if chekboxValues.__contains__('4'):
            id = 'rs' + str(varid)
            #query='http://www.mulinlab.org/vportal/portal/api?q=' + id + '&g=hg19&p=eur&f=pathogenicity'
            #query='https://myvariant.info/v1/variant/' + id

            rsIds.append(id)

        count = count + 1
    if chekboxValues.__contains__('1'):
        global grch37_filename
        grch37_filename = 'media/' + str(uuid.uuid4()) + '_GRCh37' + '.csv'
        csvfw.writeGRch37VarDataCSV(grch37vardatainstance.chr_coord_list, grch37_filename)
        print('grch37vardatainstance.chr_coord==', grch37vardatainstance.chr_coord_list)
        #excelrw.writeGRch37VarDataInExcel(grch37vardatainstance.chr_coord_list, grch37_filename)
        grch37B=True
    if chekboxValues.__contains__('2'):
        global grch38_filename
        #grch38_filename = 'media/' + str(uuid.uuid4()) + '_GRCh38' + '.xlsx'
        #excelrw.writeGRch38VarDataInExcel(grch38vardatainstance.chr_coord_dict, grch38_filename)
        grch38_filename = 'media/' + str(uuid.uuid4()) + '_GRCh38' + '.csv'
        csvfw.writeGRch38VarDataCSV(grch38vardatainstance.chr_coord_list, grch38_filename)

        grch38B=True
    if chekboxValues.__contains__('3'):
        global prot_var_data_filename
        prot_var_data_filename = 'media/' + str(uuid.uuid4()) +'_protein_data' + '.txt'
        fw = FileWriting()
        fw.writeProteinData(protData,prot_var_data.prot_var_data_dict,prot_var_data_filename)
        protDataB=True
    if chekboxValues.__contains__('4'):
        fw = FileWriting()
        global rsIds_filename
        rsIds_filename = 'media/' + str(uuid.uuid4()) + '_rsids' + '.txt'
        fw.writeFileFromList(rsIds, rsIds_filename)

        rsIdsB=True
    # fName == 'noFile':
        context = {

            'gene': givenTerm,
            'totalSNPs':totalSNPs,
            'grch37B':grch37B,
            'grch38B':grch38B,
            'rsIdsB':rsIdsB,
            'protDataB':protDataB

           }
        '''
   
        context = {

            'grch37B': grch37B,
            'grch38B': grch38B,
            'rsIdsB': rsIdsB,
            'protDataB': protDataB

        }
        '''


#annotate variants using myvariantinfo restful api

no_of_pathg_vars_fname=''
def annotateVariants():
    genId = input('provide genId')
    fabricatedTerm = genId + ' AND missense variant[Function_Class]'

    varInfoIO.parseMyVarinats(fabricatedTerm,genId)
    selectedTools=varInfoIO.chekboxValues
    print('selectedTools== ',selectedTools)

    tempvarfile=varInfoIO.variantNotFoundFile[6:]
    print('varInfoIO.variantNotFoundFile===', varInfoIO.variantNotFoundFile)
    print('tempvarfile===',tempvarfile)
    ### writing no of pathogenic and non pathogenic variants in csv
    no_of_pathg_vars_list=no_of_pahtogenic_variants(varInfoIO)
    print('no_of_pathg_vars_list==', no_of_pathg_vars_list)
    fieldnames = ['Tool','Pathogenic Variants','Non Pathogenic Variants']
    global no_of_pathg_vars_fname
    no_of_pathg_vars_fname = 'media/' + str(uuid.uuid4()) + '_no_of_pathg_vars' + '.csv'
    csvfw = CSVFileWriting()
    csvfw.writeAnnotationDateCSV(no_of_pathg_vars_list, no_of_pathg_vars_fname, fieldnames)
    optionForVarNotFound=''
    if varInfoIO.variantNotFoundFile!='':
        optionForVarNotFound=input('Please see '+"''"+ varInfoIO.variantNotFoundFile +"''"+ 'file to see variants not found by MyVariant Info. \n Would you like '
    'to get chromosome coordinates and proteomic data of these variants? press Y for Yes and n for No')
    match optionForVarNotFound:
        case 'y' | 'Y':
            vardataretrievalprocessing()
        case _:
           quit()
    ### end
    # writing no of pathogenic and non pathogenic variants in csv



def entrezIdSummaryInfo(givenId):
    print("givenId : ",givenId)
    Entrez.email = "usman.athar@gmail.com"
    handle = Entrez.esummary(db="snp", id=givenId, retmode="xml")
    records = Entrez.read(handle, validate=False)

    print("records  summary ", records)

def entrezIdfFetchInfo(givenId):
    print("givenId : ",givenId)
    Entrez.email = "usman.athar@gmail.com"
    handle = Entrez.efetch(db="snp", id=givenId)
    records = handle.read()
    print("records efetch ", records)

def ensembleIdInfo(givenId):
    print("givenId : ", givenId)
    server = "https://rest.ensembl.org"
    ext = "/variation/human/rs56116432?"

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        global noOk
        noOk+=1
     #   r.raise_for_status()
       # sys.exit()

    decoded = r.json()
    print("ensemble info===")
    print(repr(decoded))
    #print("decoded = " decoded)

def ensembleIdImpact(givenId):
    print("givenId : ", givenId)
    server = "https://rest.ensembl.org"
    ext = "/vep/human/id/" + givenId +"?"

    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    print("variant impact by VEP ", decoded)
def myvariantinfo(givenId):
    server = "https://myvariant.info/v1"
    ext = "/query/?q=" + givenId  # +rsIds[0]
    r = requests.get(server + ext)  # , headers={"Content-Type": "application/json"})
    decoded = r.json()
    print("myvarinat info===")
    print(repr(decoded))
def spdiserviceIdInfo(givenId):
    server = "https://api.ncbi.nlm.nih.gov/variation/v0/"
    ext="beta/refsnp/" + givenId
    r = requests.get(server + ext)  # , headers={"Content-Type": "application/json"})
    spdiIdInfo = r.json()
    return spdiIdInfo
    #print("variant info by spdi ", decoded)
def downloadRSIdFile(request):
    fDownload=FileDownload()
    return fDownload.download_file(request, '/'+ rsIds_filename)
def downloadGRCh37File(request):
    fDownload=FileDownload()
    #return fDownload.downloadExcelFile(request, '/'+  grch37_filename)
    return fDownload.download_file(request, '/' + grch37_filename)
def downloadGRCh38File(request):
    fDownload=FileDownload()
    #return fDownload.downloadExcelFile(request, '/'+  grch38_filename)
    return fDownload.download_file(request, '/' + grch38_filename)
def downloadProteinDataFile(request):
    fDownload=FileDownload()
    return fDownload.download_file(request, '/'+  prot_var_data_filename)

def downloadVarAnnotation(request, tool):
    global varInfoIO
    match tool:
        case 'sift':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.siftFile)
        case 'sift4g':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.sift4gFile)
        case 'provean':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.proveanFile)
        case 'mvp':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.mvpFile)
        case 'polyphen2hvar':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.polyphen2hvarFile)
        case 'polyphen2hdiv':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.polyphen2hdivFile)
        case 'primateai':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.primateaiFile)
        case 'revel':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.revelFile)
        case 'mpc':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.mpcFile)
        case 'mutpred':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.mutpredFile)
        case 'mtaster':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.mtasterFile)
        case 'massessor':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.massessorFile)
        case 'mrnn':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.mrnnFile)
        case 'msvm':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.msvmFile)
        case 'mrl':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.mlrFile)
        case 'mcap':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.mcapFile)
        case 'ls2':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.ls2File)
        case 'fathmm':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.fathmmFile)
        case 'fxf':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.fxfFile)
        case 'fkml':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.fmklFile)
        case 'bdeladdaf':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.bdeladdafFile)
        case 'bdelnoaf':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.bdelnoafFile)
        case 'vest4':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.vest4File)
        case 'dann':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.dannFile)
        case 'eigen':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.eigenFile)
        case 'eigenpc':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.eigenpcFile)
        case 'deogen2':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.doegen2File)
        case 'genocanyon':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.genocanyonFile)
        case 'genocanyon':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.genocanyonFile)
        case 'pathogenicvarsfound':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.pathogenicVarFile)
        case 'varsnotfound':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + varInfoIO.variantNotFoundFile)
        case 'table':
            fDownload = FileDownload()
            return fDownload.download_file(request, '/' + no_of_pathg_vars_fname)

        case _:
            print('No tool selected')
def selectedVarsDataRetrieval(request):
    #print('varnotfoundfile==',varFile)
    fName = request.POST.get('varNotFoundTxt')
    genId = request.POST.get('genId')
    totSNPs = request.POST.get('totsnps')
    print('fName==', fName )
    print('genId==', genId )
    print('totSNPS==', totSNPs)
    context ={
        'varNotFoundFile':fName,
        'genId': genId,
        'totSNPs': totSNPs,

    }
    return render(request, 'tools_data_retrieval.html',context)
def selVarDataRetrievalProcessing(request):
    if request.method == 'POST':
        fName = request.POST.get('varNotFoundTxt')
        fName='media/'+fName
        fName=fName.removesuffix('/')

        return vardataretrievalprocessing(request,fName)
    #return None
def no_of_pahtogenic_variants(varInfoIO):
    pathogenic_vars =0;non_pathogenic_vars=0;
    local_list=[]
    global_list=[]
    total_snps_int=int(varInfoIO.totalSNPs)
    print('type===',type(varInfoIO.siftPVars))
    if varInfoIO.isSift:
        print('inside isSift')
        non_pathogenic_vars=total_snps_int-varInfoIO.siftPVars
        local_list.append('Sift')
        local_list.append(varInfoIO.siftPVars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        print('local_list isSift',local_list)
        local_list.clear()
    if varInfoIO.isSift4G:
        non_pathogenic_vars = total_snps_int - varInfoIO.sift4gPVars
        local_list.append('Sift4G')
        local_list.append(varInfoIO.sift4gPVars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isProvean:
        non_pathogenic_vars = total_snps_int - varInfoIO.proveanPVars
        local_list.append('Provean')
        local_list.append(varInfoIO.proveanPVars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isMVP:
        non_pathogenic_vars = total_snps_int - varInfoIO.mvpPVars
        local_list.append('MVP')
        local_list.append(varInfoIO.mvpPVars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isPolyphen2HVAR:
        pathogenic_vars=varInfoIO.polyphen2hvarPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('Polyphen2_HVAR')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isPolyphen2HDIV:
        pathogenic_vars = varInfoIO.polyphen2hdivPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('Polyphen2_HDIV')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isPrimateAI:
        pathogenic_vars = varInfoIO.primateaiPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('PrimateAI')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isRevel:
        pathogenic_vars = varInfoIO.revelPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('REVEL')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isMPC:
        pathogenic_vars = varInfoIO.mpcPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('MPC')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isMutPred:
        pathogenic_vars = varInfoIO.mutpredPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('MutPred')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isMTaster:
        pathogenic_vars = varInfoIO.mtasterPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('MutationTaster')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isMAssessor:
        pathogenic_vars = varInfoIO.massessorPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('MutationAssessor')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isMRNN:
        pathogenic_vars = varInfoIO.mrnnPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('MetaRNN')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isMSVM:
        pathogenic_vars = varInfoIO.msvmPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('MetaSVM')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isMLR:
        pathogenic_vars = varInfoIO.mlrPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('MetaLR')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isMCAP:
        pathogenic_vars = varInfoIO.mcapPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('M-Cap')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isLS2:
        pathogenic_vars = varInfoIO.ls2PVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('LIST-S2')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isFATHMM:
        pathogenic_vars = varInfoIO.fathmmPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('FATHMM')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isFXF:
        pathogenic_vars = varInfoIO.fxfPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('fathmm-XF')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isFMKL:
        pathogenic_vars = varInfoIO.fmklPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('fathmm-MKL')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isBDelAddAF:
        pathogenic_vars = varInfoIO.bdeladdafPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('BayesDel-addAF')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isBDelNoAF:
        pathogenic_vars = varInfoIO.bdelnoafPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('BayesDel-noAF')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isVest4:
        pathogenic_vars = varInfoIO.vest4PVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('VEST4')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isDANN:
        pathogenic_vars = varInfoIO.dannPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('DANN')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isEigen:
        pathogenic_vars = varInfoIO.eigenPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('Eigen')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isEigenPC:
        pathogenic_vars = varInfoIO.eigenpcPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('Eigen-PC')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isDeogen2:
        pathogenic_vars = varInfoIO.doegen2PVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('DEOGEN2')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    if varInfoIO.isGenoCanyon:
        pathogenic_vars = varInfoIO.genocanyonPVars
        non_pathogenic_vars = total_snps_int - pathogenic_vars
        local_list.append('GenoCanyon')
        local_list.append(pathogenic_vars)
        local_list.append(non_pathogenic_vars)
        global_list.append(local_list.copy())
        local_list.clear()
    print('global_list== ',global_list)
    return global_list
