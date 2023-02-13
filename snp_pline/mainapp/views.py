
from django.shortcuts import render
from Bio import Entrez
import requests, sys, json
from . import grch38variantdata
from .grch37variantdata import Grch37VariantData
from .dbsnpvarientdataretrieval import DBSnpVarientDataRetrieval
from.filewriting import FileWriting
from .filedownloading import FileDownload
from .excelrw import ExcelRW
from .csvwriting import CSVFileWriting
from .proteinvariantdata import ProteinVariantData
from .refseqid_to_uniprotid import RefSeqId_to_UniProtId
from .variantinfoio import VarintInfoIO
import uuid
import myvariant
import numpy as np
# Create your views here.
noOk=0
rsIds_filename=''
grch37_filename=''
grch38_filename=''
prot_var_data_filename=''


def index(request):
    return render(request, 'index.html')
def runTools(request):

    return render(request, 'runtools.html')
def varDataRetrieval(request):
    return render(request, 'vardataretrieval.html')
results =''
def vardataretrievalprocessing(request):
    grch37B = False;
    grch38B = False;
    protDataB = False;
    rsIdsB = False
    varids=[]; givenTerm=''; totalSNPs='';
    csvfw = CSVFileWriting()
    if request.method == 'POST':
        fName= request.POST.get('varNotFoundTxt')
        print('fName...', fName)

        if fName == 'noFile':
            givenTerm = request.POST.get('genId')
            fabricatedTerm = givenTerm + ' AND missense variant[Function_Class]'
            Entrez.email = "usman.athar@gmail.com"
            handle = Entrez.esearch(db="snp", term=fabricatedTerm, retmax=1)
            variantData = Entrez.read(handle)
            totalSNPs = variantData['Count']
            varids = variantData["IdList"]
        else:
            #print('csvfw...', csvfw)
            givenTerm = str(request.POST.get('genid'))
            totalSNPs = str(request.POST.get('totsnps'))
            newfName = 'media/' + fName
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
        chekboxValues = str(request.POST.getlist('chekboxValues', ''))
        print('chekboxValues=== ', chekboxValues)
        dbsnpvardata = DBSnpVarientDataRetrieval()
        rsIds = []
        grch37vardatainstance = Grch37VariantData()
        grch38vardatainstance = grch38variantdata.Grch38VariantData()
        excelrw = ExcelRW()
        csvfw=CSVFileWriting()
        prot_var_data = ProteinVariantData()
        refseqid_to_uniprotid = RefSeqId_to_UniProtId()
        count=0; protData=''
        for varid in varids:
            print('varid==', varid)
            dbsnpvardata_str = dbsnpvardata.getvariantdata(varid)
            #print('dbsnpvardata_str==', dbsnpvardata_str)
            if chekboxValues.__contains__('grch37'):  #writing GRCh37 chromosome coordinates
                grch37vardatainstance.parsevardatabystring(dbsnpvardata_str, varid)
            if chekboxValues.__contains__('grch38'): # writing GRCh38 chromosome coordinates in csv format
                grch38vardatainstance.parsevardatabystring(dbsnpvardata_str, varid)
            if chekboxValues.__contains__('protdata'):
                prot_var_data.parsevardatabystring(dbsnpvardata_str)
                if count==0:
                    protData=refseqid_to_uniprotid.get_uniprotid_from_refseqid(givenTerm,prot_var_data.prot_var_data_dict['refSeqProtId'])

            if chekboxValues.__contains__('rsids'):
                id = 'rs' + str(varid)
                #query='http://www.mulinlab.org/vportal/portal/api?q=' + id + '&g=hg19&p=eur&f=pathogenicity'
                #query='https://myvariant.info/v1/variant/' + id

                rsIds.append(id)

            count = count + 1
        if chekboxValues.__contains__('grch37'):
            global grch37_filename
            grch37_filename = 'media/' + str(uuid.uuid4()) + '_GRCh37' + '.csv'
            csvfw.writeGRch37VarDataCSV(grch37vardatainstance.chr_coord_list, grch37_filename)
            print('grch37vardatainstance.chr_coord==', grch37vardatainstance.chr_coord_list)
            #excelrw.writeGRch37VarDataInExcel(grch37vardatainstance.chr_coord_list, grch37_filename)

            grch37B=True
        if chekboxValues.__contains__('grch38'):
            global grch38_filename
            #grch38_filename = 'media/' + str(uuid.uuid4()) + '_GRCh38' + '.xlsx'
            #excelrw.writeGRch38VarDataInExcel(grch38vardatainstance.chr_coord_dict, grch38_filename)
            grch38_filename = 'media/' + str(uuid.uuid4()) + '_GRCh38' + '.csv'
            csvfw.writeGRch38VarDataCSV(grch38vardatainstance.chr_coord_list, grch38_filename)

            grch38B=True
        if chekboxValues.__contains__('protdata'):
            global prot_var_data_filename
            prot_var_data_filename = 'media/' + str(uuid.uuid4()) +'_protein_data' + '.txt'
            fw = FileWriting()
            fw.writeProteinData(protData,prot_var_data.prot_var_data_dict,prot_var_data_filename)
            protDataB=True
        if chekboxValues.__contains__('rsids'):
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
    return render(request, 'rsids.html', context)

#annotate variants using myvariantinfo restful api
varInfoIO='';
def annotateVariants(request):

    if request.method == 'POST':
        givenTerm = request.POST.get('genId')
        fabricatedTerm = givenTerm + ' AND missense variant[Function_Class]'
        global varInfoIO
        varInfoIO=VarintInfoIO()
        varInfoIO.parseMyVarinats(request,fabricatedTerm)
        selectedTools=varInfoIO.chekboxValues
        print('selectedTools== ',selectedTools)

    tempvarfile=varInfoIO.variantNotFoundFile[6:]
    print('varInfoIO.variantNotFoundFile===', varInfoIO.variantNotFoundFile)
    print('tempvarfile===',tempvarfile)
    context = {

        'gene': givenTerm,
        'totalSNPs': varInfoIO.totalSNPs,
        'sift':varInfoIO.isSift,
        'sift4g': varInfoIO.isSift4G,
        'provean': varInfoIO.isProvean,
        'mvp':varInfoIO.isMVP,
        'polyphen2HVAR':varInfoIO.isPolyphen2HVAR,
        'polyphen2HDIV':varInfoIO.isPolyphen2HDIV,
        'primateAI':varInfoIO.isPrimateAI,
        'revel':varInfoIO.isRevel,
        'mpc':varInfoIO.isMPC,
        'mutpred':varInfoIO.isMutPred,
        'mtaster':varInfoIO.isMTaster,
        'massessor':varInfoIO.isMAssessor,
        'mrnn':varInfoIO.isMRNN,
        'msvm':varInfoIO.isMSVM,
        'mlr':varInfoIO.isMLR,
        'mcap':varInfoIO.isMCAP,
        'ls2':varInfoIO.isLS2,
        'fathmm':varInfoIO.isFATHMM,
        'fxf':varInfoIO.isFXF,
        'fkml':varInfoIO.isFMKL,
        'bdeladdaf':varInfoIO.isBDelAddAF,
        'bdelnoaf':varInfoIO.isBDelNoAF,
        'vest4':varInfoIO.isVest4,
        'dann':varInfoIO.isDANN,
        'eigen':varInfoIO.isEigen,
        'eigenpc':varInfoIO.isEigenPC,
        'deogen2':varInfoIO.isDeogen2,
        'genocanyon':varInfoIO.isGenoCanyon,
        'pathogenicVarsLen':len(varInfoIO.pathogenicVariants),
        'varsNotFoundLen': len(varInfoIO.varDataNotFound),
        'varNotFoundFile':tempvarfile,
                             #varInfoIO.variantNotFoundFile,

    }
    return render(request, 'varAnnotationDownloand.html', context)


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


