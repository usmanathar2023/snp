
from django.shortcuts import render
from Bio import Entrez
import requests, sys, json
from . import grch38variantdata
from .grch37variantdata import Grch37VariantData
from .dbsnpvarientdataretrieval import DBSnpVarientDataRetrieval
from.filewriting import FileWriting
from .filedownloading import FileDownload
from .excelrw import ExcelRW
from .proteinvariantdata import ProteinVariantData
from .refseqid_to_uniprotid import RefSeqId_to_UniProtId
import uuid
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
    if request.method == 'POST':
        givenTerm = request.POST.get('genId')
        fabricatedTerm = givenTerm + ' AND missense variant[Function_Class]'
        Entrez.email = "usman.athar@gmail.com"
        handle = Entrez.esearch(db="snp", term=fabricatedTerm, retmax=5)
        variantData = Entrez.read(handle)
        totalSNPs = variantData['Count']
        varids = variantData["IdList"]
        chekboxValues = str(request.POST.getlist('chekboxValues', ''))
        print('chekboxValues=== ', chekboxValues)
        dbsnpvardata = DBSnpVarientDataRetrieval()
        rsIds = []
        grch37vardatainstance = Grch37VariantData()
        grch38vardatainstance = grch38variantdata.Grch38VariantData()
        excelrw = ExcelRW()
        prot_var_data = ProteinVariantData()
        refseqid_to_uniprotid = RefSeqId_to_UniProtId()
        count=0; protData=''
        for varid in varids:
            dbsnpvardata_str = dbsnpvardata.getvariantdata(varid)
            if chekboxValues.__contains__('grch37'):  #writing GRCh37 chromosome coordinates
                grch37vardatainstance.parsevardatabystring(dbsnpvardata_str, varid)
            if chekboxValues.__contains__('grch38'): # writing GRCh38 chromosome coordinates in csv format
                grch38vardatainstance.parsevardatabystring(dbsnpvardata_str, varid)
            if chekboxValues.__contains__('protdata'):
                if count==0:
                    protData=refseqid_to_uniprotid.get_uniprotid_from_refseqid(givenTerm)
                prot_var_data.parsevardatabystring(dbsnpvardata_str)
            if chekboxValues.__contains__('rsids'):
                id = 'rs' + str(varid)
                rsIds.append(id)

            count = count + 1
        if chekboxValues.__contains__('grch37'):
            global grch37_filename
            grch37_filename = 'media/' + str(uuid.uuid4()) + '_GRCh37' + '.xlsx'
            excelrw.writeGRch37VarDataInExcel(grch37vardatainstance.chr_coord_dict, grch37_filename)

            grch37B=True
        if chekboxValues.__contains__('grch38'):
            global grch38_filename
            grch38_filename = 'media/' + str(uuid.uuid4()) + '_GRCh38' + '.xlsx'
            excelrw.writeGRch38VarDataInExcel(grch38vardatainstance.chr_coord_dict, grch38_filename)

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
    context = {

        'gene': givenTerm,
        'totalSNPs':totalSNPs,
        'grch37B':grch37B,
        'grch38B':grch38B,
        'rsIdsB':rsIdsB,
        'protDataB':protDataB

       }
    return render(request, 'rsids.html', context)

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
    return fDownload.downloadExcelFile(request, '/'+  grch37_filename)
def downloadGRCh38File(request):
    fDownload=FileDownload()
    return fDownload.downloadExcelFile(request, '/'+  grch38_filename)
def downloadProteinDataFile(request):
    fDownload=FileDownload()
    return fDownload.download_file(request, '/'+  prot_var_data_filename)




