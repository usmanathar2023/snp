
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
        csvfw=CSVFileWriting()
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
                #query='http://www.mulinlab.org/vportal/portal/api?q=' + id + '&g=hg19&p=eur&f=pathogenicity'
                #query='https://myvariant.info/v1/variant/' + id

                rsIds.append(id)

            count = count + 1
        if chekboxValues.__contains__('grch37'):
            global grch37_filename
            grch37_filename = 'media/' + str(uuid.uuid4()) + '_GRCh37' + '.csv'
            csvfw.writeGRch37VarDataCSV(grch37vardatainstance.chr_coord_list, grch37_filename)
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
    context = {

        'gene': givenTerm,
        'totalSNPs':totalSNPs,
        'grch37B':grch37B,
        'grch38B':grch38B,
        'rsIdsB':rsIdsB,
        'protDataB':protDataB

       }
    return render(request, 'rsids.html', context)

#annotate variants using myvariantinfo restful api
def annotateVariants(request):
    mv = myvariant.MyVariantInfo()
    mvVar = ''
    if request.method == 'POST':
        givenTerm = request.POST.get('genId')
        fabricatedTerm = givenTerm + ' AND missense variant[Function_Class]'
        Entrez.email = "usman.athar@gmail.com"
        handle = Entrez.esearch(db="snp", term=fabricatedTerm, retmax=5)
        variantData = Entrez.read(handle)
        totalSNPs = variantData['Count']
        varids = variantData["IdList"]
        chekboxValues = str(request.POST.getlist('chekboxValues', ''))
        rsIds = []
        var=''
        print('chekboxValues=== ', chekboxValues)
        for varid in varids:
            id = 'rs' + str(varid)
            var=mv.getvariant(id)

            if chekboxValues.__contains__('sift'):  #writing GRCh37 chromosome coordinates
                #mv.getvariant(id, fields=['cadd.sift.cat', 'cadd.sift.val'])
                #fld=mv.get_fields('cadd.sift.cat')
                print('fld=== ', mv.getvariant(id, fields=['dbnsfp.cadd,dbnsfp.cadd.pred', 'dbnsfp.cadd.raw_rankscore', 'dbnsfp.cadd.raw_score',
                                                           'dbnsfp.sift,dbnsfp.sift.pred', 'dbnsfp.sift.converted_rankscore', 'dbnsfp.sift.score',
                                                           'dbnsfp.sift4g,dbnsfp.sift4g.pred', 'dbnsfp.sift4g.converted_rankscore', 'dbnsfp.sift4g.score',
                                                           'dbnsfp.provean,dbnsfp.provean.pred','dbnsfp.provean.rankscore', 'dbnsfp.provean.score',
                                                           'dbnsfp.primateai,dbnsfp.primateai.pred','dbnsfp.primateai.rankscore', 'dbnsfp.primateai.score',
                                                           'dbnsfp.polyphen2.hdiv,dbnsfp.polyphen2.hdiv.pred','dbnsfp.polyphen2.hdiv.rankscore', 'dbnsfp.polyphen2.hdiv.score',
                                                           'dbnsfp.polyphen2.hvar', 'dbnsfp.polyphen2.hvar.pred', 'dbnsfp.polyphen2.hvar.rankscore', 'dbnsfp.polyphen2.hvar.score',
                                                           'dbnsfp.revel','dbnsfp.revel.rankscore', 'dbnsfp.revel.score',
                                                           'dbnsfp.mvp', 'dbnsfp.mvp.rankscore', 'dbnsfp.mvp.score',
                                                           'dbnsfp.mpc', 'dbnsfp.mpc.rankscore', 'dbnsfp.mpc.score',
                                                           'dbnsfp.mutpred,dbnsfp.mutpred.pred','dbnsfp.mutpred.rankscore', 'dbnsfp.mutpred.score',
                                                           'dbnsfp.mutationtaster,dbnsfp.mutationtaster.pred','dbnsfp.mutationtaster.converted_rankscore', 'dbnsfp.mutationtaster.score',
                                                           'dbnsfp.mutationassessor,dbnsfp.mutationassessor.pred','dbnsfp.mutationassessor.rankscore', 'dbnsfp.mutationassessor.score',
                                                           'dbnsfp.metarnn,dbnsfp.metarnn.pred','dbnsfp.metarnn.rankscore', 'dbnsfp.metarnn.score',
                                                           'dbnsfp.metasvm,dbnsfp.metasvm.pred','dbnsfp.metasvm.rankscore', 'dbnsfp.metasvm.score',
                                                           'dbnsfp.metalr,dbnsfp.metalr.pred','dbnsfp.metalr.rankscore', 'dbnsfp.metalr.score',
                                                           'dbnsfp.m-cap,dbnsfp.m-cap.pred','dbnsfp.m-cap.rankscore', 'dbnsfp.m-cap.score',
                                                           'dbnsfp.lrt,dbnsfp.lrt.pred','dbnsfp.lrt.converted_rankscore', 'dbnsfp.lrt.score',
                                                           'dbnsfp.list-s2,dbnsfp.list-s2.pred','dbnsfp.list-s2.rankscore', 'dbnsfp.list-s2.score',
                                                           'dbnsfp.linsight','dbnsfp.linsight.rankscore', 'dbnsfp.linsight.score',
                                                           'dbnsfp.gerp++,','dbnsfp.gerp++.rs', 'dbnsfp.gerp++.score',
                                                           'dbnsfp.fathmm,dbnsfp.fathmm.pred','dbnsfp.fathmm.rankscore','dbnsfp.fathmm.score',
                                                           'dbnsfp.fathmm-mkl,dbnsfp.fathmm-mkl.coding_pred','dbnsfp.fathmm-mkl.coding_rankscore', 'dbnsfp.fathmm-mkl.coding_score',
                                                           'dbnsfp.fathmm-xf,dbnsfp.fathmm-xf.coding_pred','dbnsfp.fathmm-xf.coding_rankscore', 'dbnsfp.fathmm-xf.coding_score',]))

            if chekboxValues.__contains__('polyphen'): # writing GRCh38 chromosome coordinates in csv format
                mv.getvariant(id, fields=['cadd.sift.cat', 'cadd.sift.val'])
            if chekboxValues.__contains__('provean'):
                mv.getvariant(id, fields=['cadd.sift.cat', 'cadd.sift.val'])
            if chekboxValues.__contains__('cadd'):
                id = 'rs' + str(varid)

                mv = myvariant.MyVariantInfo()
                mvVar=mv.getvariant(id,fields=['cadd.sift.cat','cadd.sift.val'])
                #siftFields=mvVar.get_fields("sift")

                print('mvVar=== ',mvVar)




            rsIdsB=True
    context = {

        'gene': givenTerm,
        'totalSNPs':totalSNPs,



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
    #return fDownload.downloadExcelFile(request, '/'+  grch37_filename)
    return fDownload.download_file(request, '/' + grch37_filename)
def downloadGRCh38File(request):
    fDownload=FileDownload()
    #return fDownload.downloadExcelFile(request, '/'+  grch38_filename)
    return fDownload.download_file(request, '/' + grch38_filename)
def downloadProteinDataFile(request):
    fDownload=FileDownload()
    return fDownload.download_file(request, '/'+  prot_var_data_filename)




