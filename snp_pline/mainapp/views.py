
from django.shortcuts import render
from Bio import Entrez
import requests, sys, json
from . import grch38variantdata
from .grch37variantdata import Grch37VariantData
from .dbsnpvarientdataretrieval import DBSnpVarientDataRetrieval
from.filewriting import FileWriting
from .filedownloading import FileDownload
import numpy as np
# Create your views here.
global noOk
def index(request):
    return render(request, 'index.html')
def routeRequest(request):
    if request.method == 'POST':
        whattodo=str(request.POST.get('actiontyperdbtn', 'run'))
        givenTerm=request.POST.get('genId')

        fabricatedTerm=givenTerm + ' AND missense variant[Function_Class]'
        print("term=== ", givenTerm)
        if whattodo=='run':
            Entrez.email = "usman.athar@gmail.com"

        else:
           Entrez.email = "usman.athar@gmail.com"
           handle = Entrez.esearch(db="snp", term=fabricatedTerm, retmax=10)
           variantData = Entrez.read(handle)

           totalSNPs=variantData['Count']
           ids=variantData["IdList"]
           #print("ids", ids)
           rsIds=[]
           string = 'rs'
           grch37vardatainstance = Grch37VariantData()
           for id in ids:
               dbsnpvardata = DBSnpVarientDataRetrieval()
               dbsnpvardata_str = dbsnpvardata.getvariantdata(id)
               # print('dbsnpvardata_str  ',dbsnpvardata_str)
               # grch38vardatainstance= grch38variantdata.Grch38VariantData()
               # grch38vardata=grch38vardatainstance.parsevardatabystring(dbsnpvardata_str)

               id = 'rs' + str(id)
               rsIds.append(id)
               grch37vardatainstance.parsevardatabystring(dbsnpvardata_str, id)
           print('grch37vardatainstance.chr_coord_dict===', grch37vardatainstance.chr_coord_dict)
           #print("rsIds", rsIds)
           fw=FileWriting()
           fw.writeFileFromList(rsIds)
           context = {
               'rsids': rsIds,
               'gene': givenTerm,
               'totalSNPs':totalSNPs,
           }


           #spdiIdInfo=spdiserviceIdInfo('4537')

           #spdiIdInfo_dict=json.loads(spdiIdInfo_str)
           #spdiIdInfo_dict=spdiIdInfo_dict.get('primary_snapshot_data')
           '''
           spdiIdInfo_str=spdiIdInfo_str[grch37Index:]
           print('spdiIdInfo_str after 37== ', spdiIdInfo_str)
           triplebracesindex=spdiIdInfo_str.find('}]}')
           varData37 = spdiIdInfo_str[:triplebracesindex]
           print('varData37== ', varData37)

           
           is38=False
           dnavarinfo37=''
           orientation=0
           vals = spdiIdInfo_dict.values()

           for value in vals:
                val=str(value)
                print('vals == \n', val)
                if val.find('38')>0:
                    is38=True
                    print('is38=== ',is38)
                elif val.find('37')>0:
                    is38=False
                if is38:
                    if str(value).__contains__('false  '):
                        orientation=1
                    if str(value).startswith('NC_') and str(value).__contains__('>'):
                        dnavarinfo38=value
                    if str(value).startswith('NP_') and str(value).endswith(''):
                        proteinvarinfo=value
                else:
                    if str(value).__contains__('false'):
                        orientation=1
                    if str(value).startswith('NC_') and str(value).__contains__('>'):
                        dnavarinfo37=value
                    if str(value).startswith('NP_') and str(value).endswith(''):
                        proteinvarinfo=value
           siftparams=str(dnavarinfo37)+str(orientation)
           print('siftparams== ',siftparams)
'''

           #entrezIdSummaryInfo(ids[0])
           #entrezIdfFetchInfo(ids[0])
           #for rsId in rsIds:
           #ensembleIdInfo('4537')

           #ensembleIdImpact(ids[0])
           #myvariantinfo('4537')


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
def downloadFile(request):
    fDownload=cpf.FileDownload()
    return fDownload.download_file(request, '/'+ acFile)