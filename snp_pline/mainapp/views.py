
from django.shortcuts import render
from Bio import Entrez
import requests, sys, json
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
           handle = Entrez.esearch(db="snp", term=fabricatedTerm, retmax=10000)
           variantData = Entrez.read(handle)
           totalSNPs=variantData['Count']
           ids=variantData["IdList"]
           #print("ids", ids)
           rsIds=[]
           string = 'rs'
           for id in ids:
               id='rs' + str(id)
               rsIds.append(id)
           #print("rsIds", rsIds)
           context = {
               'rsids': rsIds,
               'gene': givenTerm,
               'totalSNPs':totalSNPs,
           }
           spdiserviceIdInfo('4537')
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
