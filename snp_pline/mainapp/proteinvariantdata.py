import requests, sys, json
import re
import ast


class ProteinVariantData:
    def __init__(self):
        self.prot_coord_list = []
        self.prot_var_data_dict={}

    def parsevardatabystring(self, varData):
        # getting orientation
        chr_coord_list = []
        index=0; index2=0;protein_id='';position=0;original_allele='';mutant_allele='';
        index = varData.find('{"allele": {"spdi": {"seq_id": "NP_')
        while index > -1:
            #getting string having protein variant data information
            varData = varData[index:]
            index2 = varData.find('"}]}]')
            varData=varData[:index2+5]
            index=varData.find('NP_')
            index2=varData.find(',')
            protein_id=varData[index:index2-1]
            varData=varData[index2+1:]
            index=varData.find(':')
            index2 = varData.find(',')
            position=varData[index+2:index2]
            position=int(position)+1
            varData = varData[index2+1:]
            index = varData.find(':')
            index2 = varData.find(',')
            original_allele = varData[index+3:index2-1]
            varData = varData[index2+1:]
            index = varData.find(':')
            index2 = varData.find('"}')
            mutant_allele = varData[index+3:index2]
            index = varData.find('{"allele": {"spdi": {"seq_id": "NP_')
            if(original_allele==mutant_allele):
                continue
            prot_coord=original_allele+str(position)+ mutant_allele
            print('prot_coord== ', prot_coord)
            self. prot_coord_list.append(prot_coord)
            #print("varData last== ", varData)
        self.prot_var_data_dict['protCoord']=self.prot_coord_list
        self.prot_var_data_dict['refSeqProtId']=protein_id
            # return chrcooridnates38


