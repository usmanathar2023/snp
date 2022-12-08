import requests, sys, json
import re
import ast


class Grch38VariantData:
    def __init__(self):
        self.chr_coord_dict = {}

    def parsevardatabystring(self, varData, var_id):
        # getting orientation
        chr_coord_list = []
        index=0; index2=0;protein_id='';position=0;original_allele='';mutant_allele='';
        index = varData.find('{"allele":{"spdi":{"seq_id":"NP_')

        while index > -1:
            #getting string having protein variant data information
            index2 = varData.find('}}')
            varData=varData[index:index2]
            index=varData.find('NP_')
            index2=varData.find(',')
            protein_id=varData[index:index2]
            varData=varData[index2:]
            index=varData.find(':')
            index2 = varData.find(',')
            position=varData[index:index2]
            varData = varData[index2:]
            index = varData.find(':')
            index2 = varData.find(',')
            original_allele = varData[index:index2]
            varData = varData[index2:]
            index = varData.find(':')
            index2 = varData.find(',')
            mutant_allele = varData[index:index2]
            # print('altalllele== ', altalllele)
            chrcooridnates38 = str(chr) + ',' + chrpos + ',' + orientation + ',' + orgalllele + '/' + altalllele
            # print('chrcooridnates38== ', chrcooridnates38)
            chr_coord_list.append(chrcooridnates38)
            # return chrcooridnates38
        self.chr_coord_dict[var_id] = chr_coord_list

