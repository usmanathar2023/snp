import requests, sys, json
import re
import ast
class Grch37VariantData:
    def __init__(self):

        self.chr_coord_dict ={ }


    def parsevardatabystring(self,varData37, var_id):
        #getting orientation
        chr_coord_list = []

        grch37index=varData37.find('GRCh37')
        varData37 = varData37[grch37index:]
        orientation_index = varData37.find('orientation":')
        triplebracketsindex=varData37.find('}]}')
        varData37 = varData37[orientation_index+13:triplebracketsindex + 1]

        varData37=varData37.strip()
        #print('varData37== ', varData37)
        orientation_qoma_index=varData37.find(',')

        orientation= varData37[:orientation_qoma_index]
        if orientation=='false':
            orientation='1'
        else:
            orientation='0'

        # chr coordinates
        varData37 = varData37[orientation_qoma_index+1:]
        varData37Split=str(varData37).split('hgvs": "')
        for a in varData37Split:
            if a.find('>') > -1:
                #print('a==\n', a)
                greater_symbol_index=a.find('>')
                ncdata = a[:greater_symbol_index + 2]
                #separate chr coordinates as SIFT format
                #print('ncdata== ', ncdata)
                dotindex=ncdata.find('.')
                nc_no=ncdata[:dotindex]
                chr=nc_no[len(nc_no)-2:dotindex]
                chr = int(chr)
                #print('chr== ', chr)
                lastdotindex = ncdata.rfind('.')
                greater_symbol_index = ncdata.rfind('>')
                chrpos=ncdata[lastdotindex+1:greater_symbol_index-1]
                #print('chrpos== ', chrpos)
                orgalllele = ncdata[greater_symbol_index - 1:greater_symbol_index]
                #print('orgalllele== ', orgalllele)
                altalllele = ncdata[greater_symbol_index + 1:len(ncdata)]
                #print('altalllele== ', altalllele)
                chrcooridnates37=str(chr)+','+chrpos+','+orientation+','+orgalllele+'/'+altalllele
                #print('chrcooridnates37== ', chrcooridnates37)
                chr_coord_list.append(chrcooridnates37)
                #return chrcooridnates37
                #print('chr_coord_list37== ', chr_coord_list)
            self.chr_coord_dict['rs'+var_id]=chr_coord_list