import requests, sys, json
import re
import ast
class Grch37VariantData:

    def __init__(self):

        self.chr_coord_list =[]
        self.hgvs37Ids = []


    def parsevardatabystring(self,varData37, var_id):
        #getting orientation
        chr_coord_list_local = []
        hgvs37_local=[]
        #print('varData37==',varData37)
        #print("var_id", var_id)
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
                chr_coord_list_local.clear()
                hgvs37_local.clear()
                rsId=''
                hgvs37=''
                chrcooridnates37=''
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
                hgvs37='chr'+str(chr) + ':g.' + chrpos + orgalllele + '>' + altalllele
                #print('chrcooridnates37== ', chrcooridnates37)
                #print('hgvs37== ', hgvs37)
                rsId='rs'+var_id

                chr_coord_list_local.append(rsId)
                chr_coord_list_local.append(hgvs37)
                chr_coord_list_local.append(chrcooridnates37)
                self.chr_coord_list.append(chr_coord_list_local.copy())

                hgvs37_local.append(rsId)
                hgvs37_local.append(hgvs37)
                self.hgvs37Ids.append(hgvs37_local.copy())


