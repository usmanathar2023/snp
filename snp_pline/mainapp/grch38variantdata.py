import requests, sys, json
import re
import ast
class Grch38VariantData:
    def __init__(self):
        self.chr_coord_dict ={ }

    def parsevardatabystring(self,varData38,  var_id):
        #getting orientation
        chr_coord_list = []

        orientation_index=varData38.find('orientation":')
        triplebracketsindex=varData38.find('}]}')

        varData38=varData38[orientation_index+13:triplebracketsindex+1]
        varData38=varData38.strip()

        orientation_qoma_index=varData38.find(',')

        orientation= varData38[:orientation_qoma_index]
        if orientation=='false':
            orientation='1'
        else:
            orientation='0'

        # chr coordinates
        varData38 = varData38[orientation_qoma_index+1:]
        varData38Split=str(varData38).split('hgvs": "')
        for a in varData38Split:
            if a.find('>') > -1:
                print('a==\n', a)

                '''
                varData38=str(varData38).lower()
                print('varData38===after', varData38)
        
                txt = "hello planet hello"
        
                # Search for a sequence that starts with "he", followed by 1 or more  (any) characters, and an "o":
        
                x = re.findall("el.+o", txt)
                print('x  ',x)
               
        
                nc_indices=re.findall(r"\nc_.+}", varData38) #nc_.+\}$
                for a in nc_indices:
                    print('a==\n',a)
                 '''
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
                chrcooridnates38=str(chr)+','+chrpos+','+orientation+','+orgalllele+'/'+altalllele
                #print('chrcooridnates38== ', chrcooridnates38)
                chr_coord_list.append(chrcooridnates38)
                #return chrcooridnates38
                #print('chr_coord_list38 == ', chr_coord_list)
            self.chr_coord_dict['rs'+var_id] = chr_coord_list



is38=False
is37=False
dnavarinfo37=''; dnavarinfo38=''
orientation37=''; orientation38=''
proteinvarinfo37=''; proteinvarinfo38='';v_dic='';v_str =''
def iterdict(d):
    global is38, is37,v_dic,v_str , dnavarinfo37, dnavarinfo38, orientation37, orientation38, proteinvarinfo37, proteinvarinfo38;
    #print('d== ', d)
    siftparams=''
    for k, v in d.items():
        #print('type v== ', type(v))
        if isinstance(v, list):
            v = str(v).strip()
            v = v.removeprefix('[')
            v = v.removesuffix(']')
            if v.startswith('{'):
                v = json.dumps(v)
                v = json.loads(v)
        if isinstance(v, tuple):
            v = str(v).strip()
            v = v.removeprefix('(')
            v = v.removesuffix(')')
            if v.startswith('{'):
                v = json.dumps(v)
                v = json.loads(v)
        if isinstance(v, str) and v.startswith('{'):
            v_str = json.dumps(v)
            v_dic = json.loads(v_str)
            print('v type == ', type(v_dic))
            print('v== inside check string', v_dic)
        if isinstance(v, dict):
            iterdict(v)
        else:
            #print('v===',v)
            v=str(v).strip()

            if v.find('GRCh38') > -1:
                is38 = True
                is37 = False
                print('is38=== ', is38)
            elif v.find('GRCh37') > -1:
                is37 = True
                is38 = False
            if is38:
                if v.find('false')>-1:
                     orientation38 = 1
                if (v.find('NC_') >-1 and v.find('>')>-1):
                    dnavarinfo38 = v
                    print('dnavarinfo38== ', dnavarinfo38)
                if (v.startswith('NP_') and v.endswith('')):
                    proteinvarinfo38 = v
            elif is37:
                if v.find('false')>-1:
                    orientation37 = 1
                if (v.find('NC_') >-1 and v.find('>')>-1):
                    dnavarinfo37 = v
                    print('dnavarinfo37== ', dnavarinfo37)
                if (v.find('NP_')>-1 and v.endswith('')):
                    proteinvarinfo37 = v

                siftparams = str(dnavarinfo37) + str(orientation37)
                print('siftparams== ', siftparams)