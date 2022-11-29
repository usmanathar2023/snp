import requests, sys, json
class Grch38VariantData:
    def getvariantdata(self, varid):
        server = "https://api.ncbi.nlm.nih.gov/variation/v0/"
        ext = "beta/refsnp/" + varid
        r = requests.get(server + ext)  # , headers={"Content-Type": "application/json"})
        spdiIdInfo = r.json()
        spdiIdInfo_str = json.dumps(spdiIdInfo)
        grch38Index = spdiIdInfo_str.find('GRCh38')
        grch37Index = spdiIdInfo_str.find('GRCh37')
        varData38 = spdiIdInfo_str[grch38Index:grch37Index]

        print('varData38== ', varData38)