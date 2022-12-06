import requests, sys, json
class DBSnpVarientDataRetrieval:
    def getvariantdata(self, varid):
        server = "https://api.ncbi.nlm.nih.gov/variation/v0/"
        ext = "beta/refsnp/" + varid
        r = requests.get(server + ext)  # , headers={"Content-Type": "application/json"})
        spdiIdInfo = r.json()
        spdiIdInfo_str = json.dumps(spdiIdInfo)
        # spdiIdInfo_dict = json.loads(spdiIdInfo_str)
        # iterdict(spdiIdInfo_dict)
        # grch38Index = spdiIdInfo_str.find('GRCh38')
        # grch37Index = spdiIdInfo_str.find('GRCh37')
        # varData38 = spdiIdInfo_str[grch38Index:grch37Index]
        # self.parsevardatabystring(varData38)
        return spdiIdInfo_str