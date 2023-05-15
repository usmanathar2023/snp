import requests
class RefSeqId_to_UniProtId:
    def get_uniprotid_from_refseqid(self, term,refSeqProtId):

        #url = 'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=(gene:' + term +') AND (organism_id:9606) AND (xref:refseq-' + refSeqProtId[0]+ ')' #NP_000489.3)'
        url = 'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=(gene:' + term + ') AND (organism_id:9606) AND '  +  refSeqProtId[0]

        print('url\n', url)
        results = requests.get(url).text
        index=results.find('|')
        results=results[index+1:]
        index = results.find('|')
        prot_id=results[:index]
        index2 = results.find(' ')
        prot_name = results[index+1:index2]
        newlineIndex = results.find('\n')
        protein=results[newlineIndex+1:]
        protData= 'Gene Name: \n' + term + '\nUniProt ID:\n'+ prot_id + '\nUniProt Name:\n ' + prot_name + '\nProtein Data: \n' + protein
        print('prot_id\n', prot_id)
        return protData


