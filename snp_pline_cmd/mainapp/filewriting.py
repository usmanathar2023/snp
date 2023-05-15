

class FileWriting:
    def writeFileFromList(self, list, fname):
        with open(fname, 'w') as fp:
            fp.write('\n'.join(list))

    def writeProteinData(self, protData,prot_var_data_dict,prot_var_data_filename):
        with open(prot_var_data_filename, 'w') as fp:
            fp.write(protData)
            protein_ids=prot_var_data_dict['refSeqProtId']
            fp.write('\nRefSeq_Protein Id: \n'+protein_ids[0])
            fp.write('\nProtein Allelic Data:\n')
            fp.write('\n'.join(prot_var_data_dict['protCoord']))

