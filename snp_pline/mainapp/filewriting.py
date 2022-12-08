

class FileWriting:
    def writeFileFromList(self, list, fname):
        with open(fname, 'w') as fp:
            fp.write('\n'.join(list))
