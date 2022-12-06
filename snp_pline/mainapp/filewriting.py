

class FileWriting:
    def writeFileFromList(self, list):
        with open(r'media/rsids.txt', 'w') as fp:
            fp.write('\n'.join(list))
