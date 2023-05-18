import views


class StartingIPSNP:
    def processingOptions(self):
        userOptions=input('Enter 1 for annotation retrieval 2 for chromosomal or proteomic data 3 for exiting\n').strip()
        print('userOptions  ',userOptions)
        match userOptions:
            case '1':
                print('1')
                views.annotateVariants()
            case '2':
                print('2')
                views.vardataretrievalprocessing()
            case '3':
                print('3')
                quit()


            case other:
                print('No match found')
        
startingIPSNP=StartingIPSNP()
startingIPSNP.processingOptions()
