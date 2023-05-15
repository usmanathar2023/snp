from mainapp import views


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
        '''
        urlpatterns=[
        path('', views.index, name='index'),

        path('vardataretrieval/vardataretrievalprocessing', views.vardataretrievalprocessing, name='vdataretproc'),#/<varnotfoundfile>
        path('runtools', views.runTools, name='runtools'),
        path('runtools/annotatevariants', views.annotateVariants, name='annotate_variants'),
        path('runtools/selectedvarsdataretrieval', views.selectedVarsDataRetrieval, name='selectedvarsdataretrieval'),
        path('vardataretrieval', views.varDataRetrieval, name='vardataretrieval'),
        path('vardataretrieval/download', views.downloadRSIdFile, name="downloadRSId"),
        path('vardataretrieval/download_grch37data', views.downloadGRCh37File, name="downloadGRCh37"),
        path('vardataretrieval/download_grch38data', views.downloadGRCh38File, name="downloadGRCh38"),
        path('vardataretrieval/download_Protdata', views.downloadProteinDataFile, name="downloadProData"),
        path('runtools/annotatevariants/varannotdownload/<tool>', views.downloadVarAnnotation, name="downloadvarannotation"),
        path('selvardataretrievalprocessing',views.selVarDataRetrievalProcessing, name='selvdataretproc')
        ]'''
startingIPSNP=StartingIPSNP()
startingIPSNP.processingOptions()