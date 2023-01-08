from django.urls import path
from . import views

urlpatterns=[
    path('', views.index, name='index'),

    path('vardataretrieval/vardataretrievalprocessing', views.vardataretrievalprocessing, name='vdataretproc'),
    path('runtools', views.runTools, name='runtools'),
    path('runtools/annotatevariants', views.annotateVariants, name='annotate_variants'),
    path('vardataretrieval', views.varDataRetrieval, name='vardataretrieval'),
    path('vardataretrieval/download', views.downloadRSIdFile, name="downloadRSId"),
    path('vardataretrieval/download_grch37data', views.downloadGRCh37File, name="downloadGRCh37"),
    path('vardataretrieval/download_grch38data', views.downloadGRCh38File, name="downloadGRCh38"),
    path('vardataretrieval/download_Protdata', views.downloadProteinDataFile, name="downloadProData"),

]