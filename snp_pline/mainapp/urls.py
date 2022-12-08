from django.urls import path
from . import views

urlpatterns=[
    path('', views.index, name='index'),
    path('routerequest', views.routeRequest, name='routereq'),
    path('download', views.downloadFile, name="download"),
    path('download_grch37data', views.downloadGRCh37File, name="downloadGRCh37"),
    path('download_grch38data', views.downloadGRCh38File, name="downloadGRCh38"),




]