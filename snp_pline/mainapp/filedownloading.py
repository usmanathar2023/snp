# Import mimetypes module
import mimetypes
# import os module
import os

import pandas as pd
import xlwt
from io import BytesIO
# Import HttpResponse module
from django.http.response import HttpResponse
from django.http import FileResponse
from django.shortcuts import render, redirect
class FileDownload:
    def download_file(self, request, fileToDownload):
        # Define Django project base directory
        BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        # Define text file name
        filename = fileToDownload

        # Define the full file path
        filepath = BASE_DIR + filename
        #print('file to donwload===========' + filepath)
        # Open the file for reading content
        path = open(filepath, 'r')
        mime_type, _ = mimetypes.guess_type(filepath)
        response = HttpResponse(path, content_type=mime_type)
        response['Content-Disposition'] = 'attachment; filename='+filename + "'"
        #print("reponse===================== "+response.__str__())
        # Return the response value
        return response
    def downloadExcelFile(self, request, fileToDownload):
        # Define Django project base directory
        BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        # Define text file name
        filename = fileToDownload
        # Define the full file path
        filepath = BASE_DIR + filename
        f = open(filepath, 'rb')
        response = HttpResponse(f, content_type='application/vnd.ms-excel')
        response['Content-Disposition'] = 'attachment; filename=%s' % os.path.split(filepath)[-1]
        return response
        '''
        path = open(filepath, 'rb')
        response = HttpResponse(path, content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
        response['Content-Disposition'] = 'attachment; filename=' + filename + "'"
        # print("reponse===================== "+response.__str__())
        # Return the response value
        return response
        '''