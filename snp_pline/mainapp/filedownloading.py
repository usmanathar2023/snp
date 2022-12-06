# Import mimetypes module
import mimetypes
# import os module
import os
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
        print('file to donwload===========' + filepath)
        # Open the file for reading content
        path = open(filepath, 'r')
        mime_type, _ = mimetypes.guess_type(filepath)
        response = HttpResponse(path, content_type=mime_type)
        response['Content-Disposition'] = 'attachment; filename='+filename + "'"
        print("reponse===================== "+response.__str__())
        # Return the response value
        return response
