import os
import pickle
import logging
import time

import tkinter as tk
from tkinter import ttk
from threading import Thread

from google.auth.transport.requests import Request
from google_auth_oauthlib.flow import InstalledAppFlow
from google.oauth2.credentials import Credentials
from googleapiclient.discovery import build
from googleapiclient.http import MediaFileUpload

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

SCOPES = ['https://www.googleapis.com/auth/drive.file']

SharedDrive = '1L83XHhzW9ymdcMrrQYiiU0wIk9dFwbv6'

def authenticate_drive():
    """Authenticate the user and return the Drive service."""
    creds = None
    token_path = 'token.pickle'
    credentials_path = 'credentials.json'

    # The file token.pickle stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first time.
    if os.path.exists('token.pickle'):
        with open('token.pickle', 'rb') as token:
            creds = pickle.load(token)

    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
            logger.info("Credentials refreshed.")
        else:
            flow = InstalledAppFlow.from_client_secrets_file(
                'credentials.json', SCOPES)
            creds = flow.run_local_server(port=0)
            logger.info("New credentials obtained.")

        # Save the credentials for the next run
        with open('token.pickle', 'wb') as token:
            pickle.dump(creds, token)
            logger.info(f"Credentials saved to {token_path}.")

    service = build('drive', 'v3', credentials=creds)
    return service


def create_folder(service, folder_name, directory = SharedDrive):
    query = (
        f"name = '{folder_name}' and "
        f"mimeType = 'application/vnd.google-apps.folder' and "
        f"'{directory}' in parents and "
        "trashed = false"
    )
    try:
        # Search for the folder
        response = service.files().list(
            q=query,
            spaces='drive',
            fields='files(id, name)',
            pageSize=10
        ).execute()

        files = response.get('files', [])

        if files:
            # Folder exists; return the first matching folder's ID
            folder_id = files[0].get('id')
            print(f"Folder '{folder_name}' already exists with ID: {folder_id}")
            return folder_id
        else:
            # Folder does not exist; create it
            file_metadata = {
                'name': folder_name,
                'mimeType': 'application/vnd.google-apps.folder',
                'parents': [directory]
            }
            folder = service.files().create(
                body=file_metadata,
                fields='id'
            ).execute()
            folder_id = folder.get('id')
            print(f"Created folder '{folder_name}' with ID: {folder_id}")
            return folder_id

    except Exception as e:
        print(f"An error occurred: {e}")
        return None


def upload_file(service, directory, file_path, mime_type='text/csv', chunk_size=5 * 1024 * 1024, max_retries=5):
    """
    Uploads a file to Google Drive using resumable uploads to handle large files.

    :param service: Authorized Google Drive API service instance.
    :param directory: The ID of the parent directory where to upload the file.
    :param file_path: The path to the file to upload.
    :param mime_type: The MIME type of the file. Defaults to 'text/csv'.
    :param chunk_size: The size of each chunk in bytes. Default is 5MB.
    :param max_retries: Number of retries in case of transient failures.
    :return: The ID of the uploaded file, or None if failed.
    """
    file_name = os.path.basename(file_path)
    mime_type = mime_type or 'application/octet-stream'  # Default MIME type

    file_metadata = {
        'name': file_name,
        'parents': [directory]
    }

    media = MediaFileUpload(
        file_path,
        mimetype=mime_type,
        resumable=True,
        chunksize=chunk_size
    )

    request = service.files().create(
        body=file_metadata,
        media_body=media,
        fields='id',
        supportsAllDrives=True
    )

    response = None
    retry = 0

    while response is None:
        try:
            status, response = request.next_chunk()
            if status:
                progress = int(status.progress() * 100)
                logger.info(f"Upload {progress}% complete for {file_name}")
        except HttpError as error:
            if error.resp.status in [403, 500, 502, 503, 504]:
                if retry < max_retries:
                    retry += 1
                    sleep_time = 2 ** retry
                    logger.warning(f"Server error {error.resp.status}. Retrying in {sleep_time} seconds...")
                    time.sleep(sleep_time)
                    continue
                else:
                    logger.error(f"Failed to upload {file_name} after {max_retries} retries.")
                    return None
            else:
                logger.error(f"An error occurred: {error}")
                return None
        except Exception as e:
            logger.error(f"An unexpected error occurred: {e}")
            return None

    logger.info(f"Uploaded {file_name} with File ID: {response.get('id')}")
    return response.get('id')


def run_upload_with_separate_thread(service, directory, file_path, chunk_size=5 * 1024 * 1024):
    def upload_task():
        file_id = upload_file(
            service=service,
            directory=directory,
            file_path=file_path,
            chunk_size=chunk_size
        )
        if file_id:
            logger.info(f"Upload completed successfully. File ID: {file_id}")
        else:
            logger.error("Upload failed.")

    upload_thread = Thread(target=upload_task, daemon=True)
    upload_thread.start()
