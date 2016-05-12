import os

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
UPLOAD_FOLDER = os.path.join(APP_ROOT, 'uploads')
RESULTS_FOLDER = os.path.join(APP_ROOT, 'results')
ALLOWED_EXTENSIONS = ['txt', 'csv']
