import os

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
UPLOAD_FOLDER = os.path.join(APP_ROOT, 'uploads')
RESULTS_FOLDER = os.path.join(APP_ROOT, 'results')
ALLOWED_EXTENSIONS = ['txt', 'csv']
DBWEB = ['SwissProt-Compact','SwissProt-Full']
SPFULL = 'uniprot_sprot.xml.gz'
DB = {'SwissProt-Compact':'compact_'+SPFULL, 'SwissProt-Full':SPFULL}
DATABASES = {
    'default': {
        'NAME': 'upep',
        'USER': 'root',
        'PASSWORD': 'upep2016',
        'HOST': 'localhost',
        'PORT': '',
        'DB' : 'uniprotdb',
    }
}
