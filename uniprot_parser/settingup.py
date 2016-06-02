import os

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
UPLOAD_FOLDER = os.path.join(APP_ROOT, 'uploads')
RESULTS_FOLDER = os.path.join(APP_ROOT, 'results')
DB_FOLDER = os.path.join(APP_ROOT, 'database')
ALLOWED_EXTENSIONS = ['txt', 'csv']
DBWEB = ['SwissProt-Compact','SwissProt-Full']
SPFULL = 'uniprot_sprot.xml.gz'

DATABASES = {
    'default': {
        'NAME': 'schulzlab',
        'USER': 'root',
        'PASSWORD': 'schulzlab',
        'HOST': 'localhost',
        'PORT': '',
        'DB' : 'uniprotdb',
    }
}
DBXML = dict()
TAXA = ['Vertebrata','Mammalia','Eukaryota','Fungi']
TABLE_TYPE = ['ACC','FullName','EC','Organism','SubLoc','GO','FT','Seq']
TABLE_NAME = TAXA
TABLE_NAME.append('Compact')
DBMYSQL = {'SwissProt':'sprot','Trembl':'trembl'}
DBTYPEMSQL = {'SwissProt':'sp','Trembl':'tr'}
for i in TABLE_NAME:
    DBXML[i] = i+'_Uniprot.xml.gz'

