import os

APP_ROOT = os.path.dirname(os.path.abspath(__file__))
APP_STATIC = os.path.join(APP_ROOT, 'static')
APP_TEMPLATE = os.path.join(APP_ROOT, 'templates')
MAX_WORKERS = 4

BROKER_URL = 'redis://localhost:6379/0'