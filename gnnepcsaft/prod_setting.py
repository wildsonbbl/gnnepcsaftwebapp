import os

import django_heroku

# default: use settings from main settings.py if not overwritten
from .settings import *

DEBUG = False
SECRET_KEY = os.getenv("DJANGO_SECRET_KEY", SECRET_KEY)
# adjust this to the URL of your Heroku app
ALLOWED_HOSTS = ["gnnepcsaft.herokuapp.com"]
# Activate Django-Heroku.
django_heroku.settings(locals())
