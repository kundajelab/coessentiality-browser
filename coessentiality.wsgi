#!/srv/www_coessentiality/coessentiality-browser/coessentiality-browser/bin/python                                                                     
activate_this = '/srv/www_coessentiality/coessentiality-browser/coessentiality-browser/bin/activate_this.py'
with open(activate_this) as file_:
    exec(file_.read(), dict(__file__=activate_this))
import sys
print(sys.path)
sys.path.insert(0,"/srv/www_coessentiality/coessentiality-browser")
import os
cwd = os.getcwd()
print(cwd)
from app import server as application