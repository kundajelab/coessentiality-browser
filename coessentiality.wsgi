#!/usr/bin/python
import sys
sys.path.insert(0,"/var/www/coessentiality-browser")
print("path:"+str(sys.path))

import os
cwd = os.getcwd()
print("cwd:"+str(cwd))
from app import server as application
