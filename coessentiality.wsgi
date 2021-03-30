#!/usr/bin/python
import sys
sys.path.insert(0,"/var/www/coessentiality-browser")
print("path:"+sys.path)

import os
cwd = os.getcwd()
print("cwd:"+cwd)
from app import server as application
