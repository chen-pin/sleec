#! /usr/bin/env python

import shutil
import subprocess

try:
    shutil.rmtree('source')
except:
    pass
try: 
    shutil.rmtree('_build')
except:
    pass
subprocess.run(['sphinx-apidoc', '-f', '-o', 'source/', '../source/'])
subprocess.run(['make', 'html']) 
