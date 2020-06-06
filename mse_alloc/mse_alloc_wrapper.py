#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import cgi
import cgitb
from mse_alloc_dev import MseFibAlloc

cgitb.enable()

print("Content-Type: text/html")
print("")

keywords = {}

# Get the querystring
arguments = cgi.FieldStorage()

# Get the file
fileitem = arguments['file_id']
# Session ID
keywords['sessionID'] = arguments.getvalue('sessionID')
# Instrument settings
keywords['spectro'] = arguments.getvalue('spectro')
# Method parameters
keywords['meth'] = arguments.getvalue('meth')
keywords['dither'] = arguments.getvalue('dither')
keywords['iternum'] = float(arguments.getvalue('fixiternum'))
keywords['allocfrac'] = float(arguments.getvalue('allocfrac'))
# Targets settings
keywords['fovctr'] = arguments.getvalue('fov_ctr')
keywords['fovctrra'] = float(arguments.getvalue('fov_ctr_ra'))
keywords['fovctrdec'] = float(arguments.getvalue('fov_ctr_dec'))

# Check the first line of the file
chunk = fileitem.file.read(37)
if chunk != '|  RAJ2000|  DEJ2000|priority|repeat|':
    print("Your file does not have the right format.")
# If first line ok, then read the entire file and write it to the server
else:
    fileitem.file.seek(0)
    with open('tmp/'+keywords['sessionID']+'targets.dat', 'w') as f:
        f.write(fileitem.file.read())

    # Print list of keywords on the page, so we can read them when querying with the URL
#    print(keywords)

    alloc = MseFibAlloc(**keywords)
    print(" <h2> Allocation fraction </h2>")
    for i in range(len(alloc.trackfrac)):
        print(" Iteration #{0:d}: {1:.1%}".format(i+1, alloc.trackfrac[i]))
        print("<br>")

    # Delete "sessionID"_target.dat file
    os.remove('tmp/'+keywords['sessionID']+'targets.dat')