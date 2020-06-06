#!/usr/bin/python
# -*- coding: utf-8 -*-

import cgi
import cgitb

cgitb.enable()

print("Content-Type: text/html")
print("")

keywords = {}

# Get the querystring (should only be the file)
arguments = cgi.FieldStorage()

# Get the file and sessionID
fileitem = arguments['file_id']
sessionID = arguments['sessionID']
# Check the first line
chunk = fileitem.file.read(37)
if chunk != '|  RAJ2000|  DEJ2000|priority|repeat|':
    print("Your file does not have the right format.")
# If first line ok, then read the entire file and write it to the server
else:
    with open('tmp/'+sessionID+'targets.dat', 'w') as f:
        f.write(fileitem.file.read())
