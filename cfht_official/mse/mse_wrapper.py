#!/usr/bin/python
# Wrapper to parse arguments from the cgi into something understandable by python

import cgi, cgitb
from mse_itc import MseSpectrum

cgitb.enable()

print("Content-Type: text/html")
print("")

keywords = {}

# Get the querystring
arguments = cgi.FieldStorage()

# Session ID
keywords['sessionID'] = arguments.getvalue('sessionID')
# Observing conditions
keywords['seeing']  = float(arguments.getvalue('seeing'))
keywords['airmass'] = float(arguments.getvalue('airmass'))
keywords['skymag']  = float(arguments.getvalue('skymag'))
# Method parameters
keywords['meth']  = arguments.getvalue('meth')
keywords['etime'] = float(arguments.getvalue('etime'))
keywords['snr']   = float(arguments.getvalue('snr'))
# Source settings
keywords['band']     = arguments.getvalue('band')
keywords['template'] = arguments.getvalue('template')
keywords['tgtmag']   = float(arguments.getvalue('tgtmag'))
keywords['redshift'] = float(arguments.getvalue('redshift'))
keywords['src_type'] = arguments.getvalue('src_type')
# Telescope parameters
keywords['coating'] = arguments.getvalue('coating')
# Instrument settings
keywords['spectro'] = arguments.getvalue('spectro')
if keywords['spectro'] == 'HR':
    if arguments.getvalue('fibdiam') == 0:
        keywords['fibdiam'] = 0.7
    else:
        keywords['fibdiam'] = 0.8
else:
    if arguments.getvalue('fibdiam') == 0:
        keywords['fibdiam'] = 0.9
    else:
        keywords['fibdiam'] = 1.0
# Binning settings
keywords['specbin'] = float(arguments.getvalue('specbin'))
keywords['spatbin'] = float(arguments.getvalue('spatbin'))
    

# Print list of keywords on the page, so we can read them when querying with the URL
print(keywords)

spec = MseSpectrum(**keywords)
res, script, div = spec.compute_snr(doplot='online')
print(div)
print(script)