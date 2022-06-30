# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 6.14-1 replay file
# Internal Version: 2014_06_05-00.11.02 134264
# Run by oyvinpet on Fri Aug 11 13:27:53 2017
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(1.34896, 1.35), width=198.567, 
    height=133.92)
session.viewports['Viewport: 1'].makeCurrent()
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
execfile('C:/Work/Hardanger/Forceid/smoothing/ssbeam_test/exportModal.py', 
    __main__.__dict__)
#: Model: C://Work//Hardanger//Forceid//smoothing//ssbeam_test//abaqus//ssbeam.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             2
#: Number of Element Sets:       5
#: Number of Node Sets:          7
#: Number of Steps:              1
#: 0.0414372065105
#* Exit code: 0
#* File "C:/Work/Hardanger/Forceid/smoothing/ssbeam_test/exportModal.py", line 
#* 174, in <module>
#*     sys.exit(0)
