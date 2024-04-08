import wts
import glob
import pandas as pd
#a = tables.porphyr("")

pdbc = "1cyj"


#pdbc = "1ewh"

if(True):
  #  wts.TCL_Skript()
    #wts.dihedpdb(read_keep=True)        
    wts.dihedpdb()        
    wts.prepare_gaussian_logs().core()
    wts.onecsv()    
#   
print("finish   ")    