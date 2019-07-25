#from configparser import SafeConfigParser
#from src.io import dict_from_section
from src.Project import Profiles

#iniFile = "input/params.ini"
#Config = SafeConfigParser()
#Config.read(iniFile)

#constDict = dict_from_section(Config,'constants')
#cosmoDict = dict_from_section(Config,'params')

#print (constDict)

P = Profiles()

#print (P.cc)
#print (P.pp)
print (P.om)

print (P.AngDist(0.5))
