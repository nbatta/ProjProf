from configparser import SafeConfigParser
from src.io import dict_from_section

iniFile = "input/params.ini"
Config = SafeConfigParser()
Config.read(iniFile)

constDict = dict_from_section(Config,'constants')
cosmoDict = dict_from_section(Config,'params')

print (constDict)
