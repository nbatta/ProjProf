
### code from orphics https://github.com/msyriac/orphics written by Mat Madhavacheril
def dict_from_section(config,section_name):
    try:
        del config._sections[section_name]['__name__']
    except:
        pass
    return dict([a, list_from_config(config,section_name,a)[0]] for a, x in list(config._sections[section_name].items()))

def list_from_string(string):
    return [float(x) for x in string.split(',')]

def list_from_config(Config,section,name):
    return list_from_string(Config.get(section,name))
