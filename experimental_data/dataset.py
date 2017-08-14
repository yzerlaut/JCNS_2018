import sys, os, glob
sys.path.append('../../')

def get_duration(filename):
    s = filename.split('_')[-1].split('.mat')[0]
    s1, s2 = '', ''
    s1, i = s[0], 1
    while s[i]=='0':
        s1 += s[i]
        i+=1
    s2, i = s[i], i+1
    while ((i<len(s)) and (s[i]=='0')):
        s2 += s[i]
        i+=1
    return float(s1), float(s2)
    
def get_dataset():
    DATA = []
    os.chdir("../experimental_data/")
    for f in glob.glob("Monkey1"+os.path.sep+"*.mat"):
        DATA.append({'Monkey':'1', 'filename':f,
                     'duration':get_duration(f)[0]})
    for f in glob.glob("Monkey2"+os.path.sep+"*.mat"):
        DATA.append({'Monkey':'2', 'filename':f,
                     'duration':get_duration(f)[0]})
    return DATA

if __name__=='__main__':
    import pprint
    pprint.pprint(get_dataset())


