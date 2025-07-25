import sys
import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

pdbs = sys.argv[1:]
pdbs.sort(key=natural_keys)


pdb = open('multiframe.pdb', 'w')

i = 1
for frame in pdbs:
    pdb.write('MODEL %s\n' % i)
    frame_ = open(frame, 'r').readlines()
    for line in frame_:
        if line.find('END') == -1:
            pdb.write(line)
    pdb.write('ENDMDL\n')
    i += 1

pdb.close()
