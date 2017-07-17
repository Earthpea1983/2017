import os

plist = os.popen('pip list')
plist = plist.read()
plist = plist.split('\n')
[plist.pop(0) for i in range(2)]
plist.pop(-1)
for i in range(len(plist)):
    plist[i] = plist[i].split(' ')
for i in range(len(plist)):
    plist[i] = plist[i][0]
for i in range(len(plist)):
    os.system('sudo pip install --upgrade {}'.format(plist[i]))

