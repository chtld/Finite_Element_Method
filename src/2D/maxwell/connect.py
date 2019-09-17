import os
filelist=os.listdir('.')
for item in filelist:
    print item

newfile=open('/root/Music/new.txt','w')
for item in filelist:
    for txt in open(item,'r'):
        newfile.write(txt)

newfile.close()