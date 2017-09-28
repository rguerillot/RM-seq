import pandoc
import os

# pandoc.core.PANDOC_PATH = '/usr/bin/pandoc'

doc = pandoc.Document()
doc.markdown = open('README.md').read()
with open('README.txt','w+') as f:
    f.write(doc.rst)
os.system("setup.py register")
os.remove('README.txt')
