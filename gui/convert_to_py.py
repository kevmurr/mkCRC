from PyQt5 import uic
fin = open('parameters.ui','r')
fout = open('parameters.py','w')
uic.compileUi(fin,fout,execute=False)
fin.close()
fout.close()