##!/usr/bin/env python
import math
import numpy as np
import os, sys
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import re
from htmd.ui import *

def pdb_filter():
    need = open("fis1.pdb",'w')
    with open("fis.pdb",'r') as f:
        for line in f:
            if not re.match(r'^\s*$|^END|^TER|^MODEL|^CONECT',line):
                need.write(line)
    need.close()

def correctpdb():
    atom = []; atno = []; atyp = []; altyp = []; res = []; chain = [];
    resno = []; altid = []; x = []; y = []; z = []; occ = []; Bfac = [];
    ele = []; charge = []
    with open('fis1.pdb', 'r') as f:
        for line in f:
            list1 = [atom, atno, atyp, altyp, res, chain, resno, altid, x, y, z,
                     occ, Bfac, ele, charge]
            list2 = [line[0:6], line[6:11], line[12:16], line[16:17], line[17:20],
                     line[21:22], line[22:26], line[26:27], line[30:38], line[38:46],
                     line[46:54], line[54:60], line[60:66], line[76:78], line[78:80]]
            for i in range(len(list1)):
                list1[i].append(list2[i])
    df1 = pd.DataFrame({'atom': atom, 'atno': atno, 'atyp': atyp, 'altyp': altyp,
                        'res': res, 'chain': chain, 'resno': resno, 'altid': altid,
                        'X': x, 'Y': y, 'Z': z, 'occ': occ, 'B': Bfac, 'ele': ele,
                        'charge': charge})
    return df1
#df1 = pd.read_csv('fis1.pdb', sep='\s+',header=None)
#df1.columns = ['atom', 'atno', 'atyp', 'res', 'chain', 'resno', 'X', 'Y', 'Z', 'A', 'B', 'seg', 'ele']
def formatdf(df1):
    mist = ['atom', 'atno', 'atyp', 'altyp', 'res', 'chain', 'resno',
        'altid', 'X', 'Y', 'Z', 'occ', 'B', 'ele','charge']
    floats = ['X', 'Y', 'Z', 'occ', 'B']
    ints = ['atno','resno']
    for i in mist:
        df1[i]=df1[i].str.strip()
    for i in floats:
        df1[i]=df1[i].astype(float)
    for i in ints:
        df1[i]=df1[i].astype(int)
    return df1

def filter_psf_return_bonds():
    a = 0; b = 0
    temp = open("test.now",'w')
    bondpairs=[];
    with open("fis.psf", 'r') as f:
        for line in f:
            if not re.match(r'^\s*$', line):
                if re.search("!NATOM", line):
                    a+=1
                if re.search("!NBOND: bonds", line):
                    a=0
                    b+=1
                if re.search("!NTHETA: angles", line):
                    b=0
                if a>0:
                    if 'NATOM' not in line:
                        temp.write(line)
                if b>0:
                    if 'NBOND: bonds' not in line:
                        simples=[]
                        dimples=re.split("\s+",line)
                        for i in dimples:
                            if not re.match(r'^\s*$', i):
                                simples.append(i)
                        sam=range(0, len(simples), 2)
                        for t in (sam):
                            meekin=simples[t:t+2]
                            bondpairs.append(meekin)
    temp.close()
    return bondpairs

def read_psf_get_bonds():
    df2 = pd.read_csv('test.now', sep='\s+',header=None)
    df2.columns = ['atno', 'seg', 'resno', 'res', 'atyp', 'atyp1', 'deltaq', 'mass', 'restrain']
    team=[]
    for r in df2['atno']:#could be df1 too!
        point1 = df.loc[df[0]==str(r)]
        point2 = df.loc[df[1]==str(r)]
        first = point1[1].values
        second = point2[0].values
        bonds = np.concatenate((first, second), axis=0)
        check=bonds.tolist()
        print(check)
        team.append(bonds)
    bondlength=[]
    for x in team:
        bondlength.append((len(x)))
    return team, bondlength, df2

def psf_pdb_match(df1):
    dfmain = pd.concat([df1['atno'],df2['atno'],df1['atyp'],df2['atyp']], axis=1, keys=['atno','atno1','atyp','atyp1'])
    q=0
    for r in range(len(dfmain)):
        if dfmain['atno'].values[r]==dfmain['atno1'].values[r]:
            if dfmain['atyp'].values[r]==dfmain['atyp1'].values[r]:
                q+=1
    if q==len(dfmain):
        print('psf and pdb index match: OK!')
    else:
        print('psf and pdb index match: FAIL!')
    return dfmain

def run_gaussian_for_field(dflig, dfbin, nproc, method, lig, car):
    QMS0 = open("QMS0.com",'w')
    QMS0.write('%%Chk=QMS.chk\n%%Nproc=%s\n#P %s nosymm\n\n%s izol\n\n' % (nproc, method, lig))
    QMS0.close()
    QMS = open("QMS",'w')
    QMS.write('%0.f 1\n' % car)
    for index, row in dflig.iterrows():
        QMS.write(('%-2s %14.6f%14.6f%14.6f\n') %
                  (row['ele'],row['X'],row['Y'],row['Z']))
    QMS.write('\n')
    QMS.close()
    subprocess.Popen("cat QMS >>QMS0.com", shell=True).communicate()
    subprocess.Popen("G16 QMS0.com", shell=True).communicate()
    subprocess.Popen("cp QMS.chk QMS0.chk", shell=True).communicate()
    QMScom = open("QMS.com",'w')
    QMScom.write('%%Chk=QMS.chk\n%%Nproc=%s\n#P %s nosymm guess=read prop=read charge\n\n%s perm_chg\n\n'
                 % (nproc, method, lig))
    QMScom.close()
    subprocess.Popen("cat QMS >>QMS.com", shell=True).communicate()
    QMScom = open("QMS.com",'a')
    for index, row in dfbin.iterrows():
        QMScom.write(('%.6f %.6f %.6f %.5f\n') %
                     (row['X'],row['Y'],row['Z'],row['chg']))
    QMScom.write('\n')
    for index, row in dfbin.iterrows():
        QMScom.write(('%.6f %.6f %.6f\n') %
                     (row['X'],row['Y'],row['Z']))
    QMScom.write('\n')
    QMScom.close()
    print ('%s atoms saved of %s atoms' % (len(dfbin1),len(dfbin)))
    subprocess.Popen("G16 QMS.com", shell=True).communicate()


def readESP_G(field, nudge):
    c = 0; d = 0
    mead = open(field, 'w')
    with open(nudge, 'r') as f:
        for line in f:
            if re.search("Potential", line):
                c += 1
            if re.search("Gradient", line):
                c = 0
            if c > 0:
                if not re.search('----|Potential|Gradient|Atom', line):
                    d += 1
                    mead.write(line)
    mead.close()
    return d

def polasignMMmae_def(dfbin):
    alpha = []
    polka = []
    emap = set()
    for r in range(len(dfbin)):
        emap.add(dfbin['ele'][r])
        polka.append(0)
        if dfbin['ele'][r] == 'H':
            alpha.append(0.387)
        if dfbin['ele'][r] == 'S':
            alpha.append(3.000)
        if dfbin['ele'][r] == 'O':
            if dfbin['bondno'].values[r] == 1:
                alpha.append(0.569)
            else:
                alpha.append(0.637)
        if dfbin['ele'][r] == 'C':
            if dfbin['bondno'].values[r] == 4:
                alpha.append(1.061)
            if dfbin['bondno'].values[r] == 3:
                if dfbin['res'][r] == 'ARG':
                    if dfbin['atyp'][r] == 'CZ':
                        alpha.append(1.896)
                    else:
                        alpha.append(1.352)
                else:
                    alpha.append(1.352)
        if dfbin['ele'][r] == 'N':
            if dfbin['res'][r] == 'LYS':
                if dfbin['atyp'][r] == 'NZ':
                    alpha.append(0.964)
                else:
                    alpha.append(1.090)
            else:
                if dfbin['res'][r].startswith('H'):
                    if dfbin['res'][r] == 'HSP':
                        alpha.append(1.090)
                    if dfbin['res'][r] == 'HSE':
                        if dfbin['atyp'][r] == 'ND1':
                            alpha.append(1.030)
                        else:
                            alpha.append(1.090)
                    if dfbin['res'][r] == 'HSD':
                        if dfbin['atyp'][r] == 'NE2':
                            alpha.append(1.030)
                        else:
                            alpha.append(1.090)
                else:
                    alpha.append(1.090)
    return alpha, polka

def enex():
    with open("ene2", 'r') as f:
        for line in f:
            enex = re.match('^\s*([-+]?\d+\.\d+)$', line)
            enep = float(enex.group(0))
    return enep

def cycle(nproc, method, lig, num):
    subprocess.Popen("cp atoms2 atomsP%s" % num, shell=True).communicate()
    subprocess.Popen("cp com comP%s" % num, shell=True).communicate()
    subprocess.Popen("cp centr centrP%s" % num, shell=True).communicate()
    QMSP0com = open("QMSP%s.com" % num, 'w')
    QMSP0com.write('%%Chk=QMSP.chk\n%%Nproc=%s\n#P %s nosymm guess=read prop=read charge\n\n%s polar_chg %s\n\n'
                   % (nproc, method, lig, num))
    QMSP0com.close()
    subprocess.Popen("cat QMS >>QMSP%s.com" % num, shell=True).communicate()
    subprocess.Popen("cat comP%s >>QMSP%s.com" % (num, num), shell=True).communicate()
    subprocess.Popen("cp QMS0.chk QMSP.chk", shell=True).communicate()
    subprocess.Popen("G16 QMSP%s.com" % num, shell=True).communicate()

def converge(enep, nproc, method, lig):
    bstd = 1;
    num = 1
    while (bstd > convergence):
        numa = str(num - 1)
        nudge = ("QMSP%s.com.log" % numa)
        field = ("field%s" % num)
        FCS = ("FCS%s" % num)
        fiel = readESP_G(field, nudge)
        # output = subprocess.Popen('wc -l centrP%s' % numa, shell =True, stdout=subprocess.PIPE).communicate()[0]
        # ccc = int(output.split(' ')[0])
        # if fiel == ccc:
        #    print ("YES4"):
        subprocess.Popen("paste -d ' ' %s centrP%s >%s" % (field, numa, FCS), shell=True).communicate()
        subprocess.Popen("./calc_indpolF %s" % FCS, shell=True).communicate()
        subprocess.Popen("cp ene2 ene2P%s" % num, shell=True).communicate()
        ene = enex()
        cycle(nproc, method, lig, num)
        bstd = abs(ene - enep) / ene
        print("Energy %s\n" % ene)
        print("Difference %s\n" % bstd)
        enep = ene
        num += 1

def get_charges():
    bit = ['QMS', 'QMSP']
    child_processes = []
    for i in range(len(bit)):
        subprocess.Popen("formchk %s.chk" % bit[i], shell=True).communicate()
        GDATA = open("DATA%s.GDMA" % i, 'w')
        GDATA.write(
            'Title "Formamide  Gaussian 03   B3LYP_cc-pVTZ"\nverbose\nFile %s.fchk\n\nAngstrom\n\n'
            'Multipoles\n  Limit 4\n  Limit 1 H\n  Punch %s.punch\nStart\n\nFinish\n' % (bit[i], bit[i]))
        GDATA.close()
        MDATA = open("DATA%s.MULFIT" % i, 'w')
        MDATA.write(
            'Title\nPokus\n\nVerbose 3\n\nDMA %s.punch\nTotal rank 4\n\n'
            'Ranks\nC 0\nN 0\nO 0\nH 0\n\n' % bit[i])
        MDATA.close()
    for i in range(len(bit)):
        p = subprocess.Popen("gdma <DATA%s.GDMA" % i, shell=True)
        child_processes.append(p)  # start this one, and immediately return to start another
    for cp in child_processes:  # now you can join them together
        cp.wait()
    for i in range(len(bit)):
        subprocess.Popen("mulfit <DATA%s.MULFIT >%s.out" % (i, bit[i]), shell=True).communicate()
        subprocess.Popen("grep Q00 %s.out | cut -c 28- >%s.chg" % (bit[i], bit[i]), shell=True).communicate()

def cli_interface():
    parser = argparse.ArgumentParser(description='Run som.py comparison on two folders with same filenames')
    parser.add_argument('--lig', dest='ligand', type=str, help='ligand (currently 1)', required=True)
    parser.add_argument('--nproc', dest='cpuno', type=str, help='no processors for gaussian', required=False)
    parser.add_argument('--method', dest='cpuno', type=str, help='QM method and basis set', required=False)
    parser.add_argument('--convergence', dest='cpuno', type=str, help='convergence of polarization', required=False)
    args = parser.parse_args()
    return args

#########################################main so far########needs extract to NEUTRAL_fis####################
#def main():
    #args = cli_interface()
    #truedir = args.truth_dir
    #querydir = args.query_dir
lig=str('NEC')
nproc=str('4')
method=str('B3LYP 6-31G*')
convergence=0.001
os.chdir("/home/kevin/NEWT3")
CEM = Molecule('NEUTRAL_fis.pdb')
CEM.read('NEUTRAL_fis.psf')
CEM.filter('protein or resname NEC')
CEM.write('fis.pdb')
CEM.write('fis.psf')
pdb_filter()
#df1 = pd.read_csv('fis1.pdb', sep='\s+',header=None)
#df1.columns = ['atom', 'atno', 'atyp', 'res', 'chain', 'resno', 'X', 'Y', 'Z', 'A', 'B', 'seg', 'ele']
df1=correctpdb()#probably unecessary and introduces problems
df1 = formatdf(df1)#well could do
###LAME FIX######lame fix!!!
for i in range(len(df1)):
    if df1['ele'].values[i]=='HO':
        print (df1['ele'][i])
        df1['ele'].values[i]=re.sub('HO','H',df1['ele'].values[i])
###LAME FIX#######
#df1 = pd.read_csv('fis1.pdb', sep='\s+',header=None)
#df1.columns = ['atom', 'atno', 'atyp', 'res', 'chain', 'resno', 'X', 'Y', 'Z', 'A', 'B', 'seg', 'ele']
bondpairs = filter_psf_return_bonds()
df = pd.DataFrame(bondpairs)
team, bondlength, df2 = read_psf_get_bonds()
df3 = pd.DataFrame(team)
dfle = pd.DataFrame({'bondno':bondlength})
dfwoo=df3.loc[df1.index[df1['res']!=lig]]
dflen=dfle.loc[df1.index[df1['res']!=lig]]
dfwoo.to_csv('bonds',sep=' ',index=False,header=False)
blend=len(dfwoo)
dfmain=psf_pdb_match(df1)
dfbin1 = pd.concat([df1['atno'],df1['atyp'],df1['X'],df1['Y'],df1['Z'],df2['deltaq'],df2['res'],
                   df1['ele'],dflen['bondno']], axis=1, keys=['atno','atyp','X','Y','Z','chg','res','ele','bondno'])
dfbin=dfbin1.loc[dfbin1['res']!=lig]
dflig=dfbin1.loc[dfbin1['res']==lig]
car=(dflig['chg'].sum())
run_gaussian_for_field(dflig, dfbin, nproc, method, lig, car)
field=str('field0')
nudge=str('QMS.com.log')
d=readESP_G(field, nudge)
alpha, polka = polasignMMmae_def(dfbin)
dfalpha = pd.DataFrame({'alpha': alpha})
dfpolka = pd.DataFrame({'polka': polka})
dffin = pd.concat([dfbin['atno'],dfbin['atyp'],dfbin['X'],dfbin['Y'],dfbin['Z'],dfbin['chg'],dfalpha['alpha'],
                   dfpolka['polka']], axis=1, keys =['atno','atyp','X','Y','Z','chg','alpha','pol'])
f_out= open("newt",'w')
for index, row in dffin.iterrows():
    f_out.write(('%5d %-4s%13.6f%13.6f%13.6f %8.5f %6.3f %8.5f'+'\n') %
                (row['atno'],row['atyp'],row['X'],row['Y'],row['Z'],row['chg'],row['alpha'],row['pol']))
f_out.close()
if blend==len(dffin):
    print('YES2')
#return len(dffin)
################################################################
alpha, polka = polasignMMmae_def(dfbin)
sup = len(dffin)#or def
################################################################
subprocess.Popen("paste -d ' ' newt bonds >XYZS", shell=True).communicate()
if d==sup: #have to make len(dffin) output
    print('YES3')
subprocess.Popen("paste -d ' ' field0 XYZS >FCS", shell=True).communicate()
subprocess.Popen("./calc_indpolF FCS", shell=True).communicate()
subprocess.Popen("cp ene2 ene2P0", shell=True).communicate()
enep = enex()
print ("Energy %s\n" % enep)
num=0
cycle(nproc,method,lig,num)
converge(enep, nproc, method, lig)
get_charges()




############################################################
a = 0; b = 0
temp = open("NEUTRAL_fis1.psf", 'w')
bondpairs = [];
with open("NEUTRAL_fis.psf", 'r') as f:
    for line in f:
        if not re.match(r'^\s*$', line):
            if re.search("!NATOM", line):
                a += 1
            if re.search("!NBOND: bonds", line):
                a = 0
            if a > 0:
                if 'NATOM' not in line:
                    temp.write(line)
temp.close()
df4 = pd.read_csv('NEUTRAL_fis1.psf', sep='\s+',header=None)
df4.columns = ['atno', 'seg', 'resno', 'res', 'atyp', 'atyp1', 'deltaq', 'mass', 'restrain']

df5=df4.loc[df2.index]
dfmain=psf_pdb_match(df5)
###################################################################################
directory = "."
names=[]
for file in os.listdir(directory):
    if file.startswith('atomsP'):
        fil=int(re.sub('atomsP','',file))
        names.append(fil)
last = sorted(names)[-1]
file = 'atomsP' + str(last)
###################################################################################
dfaddlig=df5.loc[df1.index[df1['res']==lig]]
dfaddprot=df5.loc[df1.index[df1['res']!=lig]]

poladjust = pd.read_csv(file, sep='\s+',header=None)
poladjust.columns = ['adjust']
ligreplace = pd.read_csv('QMSP.chg', sep='\s+',header=None)
ligreplace.columns = ['QMdeltaq']

#car=(dflig['chg'].sum())
QMformal=ligreplace['QMdeltaq'].sum()

QMform = int('%0.f' % QMformal)
carform = int('%0.f' % car)

if QMform == carform:
    print ('YES5')

if len(dfaddlig) == len(ligreplace):
    for i in range(len(dfaddlig)):
       dfaddlig['deltaq'].values[i]=ligreplace ['QMdeltaq'].values[i]

if len(dfaddprot) == len(poladjust):
    for i in range(len(dfaddprot)):
       polchg=dfaddprot['deltaq'].values[i]+poladjust['adjust'].values[i]
       polcharge = float('%1.6f' % polchg)
       dfaddprot['deltaq'].values[i] = polcharge

ligind = dfaddlig['atno'].values[len(dfaddlig)-1]
protind = dfaddprot['atno'].values[len(dfaddprot)-1]
if protind < ligind:
    result = pd.concat([dfaddprot,dfaddlig])
else:
    result = pd.concat([dfaddlig,dfaddprot])
#result = pd.concat([dfaddprot,dfaddlig])



#####################################################################################
a = 0; b = 0; c=0
temp = open("NEUTRAL_fist.psf", 'w')
bondpairs = [];
with open("NEUTRAL_fis.psf", 'r') as f:
    for line in f:
        if not re.match(r'^\s*$', line):
            if re.search("!NATOM", line):
                a += 1
            if re.search("!NBOND: bonds", line):
                a = 0
            if a > 0:
                if 'NATOM' not in line:
                    dimples=re.split("\s+",line)
                    if int(dimples[1])<=len(result):
                        index=int(dimples[1])-1
                        spine=re.sub(dimples[7],str('%1.6f' % result['deltaq'].values[index]),line)
                        if dimples[7] == str('%1.6f' % result['deltaq'].values[index]):
                            c+=1
                            print (dimples[7])
                        if len(line)==len(spine):
                            b+=1
                            line=spine
                        else:
                            if len(line)>len(spine):
                                b+=1
                                line=re.sub(str(' %1.6f' % result['deltaq'].values[index]),str('  %1.6f' % result['deltaq'].values[index]), spine)
                            else:
                                b+=1
                                line=re.sub(str(' %1.6f' % result['deltaq'].values[index]),str('%1.6f' % result['deltaq'].values[index]), spine)
        temp.write(line)
temp.close()
####################################check bit##########################################
subprocess.Popen("cp NEUTRAL_fis.psf NEUTRAL_fis.psf.old", shell=True).communicate()
subprocess.Popen("cp NEUTRAL_fist.psf NEUTRAL_fis.psf", shell=True).communicate()
#b=total line no 4806
#c=same lines 32
#diff lines = from go below 4774
#4774+32=4806
#diff NEUTRAL_fist.psf NEUTRAL_fis.psf >lines
#grep lines '< ' >go
#wc -l go


######################################################################################
####done not tidied
####################################################################################
#ln src/gdma bin

#/opt/pgi-6.1/linux86-64/6.1/libso

#population Anaylsis ...
#Electrostatic-potential derived charges: Pop=Chelp, ChelpG or MK
#ChelpG = Charges from electrostatic potential (Grid Based)
#MK = Mulliken
#P B3LYP 6-31G* POP=CHELPG nosymm guess=read prop=read charge
#####################################################################################
