from table import makeTable
from table import makeNewTable
from linregress import linregress
from customFormatting import *
from bereich import bereich
from weightedavgandsem import weighted_avg_and_sem
from weightedavgandsem import avg_and_sem
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties import ufloat
import scipy.constants as const
from errorfunkt2tex import error_to_tex
from errorfunkt2tex import scipy_to_unp
from sympy import *
import random
# BackwardsVNominal = []
# BackwardsVStd = []
# for value in BackwardsV:
#     BackwardsVNominal.append(unp.nominal_values(value))
#     BackwardsVStd.append(unp.std_devs(value))
# BackwardsVNominal = np.array(BackwardsVNominal)
# BackwardsVStd = np.array(BackwardsVStd)

# einfacher:
# BackwardsVNominal = unp.nominal_values(BackwardsV)
# BackwardsVStd = unp.std_devs(BackwardsV)

# makeTable([Gaenge, ForwardsVNominal, ForwardsVStd, ], r'{Gang} & \multicolumn{2}{c}{$v_\text{v}/\si[per-mode=reciprocal]{\centi\meter\per\second}$} & ', 'name', ['S[table-format=2.0]', 'S[table-format=2.3]', ' @{${}\pm{}$} S[table-format=1.3]', ], ["%2.0f", "%2.3f", "%2.3f",])

#[per-mode=reciprocal],[table-format=2.3,table-figures-uncertainty=1]

# unp.uarray(np.mean(), stats.sem())
# unp.uarray(*avg_and_sem(values)))
# unp.uarray(*weighted_avg_and_sem(unp.nominal_values(bneuDiff), 1/unp.std_devs(bneuDiff))) achtung sum(gewichte muss gleich anzahl der Messungen sein)

# plt.cla()
# plt.clf()
# plt.plot(ForwardsVNominal*100, DeltaVForwardsNominal, 'gx', label='Daten mit Bewegungsrichtung aufs Mikrofon zu')
# plt.plot(BackwardsVNominal*100, DeltaVBackwardsNominal, 'rx', label='Daten mit Bewegungsrichtung vom Mikrofon weg')
# plt.ylim(0, line(t[-1], *params)+0.1)
# plt.xlim(0, t[-1]*100)
# plt.xlabel(r'$v/\si{\centi\meter\per\second}$')
# plt.ylabel(r'$\Delta f / \si{\hertz}$')
# plt.legend(loc='best')
# plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.savefig('build/'+'VgegenDeltaV')

# a = unp.uarray(params[0], np.sqrt(covar[0][0]))
# params = unp.uarray(params, np.sqrt(np.diag(covar)))
# makeNewTable([convert((r'$c_\text{1}$',r'$c_\text{2}$',r'$T_{\text{A}1}$',r'$T_{\text{A}2}$',r'$\alpha$',r'$D_1$',r'$D_2$',r'$A_1$',r'$A_2$',r'$A_3$',r'$A_4$'),strFormat),convert(np.array([paramsGes2[0],paramsGes1[0],deltat2*10**6,deltat1*10**6,-paramsDaempfung[0]*2,4.48*10**-6 *paramsGes1[0]/2*10**3, 7.26*10**-6 *paramsGes1[0]/2*10**3, (VierteMessung-2*deltat2*10**6)[0]*10**-6 *1410 /2*10**3, unp.uarray((VierteMessung[1]-VierteMessung[0])*10**-6 *1410 /2*10**3, 0), unp.uarray((VierteMessung[2]-VierteMessung[1])*10**-6 *2500 /2*10**3, 0),unp.uarray((VierteMessung[3]-VierteMessung[2])*10**-6 *1410 /2*10**3, 0)]),unpFormat,[[r'\meter\per\second',"",True],[r'\meter\per\second',"",True],[r'\micro\second',"",True],[r'\micro\second',"",True],[r'\per\meter',"",True],[r'\milli\meter',"",True],[r'\milli\meter',"",True],[r'\milli\meter',"",True],[r'\milli\meter',r'1.3f',True],[r'\milli\meter',r'1.3f',True],[r'\milli\meter',r'2.2f',True]]),convert(np.array([2730,2730]),floatFormat,[r'\meter\per\second','1.0f',True])+convert((r'-',r'-'),strFormat)+convert(unp.uarray([57,6.05,9.9],[2.5,0,0]),unpFormat,[[r'\per\meter',"",True],[r'\milli\meter',r'1.2f',True],[r'\milli\meter',r'1.2f',True]])+convert((r'-',r'-',r'-',r'-'),strFormat),convert(np.array([(2730-paramsGes2[0])/2730*100,(2730-paramsGes1[0])/2730*100]),unpFormat,[r'\percent','',True])+convert((r'-',r'-'),strFormat)+convert(np.array([(-paramsDaempfung[0]*2-unp.uarray(57,2.5))/unp.uarray(57,2.5)*100,(4.48*10**-6 *paramsGes1[0]/2*10**3-6.05)/6.05*100, (-7.26*10**-6 *paramsGes1[0]/2*10**3+9.90)/9.90*100]),unpFormat,[r'\percent','',True])+convert((r'-',r'-',r'-',r'-'),strFormat)],r'{Wert}&{gemessen}&{Literaturwert\cite{cAcryl},\cite{alphaAcryl}}&{Abweichung}','Ergebnisse', ['c ','c',r'c','c'])

#A, B, C = symbols('A B C')
#f = A**3 *B*cos(C)
#f2 = scipy_to_unp(f, [A, B, C])
#AW, BW = unp.uarray([1,2],[0.1,0.2])
#CW = 3
#print(f2(AW, BW, CW))
#print(error_to_tex(f,'f',[AW, BW, CW], [A, B, C],[A, B]))

a=7.31 #cm
r=2.25 #cm
omegaDurch4PI = (1-a/np.sqrt(a**2+r**2))/2
print(omegaDurch4PI)
AnzahlAnTagen=unp.uarray(6462,1) #days
HalbwertsZeit=unp.uarray(4943,5) #days
AktivitätEu=unp.uarray(4130,60) #bq
AktivitätEu=1/2**(AnzahlAnTagen/HalbwertsZeit) * AktivitätEu
print(AktivitätEu)

#Daten
#A = np.genfromtxt('scripts/a.txt',unpack=True)
#B = np.genfromtxt('scripts/b.txt',unpack=True)
#D = np.genfromtxt('scripts/d.txt',unpack=True)
#E = np.genfromtxt('scripts/e.txt',unpack=True)


def gaus(x, a, c,sigma,b):
    return a* np.exp(-(x-b)**2/(2*sigma**2))+c

xWerte = range(0,100)
E = np.array(xWerte)+0.1
E*=0
bs=[]
sigmas=[]
for number in range(0,4):
    sigma=(1+random.random())
    sigmas.append(sigma)
    a=random.random()*3+1
    b=random.random()*xWerte[-1]+xWerte[0]
    bs.append(b)
    E2=[]
    for i in xWerte:
        E2.append((random.random()+8)/9*gaus(i,a,0,sigma,b))
    E2=np.array(E2)
    E+=E2
sigmas=np.array(sigmas)
bs=np.array(bs)

def Plot(Werte, ranges, name):
    for rangeVar in ranges:
        plt.cla()
        plt.clf()
        x=np.linspace(rangeVar[0],rangeVar[1],1000)
        plt.plot(range(rangeVar[0],rangeVar[1]+1), Werte[rangeVar[0]-1:rangeVar[1]], 'gx', label='Werte')  
        #plt.xlabel(r'$v$')
        #plt.ylabel(r'$\Delta f / \si{\hertz}$')
        plt.legend(loc='best')
        plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
        plt.savefig('build/'+name+'_'+str(rangeVar[0])+'-'+ str(rangeVar[1])+'.pdf') 

def gausFitMitPlot(Werte, ranges, name):
    AllParams = []
    AllCovar = []
    for rangeVar in ranges:
        p0=[np.max(Werte[rangeVar[0]-1:rangeVar[1]])-np.min(Werte[rangeVar[0]-1:rangeVar[1]]),np.min(Werte[rangeVar[0]-1:rangeVar[1]]),(rangeVar[1]-rangeVar[0])/15,(rangeVar[1]+rangeVar[0])/2]
        params, covar = curve_fit(gaus,range(rangeVar[0],rangeVar[1]+1),Werte[rangeVar[0]-1:rangeVar[1]],maxfev=10000,p0=p0)
        AllParams.append(params)
        AllCovar.append(np.sqrt(np.diag(covar)))
        plt.cla()
        plt.clf()
        x=np.linspace(rangeVar[0],rangeVar[1],1000)
        plt.plot(range(rangeVar[0],rangeVar[1]+1), Werte[rangeVar[0]-1:rangeVar[1]], 'gx', label='Werte')  
        plt.plot(x, gaus(x,*params), 'r-', label='Fit')
        plt.plot(x, gaus(x,*p0), 'b-', label='Fit geschätzt')
        #plt.xlabel(r'$v$')
        #plt.ylabel(r'$\Delta f / \si{\hertz}$')
        plt.legend(loc='best')
        plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
        plt.savefig('build/'+name+'_'+str(rangeVar[0])+'-'+ str(rangeVar[1])+'.pdf') 
    
    return unp.uarray(AllParams,AllCovar)


#print(gausFitMitPlot(E,[[0,25]],'test'))

#for i in range(len(bs)):
#    print(bs[i])
#    print(sigmas[i])
#    print([np.max(int(bs[i]-sigmas[i]*5)-1,0),int(bs[i]+sigmas[i]*5)+1])
#    print(gausFitMitPlot(E,[[np.max(np.array([int(bs[i]-sigmas[i]*5)-1,1])),np.min(np.array([int(bs[i]+sigmas[i]*5)+1,len(E)]))]],'test'))

#print(gausFitMitPlot(E,[[0,100]],'test'))

ranges = [[100,115],[115,126],[205,230],[300,320],[605,625],[853,870],[910,930],[1020,1035],[1100,1118],[1700,1730],[1925,1950],[2145,2170],[2370,2425],[2480,2520],[2680,2720],[2740,2790],[3000,3030],[3200,3255],[3460,3540],[3605,3650]]
EU152 = np.genfromtxt('scripts/EU152',unpack=True)
print('EU152')
Plot(EU152,[[1,8192],[1,2000],[2000,4000],[4000,8192]],'EU152')
EU152Params=gausFitMitPlot(EU152,ranges,'EU152')
print(EU152Params)



Cs137 = np.genfromtxt('scripts/Cs137',unpack=True)
ranges = [[1635,1660]]
print('Cs137')
print(gausFitMitPlot(Cs137,ranges,'Cs137'))

