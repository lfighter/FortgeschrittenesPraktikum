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
#from sympy import *
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




def gaus(x, a, c,sigma,b):
    return a* np.exp(-(x-b)**2/(2*sigma**2))+c

def Line(x, a, b):
    return a* x+b

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


#######################
#print(gausFitMitPlot(E,[[0,25]],'test'))

#for i in range(len(bs)):
#    print(bs[i])
#    print(sigmas[i])
#    print([np.max(int(bs[i]-sigmas[i]*5)-1,0),int(bs[i]+sigmas[i]*5)+1])
#    print(gausFitMitPlot(E,[[np.max(np.array([int(bs[i]-sigmas[i]*5)-1,1])),np.min(np.array([int(bs[i]+sigmas[i]*5)+1,len(E)]))]],'test'))

#print(gausFitMitPlot(E,[[0,100]],'test'))

#####################################################kali


#             0         0         0        121       244       295       344       367        411         443         678         778         867         964        1005        1085        1112        1212        1299        1408         1457
ranges = [[100,115],[115,126],[205,230],[300,320],[605,625],[730,750],[853,870],[910,930],[1020,1035],[1100,1118],[1700,1730],[1925,1950],[2145,2170],[2370,2425],[2480,2520],[2680,2720],[2740,2790],[3000,3030],[3200,3255],[3460,3540],[3590,3665]]
energies= unp.uarray([0,0,0,121.7817,244.6974,295.9387,344.2785,367.7891,411.1165,443.9606,678.623,778.9045,867.380,964.057,1005.27,1085.837,1112.076,1212.948,1299.142,1408.013,1457.643],[0,0,0,0.0003,0.0008,0.0017,0.0012,0.0020,0.0012,0.0016,0.005,0.0024,0.003,0.005,0.05,0.010,0.003,0.011,0.008,0.003,0.011])
wahrscheinlichkeiten = unp.uarray([0,0,0,107.3,28.39,1.656,100.0,3.232,8.413,10.63,1.777,48.62,15.90,54.57,2.48,38.04,51.40,5.320,6.14,78.48,1.869],[0,0,0,0.4,0.10,0.015,0.6,0.015,0.026,0.03,0.012,0.22,0.09,0.13,0.04,0.10,0.23,0.021,0.03,0.13,0.014])
wahrscheinlichkeiten*= unp.uarray(0.2659,0.0013)/100
EU152 = np.genfromtxt('scripts/EU152',unpack=True)
print('EU152')
Plot(EU152,[[1,8192],[1,2000],[2000,4000],[1800,4000],[4000,8192]],'EU152')
EU152Params=gausFitMitPlot(EU152,ranges,'EU152')
print(EU152Params)
pos=[]
posStd=[]
sigma=[]
sigmaStd=[]
a=[]
aStd=[]
for params in EU152Params:
    pos.append(unp.nominal_values(params[3]))
    posStd.append(unp.std_devs(params[3]))
    sigma.append(unp.nominal_values(params[2]))
    sigmaStd.append(unp.std_devs(params[2]))
    a.append(unp.nominal_values(params[0]))
    aStd.append(unp.std_devs(params[0]))

pos=np.array(pos)
posStd=np.array(posStd)
posU=unp.uarray(pos,posStd)
sigma=np.array(sigma)
sigmaStd=np.array(sigmaStd)
sigmaU=unp.uarray(sigma,sigmaStd)
a=np.array(a)
aStd=np.array(aStd)
aU=unp.uarray(a,aStd)

x=np.linspace(1,8192,1000)
params, error, sigmay = linregress(pos[wahrscheinlichkeiten!=0], unp.nominal_values(energies[wahrscheinlichkeiten!=0]))
umrechnungsParams=unp.uarray(params,error)
print(umrechnungsParams)
plt.cla()
plt.clf()
plt.plot(pos[wahrscheinlichkeiten!=0], unp.nominal_values(energies[wahrscheinlichkeiten!=0]), 'gx', label='Werte') 
plt.plot(x, Line(x, *params), 'b-', label='Fit') 
plt.xlabel(r'Kanal')
plt.ylabel(r'$E_\gamma/\si{\kilo\electronvolt}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/EnergieKali.pdf') 


#def polynom3(x,a,b,c,d,f,g):
#    return g*x**5+f*x**4+d*x**3+a*x**2+b*x+c

#def polynom3(x,a,b,c,d):
#    return d*x**3+a*x**2+b*x+c

#def polynom3(x,a,b,c):
#    return a*x**2+b*x+c

def polynom3(x,a,b,c,d,f):
    return f*x**4+d*x**3+a*x**2+b*x+c

def expC(x,a,b):
    return a*np.exp(-(x/1000-b))

a=7.31 #cm
r=2.25 #cm
print('omega/4pi')
omegaDurch4PI = (1-a/np.sqrt(a**2+r**2))/2
print(omegaDurch4PI)
AnzahlAnTagen=unp.uarray(6462,1) #days
HalbwertsZeit=unp.uarray(4943,5) #days
AktivitätEu=unp.uarray(4130,60) #bq
AktivitätEu=1/2**(AnzahlAnTagen/HalbwertsZeit) * AktivitätEu
print(AktivitätEu)

print(np.sum(EU152)/(AktivitätEu*omegaDurch4PI*4740))
x=np.linspace(0,3700,1000)
inhalt=aU*(np.sqrt(2*np.pi)*sigmaU)
xA=posU[wahrscheinlichkeiten!=0]
yA=inhalt[wahrscheinlichkeiten!=0]/(AktivitätEu*omegaDurch4PI*4740*wahrscheinlichkeiten[wahrscheinlichkeiten!=0])
xA=Line(xA[1:],*unp.nominal_values(umrechnungsParams))
yA=yA[1:]
params, covar = curve_fit(polynom3, unp.nominal_values(xA), unp.nominal_values(yA),sigma=unp.std_devs(yA))
paramsEQU=unp.uarray(params, np.sqrt(np.diag(covar)))
print(paramsEQU)
plt.cla()
plt.clf()
plt.plot(x, polynom3(x, *params), 'r-', label='Fit') 
plt.errorbar(unp.nominal_values(xA), unp.nominal_values(yA), yerr=unp.std_devs(yA), xerr=unp.std_devs(xA), label='Werte',fmt='x', capthick=0.5, linewidth='0.5',ecolor='b',capsize=1,markersize=1.5) 
plt.xlim(150,1500)
plt.ylim(0,0.3)
plt.xlabel(r'$E_\gamma/\si{\kilo\electronvolt}$')
plt.ylabel(r'$Q$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Q.pdf') 

######################################################Cs137

Cs137 = np.genfromtxt('scripts/Cs137',unpack=True)
ranges = [[1635,1660]]
print('Cs137')
peakCs137=gausFitMitPlot(Cs137,ranges,'Cs137')

print(peakCs137)

def diffWirkung(E,c,h):
    Egamma=Line(unp.nominal_values(peakCs137[0][3]),*unp.nominal_values(umrechnungsParams))
    r=const.elementary_charge/(4*np.pi * const.epsilon_0 * const.electron_mass*const.c**2)
    m0=const.electron_mass*const.c**2 / (1000*const.electron_volt)
    e=Egamma/m0
    th=8/3 * np.pi * r**2
    return 3/8 * th* 1/(m0*e**2) * (2+(E/(Egamma-E))**2 * (1/e**2 + (1-2/e)* (Egamma-E)/Egamma )) * c+h

def Wirkungin(E,c,h):
    Egamma=Line(unp.nominal_values(peakCs137[0][3]),*unp.nominal_values(umrechnungsParams))
    r=const.elementary_charge/(4*np.pi * const.epsilon_0 * const.electron_mass*const.c**2)
    m0=const.electron_mass*const.c**2 / (1000*const.electron_volt)
    e=Egamma/m0
    th=8/3 * np.pi * r**2
    return 3*E*th *(E**2 *(((Egamma-E)*(e-2)*e)+Egamma)/(Egamma*(E-Egamma)**2*e**2)+2)/(8*e**2 *m0)* c+h

def Wirkung(E,c,h):
    return Wirkungin(E,c,h)

rangeVar=[640,1170]
rangeVar=[750,1170]
xA=Line(np.array(range(rangeVar[0],rangeVar[1]+1)),*unp.nominal_values(umrechnungsParams))
yA=Cs137[rangeVar[0]-1:rangeVar[1]]
params, covar = curve_fit(diffWirkung, xA, yA)
params2, covar2 = curve_fit(Wirkung, xA, yA)
paramsKon=unp.uarray(params, np.sqrt(np.diag(covar)))
print(paramsKon)
paramsKon2=unp.uarray(params2, np.sqrt(np.diag(covar2)))
print(paramsKon2)
rangeVar=[1,1250]
xA2=Line(np.array(range(rangeVar[0],rangeVar[1]+1)),*unp.nominal_values(umrechnungsParams))
yA2=Cs137[rangeVar[0]-1:rangeVar[1]]
x=np.linspace(rangeVar[0],rangeVar[1],1000)
x=Line(x,*unp.nominal_values(umrechnungsParams))
plt.cla()
plt.clf()
plt.plot(x, Wirkung(x, *params2), 'b-', label='Fit2') 
plt.plot(x, diffWirkung(x, *params), 'r-', label='Fit') 
plt.plot(xA2, yA2, 'gx', label='Werte') 
plt.plot(xA, yA, 'bx', label='Werte fit')
Egamma=Line(unp.nominal_values(peakCs137[0][3]),*unp.nominal_values(umrechnungsParams))
print('Egamma', Egamma)
m0=const.electron_mass*const.c**2 / (1000*const.electron_volt)
e=Egamma/m0
Emax=Egamma * 2*e/(1+2*e)
plt.plot(np.array([Emax,Emax]), np.array([0,100]), 'b-', label='Werte fit')  
plt.xlabel(r'$E_\gamma/\si{\kilo\electronvolt}$')
plt.ylabel(r'$N$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Cs137Kon.pdf')



########################################################Ba
ranges = [[1,8192]]
D = np.genfromtxt('scripts/D',unpack=True)
print(Plot(D,ranges,'D'))


########################################################?c060?
ranges = [[1,8192]]
D = np.genfromtxt('scripts/unbekannt',unpack=True)
print(Plot(D,ranges,'unbekannt'))