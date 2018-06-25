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
import uncertainties 
import scipy.constants as const
from errorfunkt2tex import error_to_tex
from errorfunkt2tex import scipy_to_unp
from matplotlib.legend_handler import (HandlerLineCollection,HandlerTuple)
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
# params =  uncertainties.correlated_values(params, covar)
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
        plt.xlabel(r'Kanal')
        plt.ylabel(r'$N$')
        plt.legend(loc='best')
        plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
        plt.savefig('build/'+name+'_'+str(rangeVar[0])+'-'+ str(rangeVar[1])+'.pdf') 

def gausFitMitPlot(Werte, ranges, name, plotF=False):
    AllParams = []
    for rangeVar in ranges:
        p0=[np.max(Werte[rangeVar[0]-1:rangeVar[1]])-np.min(Werte[rangeVar[0]-1:rangeVar[1]]),np.min(Werte[rangeVar[0]-1:rangeVar[1]]),(rangeVar[1]-rangeVar[0])/15,(rangeVar[1]+rangeVar[0])/2]
        params, covar = curve_fit(gaus,range(rangeVar[0],rangeVar[1]+1),Werte[rangeVar[0]-1:rangeVar[1]],maxfev=10000,p0=p0)
        AllParams.append(uncertainties.correlated_values(params, covar))
        if plotF:
            plt.cla()
            plt.clf()
            x=np.linspace(rangeVar[0]-0.02*(rangeVar[1]-rangeVar[0]),rangeVar[1]+0.02*(rangeVar[1]-rangeVar[0]),1000)
            plt.plot(range(rangeVar[0],rangeVar[1]+1), Werte[rangeVar[0]-1:rangeVar[1]], 'gx', label='Werte')  
            plt.plot(x, gaus(x,*params), 'r-', label='Fit')
            #plt.plot(x, gaus(x,*p0), 'b-', label='Fit geschätzt')
            plt.xlim(x[0],x[-1]) 
            plt.xlabel(r'Kanal')
            plt.ylabel(r'$N$')
            plt.legend(loc='best')
            plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
            plt.savefig('build/'+name+'_'+str(rangeVar[0])+'-'+ str(rangeVar[1])+'.pdf') 
    
    return np.array(AllParams)


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
sigma=[]
a=[]
for params in EU152Params:
    pos.append(params[3])
    sigma.append(params[2])
    a.append(params[0])

posU=np.array(pos)
pos=unp.nominal_values(pos)
posStd=unp.std_devs(pos)
sigmaU=np.array(sigma)
sigma=unp.nominal_values(sigma)
sigmaStd=unp.std_devs(sigma)
aU=np.array(a)
a=unp.nominal_values(a)
aStd=unp.std_devs(a)


x=np.linspace(1,8192,1000)
params, covar = curve_fit(Line,pos[wahrscheinlichkeiten!=0], unp.nominal_values(energies[wahrscheinlichkeiten!=0]))
umrechnungsParams=uncertainties.correlated_values(params, covar)
print(umrechnungsParams)
plt.cla()
plt.clf()
xA=posU[wahrscheinlichkeiten!=0]
yA=energies[wahrscheinlichkeiten!=0]
plt.errorbar(unp.nominal_values(xA), unp.nominal_values(yA), yerr=unp.std_devs(yA), xerr=unp.std_devs(xA), label='Werte',fmt='x', capthick=0.5, linewidth='0.5',ecolor='b',capsize=1,markersize=1.5) 
plt.plot(x, Line(x, *params), 'r-', label='Fit')
plt.xlim(x[0],x[-1]) 
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

def potenzFunktion(x, a, b, p):
    return a*(x-b)**p

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
params, covar = curve_fit(potenzFunktion, unp.nominal_values(xA), unp.nominal_values(yA),maxfev=10000,sigma=unp.std_devs(yA))
paramsEQU=uncertainties.correlated_values(params, covar)
print(paramsEQU)
plt.cla()
plt.clf()
plt.plot(x, potenzFunktion(x, *params), 'r-', label='Fit') 
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
ranges = [[1635,1660],[400,550]]
print('Cs137')
peakCs137=gausFitMitPlot(Cs137,ranges,'Cs137')
print(peakCs137)
A0=(range(1635,1660+1),Cs137[1635-1:1660])
A1=(range(1642,1644+1),Cs137[1642-1:1644])
A2=(range(1644,1647+1),Cs137[1644-1:1647])
A3=(range(1649,1652+1),Cs137[1649-1:1652])
A4=(range(1652,1654+1),Cs137[1652-1:1654])
print('paramsHalb')
params, covar = curve_fit(Line, *A1,maxfev=10000)
params1=uncertainties.correlated_values(params, covar)
params, covar = curve_fit(Line, *A2,maxfev=10000)
params2=uncertainties.correlated_values(params, covar)
params, covar = curve_fit(Line,*A3,maxfev=10000)
params3=uncertainties.correlated_values(params, covar)
params, covar = curve_fit(Line, *A4,maxfev=10000)
params4=uncertainties.correlated_values(params, covar)

halbeHöhe=np.max(Cs137[1635-1:1660])/2
zehntelHöhe=np.max(Cs137[1635-1:1660])/10
zehntelBreite=-(zehntelHöhe-params1[1])/params1[0]+(zehntelHöhe-params4[1])/params4[0]
halbeBreite=-(halbeHöhe-params2[1])/params2[0]+(halbeHöhe-params3[1])/params3[0]
print('zehntelBreite', zehntelBreite*umrechnungsParams[0])
print('halbeBreite', halbeBreite*umrechnungsParams[0])
print('jk', zehntelBreite/halbeBreite)
print('halbeBreiteBerrechnet', unp.sqrt(np.log(2)*8 *0.1 * Line(peakCs137[0][3],*umrechnungsParams)* 0.0029))

x=np.linspace(1635,1660)
plt.cla()
plt.clf()
mm1, = plt.plot(*A0, 'gx', label='Werte0')  
mm2, = plt.plot(*A1, 'bx', label='Werte1')  
mm3, = plt.plot(*A4, 'yx', label='Werte4') 
mm4, = plt.plot(x, x*0+zehntelHöhe, 'r-', label='Fit')
mm5, = plt.plot(x, Line(x,*unp.nominal_values(params1)), 'b-', label='Fit')
mm6, = plt.plot(x, Line(x,*unp.nominal_values(params4)), 'y-', label='Fit')
plt.ylim(0,max(Cs137[1635-1:1660])/4)
#plt.plot(x, gaus(x,*p0), 'b-', label='Fit geschätzt')
plt.xlabel(r'Kanal')
plt.ylabel(r'$N$')
plt.legend([(mm1, mm2, mm3), mm4, mm5, mm6], ['Wertepaare','Zehntel der Höhe','Fit der linken Flanke','Fit der rechten Flanke'],handler_map={tuple: HandlerTuple(ndivide=None)},loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Cs137Zehntel.pdf') 

plt.cla()
plt.clf()
mm1, = plt.plot(*A0, 'gx', label='Werte0') 
mm2, = plt.plot(*A2, 'yx', label='Werte2') 
mm3, = plt.plot(*A3, 'rx', label='Werte3') 
mm4, = plt.plot(x, x*0+halbeHöhe, 'b-', label='Fit')
mm5, = plt.plot(x, Line(x,*unp.nominal_values(params2)), 'y-', label='Fit')
mm6, = plt.plot(x, Line(x,*unp.nominal_values(params3)), 'r-', label='Fit')
plt.ylim(0,max(Cs137[1635-1:1660])+100)
plt.xlabel(r'Kanal')
plt.ylabel(r'$N$')
plt.legend([(mm1, mm2, mm3), mm4, mm5, mm6], ['Wertepaare','Hälfte der Höhe','Fit der linken Flanke','Fit der rechten Flanke'],handler_map={tuple: HandlerTuple(ndivide=None)},loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Cs137Halb.pdf') 

##########################################################add

def diffWirkung(E,c,h):
    Egamma=Line(unp.nominal_values(peakCs137[0][3]),*unp.nominal_values(umrechnungsParams))
    r=const.elementary_charge/(4*np.pi * const.epsilon_0 * const.electron_mass*const.c**2)
    m0=const.electron_mass*const.c**2 / (1000*const.electron_volt)
    e=Egamma/m0
    th=8/3 * np.pi * r**2
    return 3/8 * th* 1/(m0*e**2) * (2+(E/(Egamma-E))**2 * (1/e**2 + (1-2/e)* (Egamma-E)/Egamma )) * c+h

def Wirkungin(c,h):
    Egamma=Line(unp.nominal_values(peakCs137[0][3]),*unp.nominal_values(umrechnungsParams))
    r=const.elementary_charge/(4*np.pi * const.epsilon_0 * const.electron_mass*const.c**2)
    m0=const.electron_mass*const.c**2 / (1000*const.electron_volt)
    e=Egamma/m0
    th=8/3 * np.pi * r**2
    term1 = 3*Egamma*c*th
    term2 = 2*e*(2+e*(8+e*(11+e)))
    term3 = (1+2*e)**2 * ((e-2)*e-2)*np.log(1+2*e)
    term4 = 8*e**4* (1+2*e)**2 *m0
    term5 = (2*e*Egamma*h)/(1+2*e)
    return term1*(term2+term3)/term4 + term5



rangeVar=[640,1170]
rangeVar=[750,1170]
xA=Line(np.array(range(rangeVar[0],rangeVar[1]+1)),*unp.nominal_values(umrechnungsParams))
yA=Cs137[rangeVar[0]-1:rangeVar[1]]
params, covar = curve_fit(diffWirkung, xA, yA)
paramsKon=uncertainties.correlated_values(params, covar)
print(paramsKon)
rangeVar=[1,1250]
xA2=Line(np.array(range(rangeVar[0],rangeVar[1]+1)),*unp.nominal_values(umrechnungsParams))
yA2=Cs137[rangeVar[0]-1:rangeVar[1]]
x=np.linspace(rangeVar[0],rangeVar[1],1000)
x=Line(x,*unp.nominal_values(umrechnungsParams))
plt.cla()
plt.clf()
plt.plot(x, diffWirkung(x, *params), 'r-', label='Fit') 
plt.plot(xA2, yA2, 'gx', label='Werte') 
plt.plot(xA, yA, 'bx', label='Werte fit')
Egamma=Line(unp.nominal_values(peakCs137[0][3]),*umrechnungsParams)
print('Egamma', Egamma)
print('ERück', Line(unp.nominal_values(peakCs137[1][3]),*umrechnungsParams))
print('Integral', Wirkungin(*paramsKon)/unp.nominal_values(umrechnungsParams[0]))
print('IntegralPeak', peakCs137[0][0]*(np.sqrt(2*np.pi)*peakCs137[0][2]))

m0=const.electron_mass*const.c**2 / (1000*const.electron_volt)
e=Egamma/m0
Emax=unp.nominal_values(Egamma * 2*e/(1+2*e))
plt.plot(np.array([Emax,Emax]), np.array([0,180]), 'y-', label='Komptenkante')  
plt.xlabel(r'$E_\gamma/\si{\kilo\electronvolt}$')
plt.ylabel(r'$N$')
plt.ylim(0,175)
plt.xlim(0,500)
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Cs137Kon.pdf')







########################################################Ba
print('Ba')
ranges = [[1,8192]]
D = np.genfromtxt('scripts/D',unpack=True)
Plot(D,ranges,'D')
ranges = [[200,215],[650,740],[750,770],[880,900],[950,970]]
DParams=gausFitMitPlot(D,ranges,'D',True)
print(DParams)
print('E1', Line(DParams[0][3],*umrechnungsParams))
print('E2', Line(DParams[1][3],*umrechnungsParams))
print('E3', Line(DParams[2][3],*umrechnungsParams))
print('E4', Line(DParams[3][3],*umrechnungsParams))
print('E5', Line(DParams[4][3],*umrechnungsParams))
Pos=[]
Sigma=[]
a=[]
for param in DParams:
    Pos.append(param[3])
    Sigma.append(param[2])
    a.append(param[0])
Pos=np.array(Pos)
Pos=Pos[1:]
Sigma=np.array(Sigma)
Sigma=Sigma[1:]
a=np.array(a)
a=a[1:]
inhalt=a*(np.sqrt(2*np.pi)*Sigma)
t=3669
wahrscheinlichkeitenBa=
print(inhalt/potenzFunktion(Line(Pos,*umrechnungsParams),*paramsEQU))



########################################################?c060?

print('c060')
ranges = [[1,8192]]
D = np.genfromtxt('scripts/unbekannt',unpack=True)
Plot(D,ranges,'unbekannt')
ranges = [[2900,2930],[3280,3340]]
DParams=gausFitMitPlot(D,ranges,'unbekannt')
print(DParams)
print('E1', Line(DParams[0][3],*umrechnungsParams))
print('E2', Line(DParams[1][3],*umrechnungsParams))