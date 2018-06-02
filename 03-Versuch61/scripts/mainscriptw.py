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

#Daten
Winkelpol, Intenspol  = np.genfromtxt('scripts/polarisation.txt',unpack=True)
T00x, T00I  = np.genfromtxt('scripts/T00mode.txt',unpack=True)
T01x, T01I  = np.genfromtxt('scripts/T01mode.txt',unpack=True)
wellenlaenge = np.genfromtxt('scripts/wellenlaenge.txt',unpack=True)
T00I = T00I-100
print("T00x:", T00x)
print("T00I:", T00I)
print("T01x:", T01x)
print("T01I:", T01I)
#T01I = T01I -400
#hier tabellen erzeugen



#alle Angaben in si basiseinheiten, abhängigkeiten
T00I = T00I/(T00I[0])
T01I = T01I/(T01I[0])
wellenlaenge = 0.001*wellenlaenge
wellenlaengeabstandgitter = 5.5*0.01
gitterkonstante  = 80 *1000
#Stabilitätsprüfung

#def gi(L,r):
#	return 1-L/r


#L=np.linspace(0,10,1000)
#plt.cla()
#plt.clf()
#plt.plot(L, gi(L,1)*gi(L,1), 'g-', label='1')
#plt.plot(L, gi(L,3)*gi(L,5), 'r-', label='3,5')
#plt.plot(L, gi(L,4)*gi(L,5), 'b-', label='4,5')
#plt.plot(L, 1+0*gi(L,3)*gi(L,5), 'y-', label='1konst')
#plt.plot(L, 0+0*gi(L,3)*gi(L,5), 'k-', label='0konst')
#plt.ylim(-0.2, 1.2)
# plt.xlim(0, t[-1]*100)
# plt.xlabel(r'$v/\si{\centi\meter\per\second}$')
# plt.ylabel(r'$\Delta f / \si{\hertz}$')
#plt.legend(loc='best')
#plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
#plt.savefig('build/'+'unnötigerKram')



#T00 mode fit
def T00(x,a,b,c):
	return a*np.exp(-2*((x-c)**2)/(b**2))
	

params, covariance_matrix = curve_fit(T00,T00x,T00I,p0 = [1500,1000,5200])
errors = np.sqrt(np.diag(covariance_matrix))
print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])
print('c =', params[2], '±', errors[2])
#print('d =', params[3], '±', errors[3])

#der plot
x = np.linspace(0, T00x[-1],10000)
plt.cla()
plt.clf()
plt.plot(x, T00(x,params[0],params[1],params[2]), '-', label='Daten mit Bewegungsrichtung aufs Mikrofon zu')
plt.plot(T00x, T00I, 'rx', label='Daten mit Bewegungsrichtung vom Mikrofon weg')
#plt.ylim(0, line(t[-1], *params)+0.1)
#plt.xlim(0, t[-1]*100)
#plt.xlabel(r'$v/\si{\centi\meter\per\second}$')
#plt.ylabel(r'$\Delta f / \si{\hertz}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/'+'T00')

def T01(x,a,b,c):
	return ((x-c)**2)*a*np.exp(-2*((x-c)**2)/(b**2))

params, covariance_matrix = curve_fit(T01,T01x,T01I,p0 = [40,1000,5200])
errors = np.sqrt(np.diag(covariance_matrix))
print('a =', params[0], '±', errors[0])
print('b =', params[1], '±', errors[1])
print('c =', params[2], '±', errors[2])
#print('d =', params[3], '±', errors[3])

#der plot
x = np.linspace(0, T01x[-1],10000)
plt.cla()
plt.clf()
plt.plot(x, T01(x,params[0],params[1],params[2]), '-', label='Daten mit Bewegungsrichtung aufs Mikrofon zu')
plt.plot(T01x, T01I, 'rx', label='Daten mit Bewegungsrichtung vom Mikrofon weg')
#plt.ylim(0, line(t[-1], *params)+0.1)
#plt.xlim(0, t[-1]*100)
#plt.xlabel(r'$v/\si{\centi\meter\per\second}$')
#plt.ylabel(r'$\Delta f / \si{\hertz}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/'+'T01')

#Wellenlängenbestimmung


#def welle(x,c,g,b):


pos = np.genfromtxt('scripts/wellenlaenge.txt', unpack=True)

min = np.linspace(0,0,5)

g=0.0125
b=55

def wellenlaenge(x,n):
	return g*np.sin(np.arctan(np.abs(x)/b))/n

def fitFunkt(x,pos0):
	index0 = np.where(np.abs(x-pos0)==np.min(np.abs(x-pos0)))
	indexall = np.where(np.abs(x-pos0)>np.min(np.abs(x-pos0)))
	indexall=np.array(indexall[0])
	index0=np.array(index0[0])
	n = np.abs(indexall-index0)
	print(n)
	wellenlaenge2=wellenlaenge((x-pos0)[np.abs(x-pos0)>np.min(np.abs(x-pos0))],n)
	return wellenlaenge2-np.mean(wellenlaenge2)

def ErgebnisFunkt(x,pos0):
	index0 = np.where(np.abs(x-pos0)==np.min(np.abs(x-pos0)))
	indexall = np.where(np.abs(x-pos0)>np.min(np.abs(x-pos0)))
	indexall=np.array(indexall[0])
	index0=np.array(index0[0])
	n = np.abs(indexall-index0)
	wellenlaenge2=wellenlaenge((x-pos0)[np.abs(x-pos0)>np.min(np.abs(x-pos0))],n)
	return wellenlaenge2

params, covar = curve_fit(fitFunkt,pos,min,maxfev=10000,bounds=(0,1))
print(unp.uarray(params, np.sqrt(np.diag(covar))))
print(np.mean(ErgebnisFunkt(pos,*params)*10**6))

plt.cla()
plt.clf()
plt.plot((pos)[np.abs(pos-params[0])>np.min(np.abs(pos-params[0]))], ErgebnisFunkt(pos,*params)*10**6, 'g-', label='1')

# plt.xlim(0, t[-1]*100)
# plt.xlabel(r'$v/\si{\centi\meter\per\second}$')
# plt.ylabel(r'$\Delta f / \si{\hertz}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/'+'unnötigerKram2')





