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

lambdavac = 632.99 *10**-6 #in mm
L = unp.uarray(100,0.1) #in mm
T = 1 # in mm

def kontrastf(phi,phimax,a):
	return 2*np.abs(np.sin(phi-phimax)*np.cos(phi-phimax))*a
v0=406
winkel, maximum, minimum  = np.genfromtxt('scripts/kontrast.txt',unpack=True)
winkelinrad = 2*np.pi*winkel/360
kontrast = (maximum-minimum)/(maximum+minimum-2*406)
makeNewTable([winkel, maximum, minimum, kontrast],r'{$\phi/\si{\degree}$} & {$U_\text{max}/\si{\milli\volt}$} & {$U_\text{min}/\si{\milli\volt}$} & {$K$}','Kontrast',['S[table-format=2.0]','S[table-format=4.0]','S[table-format=3.0]','S[table-format=1.2]'],['{:1.0f}','{:1.0f}','{:1.0f}','{:1.2f}'])
params, covar = curve_fit(kontrastf,winkelinrad,kontrast)
fitergebniss = unp.uarray(params, np.sqrt(np.diag(covar)))
winkelf=np.linspace(-1,90,1000)
plt.cla()
plt.clf()
plt.plot(winkelf, kontrastf(winkelf*2*np.pi/360,*params), 'b-', label='Fit')
plt.plot(winkel, kontrast, 'rx', label='Messwerte')
print(fitergebniss)
print('Maximum bei:', fitergebniss[0]/(2*np.pi)*360%90+45)
#plt.ylim(0, line(t[-1], *params)+0.1)
plt.xlim(-1, 90)
plt.xlabel(r'$\phi/\si{\degree}$')
plt.ylabel(r'$K$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/'+'kontrast')


#gglllaaaaaaaaaaaaassssssss
def n(a,phi1,phi2):
	return (a*a+2*(unp.cos(phi1)-unp.cos(phi2))*(1-a))/(2*(unp.cos(phi1)-unp.cos(phi2)-a))

def nalt(M, phi):
	return (phi*T)/(phi*T-M*lambdavac)
zehnGradinRad = 10*2*np.pi/360

def MFastExakt(phi,n):
	return 
def MFastExakt(phi,n):
	vorfaktor = T/lambdavac
	ersterTerm = (n-np.cos(phi+zehnGradinRad-np.arcsin(1/n*np.sin(phi+zehnGradinRad))))/np.cos(np.arcsin(1/n*np.sin(phi+zehnGradinRad)))
	zweiterTerm = (n-np.cos(phi-zehnGradinRad-np.arcsin(1/n*np.sin(phi-zehnGradinRad))))/np.cos(np.arcsin(1/n*np.sin(phi-zehnGradinRad)))
	return vorfaktor*(ersterTerm-zweiterTerm)

def MFastExakt2(phi,n):
	vorfaktor = T/lambdavac
	ersterTerm = (n-np.cos(phi+zehnGradinRad-np.arcsin(1/n*np.sin(phi+zehnGradinRad))))/np.cos(phi+zehnGradinRad)
	return vorfaktor*(ersterTerm)

Mglass = np.genfromtxt('scripts/BIglas.txt',unpack=True)
makeNewTable([Mglass],r'{$M$}','Glas',['S[table-format=2.0]'],['{:1.0f}'])
params, covar = curve_fit(MFastExakt,(Mglass*0+zehnGradinRad),Mglass,bounds=(1,3))
fitergebniss = unp.uarray(params, np.sqrt(np.diag(covar)))
MglassU = unp.uarray(*avg_and_sem(Mglass))
print(Mglass)
a=MglassU*lambdavac/(2*T)
N=n(a,0*2*np.pi/360,20*2*np.pi/360)
print(10*2*np.pi/360)
print(N)
print(nalt(MglassU,10*2*np.pi/360 * 20*2*np.pi/360))
print(MFastExakt(0,1.3))
print(fitergebniss)
phi=np.linspace(-zehnGradinRad*6,zehnGradinRad*6,1000000)
plt.cla()
plt.clf()
plt.plot(phi, MFastExakt(phi, *params), 'b-', label='Fit')
plt.plot(phi, MFastExakt(phi, 1), '-', label='n=1')
plt.plot(phi, MFastExakt(phi, 1.2), '-', label='n=1.2')
plt.plot(phi, MFastExakt(phi, 1.5), '-', label='n=1.5')
plt.plot(phi, MFastExakt(phi, 2), '-', label='n=2')
plt.plot(phi, MFastExakt(phi, 2.5), '-', label='n=2.5')
plt.plot((Mglass*0+zehnGradinRad),Mglass, 'rx', label='Messwerte')
plt.ylim(-180, 180)
plt.xlim(-zehnGradinRad*6, zehnGradinRad*6)
plt.xlabel(r'$\phi$')
plt.ylabel(r'$M$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/'+'glas')
#gasssssssssssss
def nhoch2(x,a,b):
	return a*x +b


p, m1, m2, m3 = np.genfromtxt('scripts/BIgasMA.txt',unpack=True)
makeNewTable([p, m1, m2, m3],r'{$p/\si{\milli\bar}$} & {$M_1$} & {$M_2$} & {$M_3$}','Luft',['S[table-format=3.0]','S[table-format=2.0]','S[table-format=2.0]','S[table-format=2.0]'],['{:1.0f}','{:1.0f}','{:1.0f}','{:1.0f}'])
m=m1.tolist()+m2.tolist()+m3.tolist()
m=np.array(m)
p=p.tolist()
p*=3
p=np.array(p)
n2 = unp.nominal_values(m*lambdavac/L  + 1)
#print(n2)
params, error, sigmay= linregress(p,n2**2)
fitparams = unp.uarray(params,error)
print(fitparams)
pres=np.linspace(0,1100,1000)
plt.cla()
plt.clf()
plt.plot(pres, nhoch2(pres, *params), 'b-', label='Fit')
plt.plot(p, n2**2, 'rx', label='Messwerte')
#plt.ylim(1, 1.0003)
plt.xlim(0, 1100)
plt.xlabel(r'$p/\si{\milli\bar}$')
plt.ylabel(r'$n^2$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/'+'luft')
nluft=unp.sqrt(nhoch2(1013,fitparams[0],fitparams[1]))
nlufttheo=unp.uarray(1.000265205,0.000000032 )
print('n_luft:', nluft)
print('n_luftlit:', nlufttheo)
print('n_luftlitabweichung:', unp.sqrt((nlufttheo-nluft)**2)/nlufttheo)










