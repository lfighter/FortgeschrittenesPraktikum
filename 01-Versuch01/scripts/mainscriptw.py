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
# unp.uarray(*weighted_avg_and_sem(unp.nominal_values(bneuDiff), 1/unp.std_devs(bneuDiff)))

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

def f(x,a,b,c):
	return np.exp(-a*x+b)+c

Counts=np.genfromtxt('scripts/Messwerte.txt')
channel=np.linspace(1/100,512/100,512)
params, covar = curve_fit(f,channel[17:],Counts[17:],maxfev=100000,sigma=np.sqrt(Counts[17:])+1)
fitparams=unp.uarray(params, np.sqrt(np.diag(covar)))
print(fitparams)
channelplot = np.linspace(0,513/100,1000)
plt.cla()
plt.clf()
plt.plot(channel[0:17], Counts[0:17], 'bx', label='Im Fit nicht mit einbezogene Messwerte',linewidth='0.1')
plt.plot(channel[17:], Counts[17:], 'g+', label='Im Fit mit einbezogene Messwerte',linewidth='0.1')
plt.plot(channelplot, f(channelplot,*params), 'r-', label='Fit')
# plt.ylim(0, line(t[-1], *params)+0.1)
plt.xlim(0, 513/100)
#plt.ylabel(r'$N\si{\per\second}$')
#plt.xlabel(r'$Channel$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/'+'Fit')
#plt.plot(channel[Counts>0], Counts[Counts>0], 'g+', label='Messwerte',linewidth='0.1')
#plt.plot(channelplot, f(channelplot,*params), 'r-', label='Fit')
#plt.xlim(0, 513)
#plt.yscale('log')
#plt.ylabel(r'$\log(N\si{\per92164\second})$')
#plt.legend(loc='best')
#plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
#plt.savefig('build/'+'FitLin')
lambdas = (1/fitparams[0])/unp.uarray(22,0.5)*100
print('Mittlere Lebensdauer in us: ',lambdas)

Counts2 = []
channel2 = []
for i in range(0,511):
	Counts2.append(Counts[i]+Counts[i+1])
	if(Counts[i]!=0 or Counts[i+1]!=0):
		channel2.append((Counts[i]*channel[i]+Counts[i+1]*channel[i+1])/(Counts[i]+Counts[i+1]))
	else:
		channel2.append((channel[i]+channel[i+1])/2.0)
channel2=np.array(channel2)
Counts2=np.array(Counts2)
params2, covar2 = curve_fit(f,channel2[15:],Counts2[15:],maxfev=100000)
fitparams2=unp.uarray(params2, np.sqrt(np.diag(covar2)))
channelplot2 = np.linspace(0,513/100,1000)
plt.cla()
plt.clf()
plt.plot(channel2[0:15], Counts2[0:15], 'bx', label='Messwerte',linewidth='0.1')
plt.plot(channel2[15:], Counts2[15:], 'g+', label='Messwerte',linewidth='0.1')
plt.plot(channelplot2, f(channelplot2,*params2), 'r-', label='Fit')
# plt.ylim(0, line(t[-1], *params)+0.1)
plt.xlim(0, 513/100)
#plt.ylabel(r'$N\si{\per\second}$')
#plt.xlabel(r'$Channel$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/'+'Fit2')
lambdas2 = (1/fitparams2[0])/22*100
print('Mittlere Lebensdauer in us: ',lambdas2)


#def f2(x,a,b):
#	return a*x+b
#CountsDown=Counts[Counts>0]
#channel2=channel[Counts>0]
#CountsDown=np.log(CountsDown)
#params2, covar2 = curve_fit(f2,channel2[0:100],CountsDown[0:100],maxfev=100000)
#fitparams2=unp.uarray(params2, np.sqrt(np.diag(covar2)))
#print(fitparams2)
#channelplot2 = np.linspace(0,513,1000)
#plt.cla()
#plt.clf()
#plt.plot(channel2, CountsDown, 'g+', label='Messwerte',linewidth='0.1')
#plt.plot(channelplot2, f2(channelplot2,*params2), 'r-', label='Fit')
# plt.ylim(0, line(t[-1], *params)+0.1)
#plt.xlim(0, 513)
# plt.xlabel(r'$v/\si{\centi\meter\per\second}$')
# plt.ylabel(r'$\Delta f / \si{\hertz}$')
#plt.legend(loc='best')
#plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
#plt.savefig('build/'+'Fit2')
#lambdas2 = (1/(fitparams2[0]))/22
#print(lambdas2)

#1us ~ 22 ch 


