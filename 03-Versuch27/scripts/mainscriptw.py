﻿from table import makeTable
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
I, B  = np.genfromtxt('scripts/Bfeldkali.txt',unpack=True)

#x werte nun in mm
print("I in ampere:", I)
print("B in tesla:", B)
#hier tabellen erzeugen
makeTable([I[0:20],B[0:20]], r'{$ I / \si{\ampere}$} & {$ B/ \si{\tesla}$}','tabIB1' , ['S[table-format=2.1]' , 'S[table-format=1.3]'] ,  ["%2.1f", "%1.3f"])
makeTable([I[20:39],B[20:39]], r'{$ I/ \si{\ampere}$} & {$ B/ \si{\tesla}$}','tabIB2' , ['S[table-format=2.1]' , 'S[table-format=1.3]'] ,  ["%2.1f", "%1.3f"])

#bfeld kalibrierung


def BvonI(x,a,b):
	return a*x+b
	

params, covariance_matrix = curve_fit(BvonI,I,B)
errors = unp.uarray(params, np.sqrt(np.diag(covariance_matrix)))
#errors =  unp.uarray(uncertainties.correlated_values(params, covariance_matrix))
print('Die Parameter der B feld kalibrierung:')
print('a =' , errors[0])
print('b =',  errors[1])
tesla = errors

#der plot
x = np.linspace(0, I[-1],10000)
plt.cla()
plt.clf()
plt.plot(x, BvonI(x,unp.nominal_values(tesla[0]),unp.nominal_values(tesla[1])), '-', label='Die gefittete Kurve')
plt.plot(I, B, 'rx', label='Die Messdaten')
#plt.ylim(0, line(t[-1], *params)+0.1)
#plt.xlim(0, t[-1]*100)
plt.xlabel(r'$I/\si{\ampere}$')
plt.ylabel(r'$B /\si{\tesla} $')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/'+'Bfeldkali')

#konstanten
lambdan = ufloat(643.8*10**(-9),0)
lambdaa = 480*10**(-9)
lambdaDn = 4.89*10**(-11)
lambdaDa = 2.695*10**(-11)
h = const.h
clight = const.c
ub = const.physical_constants["Bohr magneton"][0]


#normaler Zeemaneffekt
ln, rn  = np.genfromtxt('scripts/normalB1.txt',unpack=True)
mn = np.genfromtxt('scripts/normalB0.txt',unpack=True)
#ma1 und ma2 müssen noch in die datei

#berechne sigma_s
sigma_sn = rn - ln
delta_sn =[]
i = 0
while i < len(mn)-1:
	delta_sn.append(mn[i+1]-mn[i])
	i = i+1



print("delta:", delta_sn)

#berechne delta lambda

def sigmalambda(sigma_s,delta_s,deltalambda):
	return 0.5*(sigma_s / delta_s)*deltalambda


sigmalambdan = sigmalambda(delta_sn,sigma_sn,lambdaDn)
print()


#mittlere sigmalambda

def C(x,c):
	return c
x = range(0,10)
params, covariance_matrix = curve_fit(C,x,sigmalambdan)
errors = unp.uarray(params, np.sqrt(np.diag(covariance_matrix)))
sigmalambdanm = errors
#errors2 =  ufloat(uncertainties.correlated_values(params, covariance_matrix))
print('sigmalamdanm:')
print('sigmalambdanm =' , sigmalambdanm)

#hilfsplot der sigmalambdas

plt.cla()
plt.clf()
plt.plot(x, sigmalambdan, 'rx', label='Die Messdaten')
#plt.ylim(0, line(t[-1], *params)+0.1)
#plt.xlim(0, t[-1]*100)
plt.xlabel(r'$I/\si{\ampere}$')
plt.ylabel(r'$B /\si{\tesla} $')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/'+'normal')



#berechne g faktor von nomalem zeeman
g = ((h*clight /lambdan) -(h*clight /(lambdan+sigmalambdanm)))

#1/(ub* BvonI(9.5,tesla[0],tesla[1]))
print('g =' , g)