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

def line(x,a,b):
	return a*x+b

CountsKallibrierung=np.genfromtxt('scripts/Kalibrierung.txt')
print('Anzahl an Stopereignissen:',np.sum(CountsKallibrierung))
print('Anzahl an Stopereignissen (Impulszähler):', 17880)

PeakPos = []
PeakPosStd=[]
for i in range(0,512):
	if CountsKallibrierung[i] != 0:
		if CountsKallibrierung[i-1] != 0:
			PeakPos.pop()
			PeakPosStd.pop()
			nom,std=weighted_avg_and_sem([i-1,i],[CountsKallibrierung[i-1],CountsKallibrierung[i]])
			PeakPos.append(nom)
			PeakPosStd.append(std)
			#print('PeakPos1',i,':', PeakPos[-1])
			#print('PeakPos2',i,':',unp.uarray(*weighted_avg_and_sem([i-1,i],[CountsKallibrierung[i-1],CountsKallibrierung[i]])))
		else:
			nom,std=weighted_avg_and_sem([i],[CountsKallibrierung[i]])
			PeakPos.append(nom)
			PeakPosStd.append(std)
PeakPos=unp.uarray(PeakPos,PeakPosStd)
makeNewTable([convert(PeakPos,unpFormat,[r'','',True])],r'\multicolumn{1}{c}{Kanal}','tab1', [r'S'])
#print(PeakPos)
time = np.linspace(0.9,1*len(PeakPos),len(PeakPos))
#print(time)
params,std,sigmay = linregress(unp.nominal_values(PeakPos),time)
uparams=unp.uarray(params,std)
print('a:',uparams[0])
print('b:',uparams[1])

timeOffset=params[1]
PeakPos2=np.linspace(0,unp.nominal_values(PeakPos)[-1]*1.02+1)
plt.cla()
plt.clf()
plt.errorbar(unp.nominal_values(PeakPos), time,fmt='x',xerr=unp.std_devs(PeakPos), label='Messwerte')
plt.plot(PeakPos2, line(PeakPos2,*params), 'r-', label='Fit')
# plt.ylim(0, line(t[-1], *params)+0.1)
plt.xlim(0, PeakPos2[-1])
plt.ylabel(r'$T/\si{\second}$')
plt.xlabel(r'$\text{Kanal}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/'+'LinFit')
EineMicroSekInChan = 1/unp.uarray(params[0],std[0])
print('EineMicroSekInChan: ',EineMicroSekInChan)






def f(x,a,b,c):
	return (np.exp(-a*x+b)+c)


lambdasTheorie = unp.uarray(2.1969811,0.0000022)
gelaufeneZeitInSec = 92164
AnzahlAnEreignissenInsgesammt = unp.uarray(2212943,np.sqrt(2212943))
print('AnzahlAnEreignissenInsgesammt: ',AnzahlAnEreignissenInsgesammt)
Erwartungswert = AnzahlAnEreignissenInsgesammt/gelaufeneZeitInSec*(10**-6)*20
WahrscheinlichkeitnachfolgendesTeilchen = Erwartungswert*unp.exp(-Erwartungswert)
AnzahlUntergrundEreignisse = WahrscheinlichkeitnachfolgendesTeilchen * AnzahlAnEreignissenInsgesammt
AnzahlUntergrundProKanal = AnzahlUntergrundEreignisse /(20*EineMicroSekInChan)
print('Erwartungswert: ', Erwartungswert)
print('WahrscheinlichkeitnachfolgendesTeilchen: ',WahrscheinlichkeitnachfolgendesTeilchen)
print('AnzahlUntergrundEreignisse: ',AnzahlUntergrundEreignisse)
print('AnzahlUntergrundProKanal: ',AnzahlUntergrundProKanal)

AnzahlGestoppt = 17880
print("anzahlChannels:", 20*EineMicroSekInChan)


Counts=np.genfromtxt('scripts/Messwerte.txt')
Counts=unp.uarray(Counts,np.sqrt(Counts))
channel=np.linspace(1/100,512/100,512)
#makeNewTable([convert(channel*100,floatFormat,[r'','1.0f',True]),convert(Counts,unpFormat,[r'','',False])],r'{Kanal}&\multicolumn{2}{c}{Anzahl an Ereignissen}','tab2', ['S[table-format=2.0]', 'S[table-format=2.3]', ' @{${}\pm{}$} S[table-format=1.3]'])
Counts2 = []
channel2 = []
for i in range(0,511,2):
	Counts2.append(Counts[i]+Counts[i+1])
	channel2.append((channel[i]+channel[i+1])/2.0)

	
channel2=np.array(channel2)
Counts2Nominal = []
Counts2Std = []
for value in Counts2:
     Counts2Nominal.append(unp.nominal_values(value))
     Counts2Std.append(unp.std_devs(value))
Counts2=unp.uarray(Counts2Nominal,Counts2Std)



params2, covar2 = curve_fit(f,channel2[9:-37],unp.nominal_values(Counts2[9:-37]),maxfev=100000,sigma=unp.std_devs(Counts2[9:-37]))
fitparams2=unp.uarray(params2, np.sqrt(np.diag(covar2)))

channelplot2 = np.linspace(0,513/100,1000)
plt.cla()
plt.clf()
#plt.plot(channel2[0:9], Counts2[0:9], 'bx', label='Messwerte',linewidth='0.1')
plt.errorbar(channel2[9:-37]*100/EineMicroSekInChan.nominal_value+timeOffset, unp.nominal_values(Counts2[9:-37]),yerr=unp.std_devs(Counts2[9:-37]), label='In den Fit einbezogene Messwerte',fmt='x', capthick=0.5, linewidth='0.5',ecolor='b',capsize=1,markersize=1.5)
plt.plot(channelplot2*100/EineMicroSekInChan.nominal_value+timeOffset, f(channelplot2,*params2), 'r-', label='Fit',linewidth='0.5')
##plt.ylim(0, line(t[-1], *params)+0.1)
plt.xlim(0, 513/EineMicroSekInChan.nominal_value)
plt.yscale('log')
plt.ylabel(r'$N$')
plt.xlabel(r'$T/\si{\second}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/'+'Fit')
lambdas2 = (1/fitparams2[0])/EineMicroSekInChan*100
print('a in 1/us: ', fitparams2[0]*EineMicroSekInChan/100)
print('b: ', fitparams2[1])
print('c: ', fitparams2[2])
print('Mittlere Lebensdauer2 in us: ',lambdas2)
print('Mittlere Lebensdauer relative Abweichung: ',abs(1-lambdas2/lambdasTheorie))
print('Untergrund pro Channel: ', fitparams2[2]/2)
print('Untergrund relative Abweichung: ', abs(1-(fitparams2[2]/2)/AnzahlUntergrundProKanal))




