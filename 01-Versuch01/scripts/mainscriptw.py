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


def line(x,a,b):
    return a*x+b

CountsKallibrierung=np.genfromtxt('scripts/Kalibrierung.txt')
counts=CountsKallibrierung
peakPosnom=[]
peakPosstd=[]
for i in range(0,len(CountsKallibrierung)):
    if(counts[i]!=0):
        samePeak=[]
        for j in range(i,len(CountsKallibrierung)):
            if counts[j]==0:
                break
            samePeak.append(counts[j])
            counts[j]=0
        Kanaele=np.linspace(i+1,i+len(samePeak),len(samePeak))
        #print(Kanaele)
        #print(np.array(samePeak))
        nom,std=weighted_avg_and_sem(Kanaele,np.array(samePeak))
        peakPosnom.append(nom)
        peakPosstd.append(std)
          
peakPosnom=np.array(peakPosnom)
peakPosstd=np.array(peakPosstd)
peakPos=unp.uarray(peakPosnom,peakPosstd)  
time=np.linspace(0.9,0.9+len(peakPosnom)-1,len(peakPosnom)) #in microsek
makeNewTable([convert(peakPos,unpFormat,[r'','1.2f',True]),time],r'\multicolumn{1}{c}{Kanal} & {T/\si{\micro\second}}','tab1', [r'S', r'S'])
params, error, sigmay = linregress(peakPosnom,time)
fitparamsKal=unp.uarray(params,error)
EineMicroSekInChan = 1/fitparamsKal[0]
print('Parameter Kalibrierung:')
print(fitparamsKal)
print('EineMicroSekInChan: ',EineMicroSekInChan)

x=np.linspace(1,512,1000)
plt.cla()
plt.clf()
plt.errorbar(peakPosnom, time,fmt='x',xerr=peakPosstd, label='Messwerte')
plt.plot(x, line(x,*params), 'r-', label='Fit')
# plt.ylim(0, line(t[-1], *params)+0.1)
plt.xlim(1,x[-1])
plt.ylabel(r'$T/\si{\micro\second}$')
plt.xlabel(r'$\text{Kanal}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/'+'LinFit')




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
CountsZusammen=[]
Channel=[]
for i in range(1,512,2):
    CountsZusammen.append(Counts[i]+Counts[i-1])
    Channel.append(i+0.5)
CountsZusammen=np.array(CountsZusammen)
Channel=np.array(Channel)

def expC(x,a,b,c):
	return np.exp(-a*x+b)+c

params2, covar2 = curve_fit(expC, line(Channel[9:-37], *params), CountsZusammen[9:-37])
params3, covar3 = curve_fit(expC, line(Channel[9:-37], *params), CountsZusammen[9:-37],sigma=np.sqrt(CountsZusammen[9:-37]))
params4, covar4 = params3, covar3
for i in range(1,1000):
    params4, covar4 = curve_fit(expC, line(Channel[9:-37], *params), CountsZusammen[9:-37],sigma=1/np.abs(CountsZusammen[9:-37]-expC(line(Channel[9:-37], *params),*params4)),maxfev=100000)
fitparamsLebens1=unp.uarray(params2,np.sqrt(np.diag(covar2)))
fitparamsLebens2=unp.uarray(params3,np.sqrt(np.diag(covar3)))
fitparamsLebens3=unp.uarray(params4,np.sqrt(np.diag(covar4)))
print('Parameter Lebensdauer:')
print(fitparamsLebens1)
print('Lebensdauer in microsekunden:', 1/fitparamsLebens1[0])
print('Parameter Lebensdauer gewichtet:')
print(fitparamsLebens2)
print('Lebensdauer gewichtet in microsekunden:', 1/fitparamsLebens2[0])
print('Parameter Lebensdauer iterativ gewichtet:')
print(fitparamsLebens3)
print('Lebensdauer iterativ gewichtet in microsekunden:', 1/fitparamsLebens3[0])
print('Anzahl an Werten insgesammt:', np.sum(CountsZusammen))
print('Anzahl an Werten die rausgenommen wurden:', np.sum(CountsZusammen[1:9]))


x=np.linspace(0,(Channel[9:-37])[-1]*1.02,1000)
plt.cla()
plt.clf()
plt.errorbar(line(Channel[9:-37], *params), CountsZusammen[9:-37],yerr=np.sqrt(CountsZusammen[9:-37]), label='In den Fit einbezogene Messwerte',fmt='x', capthick=0.5, linewidth='0.5',ecolor='b',capsize=1,markersize=1.5)
plt.errorbar(line(Channel[1:9], *params), CountsZusammen[1:9],yerr=np.sqrt(CountsZusammen[1:9]), label='Messwerte mit Rauschen',fmt='x', capthick=0.5, linewidth='0.5',ecolor='g',color='darkgreen',capsize=1,markersize=1.5)
plt.plot(line(x, *params), expC(line(x, *params),*params2), 'r-', label='Fit ungewichtet',linewidth='0.5')
plt.plot(line(x, *params), expC(line(x, *params),*params3), 'k-', label='Fit gewichtet',linewidth='0.5')
plt.plot(line(x, *params), expC(line(x, *params),*params4), 'y-', label='Fit iterativ gewichtet',linewidth='0.5')
plt.xlim(0, line(x[-1], *params))
plt.yscale('log')
plt.ylabel(r'$N$')
plt.xlabel(r'$T/\si{\micro\second}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/'+'Fit')

def constante(x,c):
	return x*0+c

koinz=np.genfromtxt('scripts/koinz',unpack=True)
makeNewTable(convert([*koinz],floatFormat,[r'','3.0f',False]),r'{$\vardelta t/\si{\nano\second}$} & {$N$}','tab2', [r'S[table-format=3.0]', r'S[table-format=3.0]'])
koinz1=[koinz[0][3:9],koinz[1][3:9]]
koinz2=[koinz[0][11:17],koinz[1][11:17]]
koinz3=[koinz[0][19:24],koinz[1][19:24]]
params1, error1, sigmay1 = linregress(*koinz1)
params2, error2 = avg_and_sem(koinz2[1])
params3, error3, sigmay3 = linregress(*koinz3)
params1U=unp.uarray(params1,error1)
params2U=unp.uarray(params2,error2)
params3U=unp.uarray(params3,error3)
print('Params1:',params1U)
print('Params2:',params2U)
print('Params3:',params3U)
x1=(params2U/2-params1U[1])/params1U[0]
x2=(params2U/2-params3U[1])/params3U[0]
print('x1:',x1)
print('x2:',x2)
print('deltatk:',x2-x1-40)
x=np.linspace(-30,30,1000)
plt.cla()
plt.clf()
plt.plot(*koinz, 'yx', label='Wertepaare',linewidth='0.5')
plt.plot(*koinz1, 'rx', label='Wertepaare',linewidth='0.5')
plt.plot(*koinz2, 'bx', label='Wertepaare',linewidth='0.5')
plt.plot(*koinz3, 'gx', label='Wertepaare',linewidth='0.5')
plt.plot(x,line(x,*params1), 'r-', label='Fit',linewidth='0.5')
plt.plot(x,constante(x,params2), 'b-', label='Fit',linewidth='0.5')
plt.plot(x,constante(x,params2/2), 'k-', label='Fit',linewidth='0.5')
plt.plot(x,line(x,*params3), 'g-', label='Fit',linewidth='0.5')
plt.xlim(-30, 30)
plt.ylim(0, 300)
#plt.yscale('log')
plt.ylabel(r'$N$')
plt.xlabel(r'$T/\si{\nano\second}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/'+'koinz')
