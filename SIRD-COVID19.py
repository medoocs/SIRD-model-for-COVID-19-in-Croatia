# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 14:20:53 2021

@author: NIKOLA
"""

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from sklearn.metrics import mean_squared_error

def plotPred(sve):
    # Lockdown koji pomičemo
    tren = 17
    
    # Za koliko ga pomičemo
    od = tren - 13
    do = tren + 14
    
    # Sa kojim danom uspoređujemo
    #end = 344 # Zadnji dan u datasetu
    end = 276 # Max zaraženih
    
    sus = []
    inf = []
    rec = []
    ded = []
    
    for d in range(od, do):
        #                         t     beta   gamma  sigma
        KOEFICIJENTI = np.array([[0,    0.7,   0.69,  0.049 ], #28-02-2020 Poziva se na nošenje maski
                                 [3,    0.46,  0.3,   0.01  ], #03-03-2020 Poziv rizičnim grupama da ostanu kod kuće
                                 [d,   0.79,  0.69,  0.044 ], #19-03-2020 Lockdown 1
                                 [72,   0.2,   0.69,  0.01  ], #11-05-2020 Zatvaranje osnovnih škola
                                 [79,   0.73,  0.69,  0.013 ], #18-05-2020 Ograničavanje javnog druženja na 50 osoba
                                 [117,  0.4,   0.3,   0.01  ], #25-06-2020 Obavezne maske unutra
                                 [135,  0.71,  0.69,  0.012 ], #13-07-2020 Ograničavanje privatnih skupova na 50 osoba, zabrana javih okupljanja i restrikcije na privatna okupljanja
                                 [167,  0.6,   0.54,  0.013 ], #14-08-2020 Zatvaranje kafića i restorana
                                 [242,  0.675, 0.64,  0.01  ], #28-10-2020 Obavezne maske svugdje
                                 [273,  0.725, 0.69,  0.018 ], #28-11-2020 Poziva se na rad od doma, zatvaranje kafića, restorana i teretana
                                 [300,  0.73,  0.717, 0.03  ], #25-12-2020 Lockdown 2
                                 [350,  0.73,  0.69,  0.018 ]])
        
        T = KOEFICIJENTI[:, 0]
        BETA = KOEFICIJENTI[:, 1]
        GAMMA = KOEFICIJENTI[:, 2]
        SIGMA = KOEFICIJENTI[:, 3]
        beta = interp1d(T, BETA, kind=0)
        gamma = interp1d(T, GAMMA, kind=0)
        sigma = interp1d(T, SIGMA, kind=0)
        
        def SIRD(y, t):
            S, I, R, D = y
            dSdt = -beta(t) * I * S
            dIdt = (beta(t) * I * S) - (sigma(t) * I + gamma(t) * I)
            dRdt = gamma(t) * I
            dDdt = sigma(t) * I
            return dSdt, dIdt, dRdt, dDdt
        
        days = 345
        y0 = 1.0, 1/pop, 0.0, 0.0
        t = np.linspace(0, days-1, days)
        REZ = odeint(SIRD, y0, t)
        S = REZ[: ,0]
        I = REZ[: ,1]
        R = REZ[: ,2]
        D = REZ[: ,3]
        sus.append(S[end]*pop)
        inf.append(I[end]*pop)
        rec.append(R[end]*pop)
        ded.append(D[end]*pop)
    
    sus = pd.Series(sus)
    inf = pd.Series(inf)
    rec = pd.Series(rec)
    ded = pd.Series(ded)
    
    if sve:
        f, [ax, ax1, ax2, ax3] = plt.subplots(4,1,figsize=(10, 15), sharex=True) 
    else:
        f, ax = plt.subplots(1,1,figsize=(10, 10), sharex=True)
        f, ax1 = plt.subplots(1,1,figsize=(10, 10), sharex=True)
        f, ax2 = plt.subplots(1,1,figsize=(10, 10), sharex=True)
        f, ax3 = plt.subplots(1,1,figsize=(10, 10), sharex=True)
        
    od += 3
    do -= 3
        
    ax.plot(np.arange(od, do), sus.values[3:-3], 'b', label='Susceptible')
    ax.plot(np.arange(od, do), [susceptible[end]]*(do-od), c='k', label='Susceptible-real', ls='--')           
    ax1.plot(np.arange(od, do), inf.values[3:-3], 'r', label='Infected')
    ax1.plot(np.arange(od, do), [infected[end]]*(do-od), c='k', label='Infected-real', ls='--')     
    ax2.plot(np.arange(od, do), rec.values[3:-3], 'g', label='Recovered')
    ax2.plot(np.arange(od, do), [recovered[end]]*(do-od), c='k', label='Recovered-real', ls='--')     
    ax3.plot(np.arange(od, do), ded.values[3:-3], 'k', label='Dead')
    ax3.plot(np.arange(od, do), [dead[end]]*(do-od), c='k', label='Dead-real', ls='--')
        
    ax.title.set_text(f'Susceptible prediction vs real on day {end}')
    ax1.title.set_text(f'Infected prediction vs real on day {end}')
    ax2.title.set_text(f'Recovered prediction vs real on day {end}')
    ax3.title.set_text(f'Dead prediction vs real on day {end}')
            
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    legend1 = ax1.legend()
    legend1.get_frame().set_alpha(0.5)
    legend2 = ax2.legend()
    legend2.get_frame().set_alpha(0.5)
    legend3 = ax3.legend()
    legend3.get_frame().set_alpha(0.5)
    xos = []
    for i in np.arange(od-tren,do-tren):
        if i < 0:
            xos.append(f"LC {i}")
        else:
            xos.append(f"LC +{i}")
                
    plt.setp(ax, xticks=np.arange(od,do), xticklabels = xos)
    ax.tick_params(axis='x', rotation=90)
    plt.setp(ax1, xticks=np.arange(od,do), xticklabels = xos)
    ax1.tick_params(axis='x', rotation=90)
    plt.setp(ax2, xticks=np.arange(od,do), xticklabels = xos)
    ax2.tick_params(axis='x', rotation=90)
    plt.setp(ax3, xticks=np.arange(od,do), xticklabels = xos)
    ax3.tick_params(axis='x', rotation=90)
        
    ax.set(xlabel='Pomak u vremenu', ylabel='Vrijednost')
    ax1.set(xlabel='Pomak u vremenu', ylabel='Vrijednost')
    ax2.set(xlabel='Pomak u vremenu', ylabel='Vrijednost')
    ax3.set(xlabel='Pomak u vremenu', ylabel='Vrijednost')
        
    plt.show()

def plotK(beta, gamma, sigma):
    beta2 = list(beta(np.arange(0,350)))
    gamma2 = list(gamma(np.arange(0,350)))
    sigma2 = list(sigma(np.arange(0,350)))
    f, ax = plt.subplots(1,1,figsize=(10, 10), sharex=True)
    f, ax1 = plt.subplots(1,1,figsize=(10, 10), sharex=True)
    f, ax2 = plt.subplots(1,1,figsize=(10, 10), sharex=True)
    
    ax.plot(np.arange(0, 350),beta2, 'blue', label='BETA')
    ax1.plot(np.arange(0, 350),gamma2, 'blue', label='GAMMA')
    ax2.plot(np.arange(0, 350),sigma2, 'blue', label='SIGMA')
    
    ax.title.set_text('Tablična interpolacija BETA koeficijenta')
    ax1.title.set_text('Tablična interpolacija GAMMA koeficijenta')
    ax2.title.set_text('Tablična interpolacija SIGMA koeficijenta')
    
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    legend1 = ax1.legend()
    legend1.get_frame().set_alpha(0.5)
    legend2 = ax2.legend()
    legend2.get_frame().set_alpha(0.5)
    
    ax.set(xlabel='Vrijeme', ylabel='Vrijednost')
    ax1.set(xlabel='Vrijeme', ylabel='Vrijednost')
    ax2.set(xlabel='Vrijeme', ylabel='Vrijednost')
    
    plt.show()

def plotter(t, S, I, R, D):

    f, [ax, ax1, ax2, ax3] = plt.subplots(4,1,figsize=(10, 10), sharex=True)
    
    ax.plot(t, S, 'b', alpha=0.7, linewidth=2, label='Susceptible')
    ax.plot(susceptible, c='k', ls='--', label='Susceptible-real')

    ax1.plot(t, I, 'r', alpha=0.7, linewidth=2, label='Infected')
    ax1.plot(infected, c='k', ls='--', label='Infected-real')
    
    ax2.plot(t, R, 'g', alpha=0.7, linewidth=2, label='Recovered')
    ax2.plot(recovered, c='k', ls='--', label='Recovered-real')
        
    ax3.plot(t, D, 'k', alpha=0.7, linewidth=2, label='Dead')
    ax3.plot(dead, c='k', ls='--', label='Dead-real')

    ax.title.set_text('SIRD-Model')
    
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    legend1 = ax1.legend()
    legend1.get_frame().set_alpha(0.5)
    legend2 = ax2.legend()
    legend2.get_frame().set_alpha(0.5)
    legend3 = ax3.legend()
    legend3.get_frame().set_alpha(0.5)
    
    plt.setp(ax3, xticks=[3,17,72,79,117,135,167,242,273,300], 
    xticklabels=["03-03-2020","19-03-2020","11-05-2020","18-05-2020","25-06-2020","13-07-2020","14-08-2020","28-10-2020","28-11-2020","25-12-2020"])
    ax3.tick_params(axis='x', rotation=90)
    
    ax.set(ylabel='Vrijednost')
    ax1.set(ylabel='Vrijednost')
    ax2.set(ylabel='Vrijednost')
    ax3.set(xlabel='Pomak u vremenu', ylabel='Vrijednost')
    plt.show()

def SIRD(y, t):
    S, I, R, D = y
    dSdt = -beta(t) * I * S
    dIdt = (beta(t) * I * S) - (sigma(t) * I + gamma(t) * I)
    dRdt = gamma(t) * I
    dDdt = sigma(t) * I
    return dSdt, dIdt, dRdt, dDdt

# Inicijalizacija pravih podataka
pop = 4089636
data = pd.read_csv('download.txt', sep=",", header=0)
date = data["Datum"][::-1]

susceptible = pd.Series(data["SlucajeviHrvatska"][::-1]).rolling(window=7, center=True).mean()
susceptible = pop - susceptible.values[3:348]
susceptible[np.isnan(susceptible)] = 0

infected = pd.Series(np.diff(data["SlucajeviHrvatska"][::-1], n=1)).rolling(window=7, center=True).mean()
infected = infected.values[:347]
infected[np.isnan(infected)] = 0

recovered = pd.Series(data["IzlijeceniHrvatska"]).rolling(window=7, center=True).mean()#.iloc[::-1]
recovered = recovered.values[:2:-1]
recovered[np.isnan(recovered)] = 0

dead = pd.Series(data["UmrliHrvatska"]).rolling(window=7, center=True).mean()
dead = dead.values[:2:-1]
dead[np.isnan(dead)] = 0

#                         t     beta   gamma  sigma
KOEFICIJENTI = np.array([[0,    0.7,   0.69,  0.049 ], #28-02-2020 Poziva se na nošenje maski
                         [3,    0.46,  0.3,   0.01  ], #03-03-2020 Poziv rizičnim grupama da ostanu kod kuće
                         [17,   0.79,  0.69,  0.044 ], #19-03-2020 Lockdown 1
                         [72,   0.2,   0.69,  0.01  ], #11-05-2020 Zatvaranje osnovnih škola
                         [79,   0.73,  0.69,  0.013 ], #18-05-2020 Ograničavanje javnog druženja na 50 osoba
                         [117,  0.4,   0.3,   0.01  ], #25-06-2020 Obavezne maske unutra
                         [135,  0.71,  0.69,  0.012 ], #13-07-2020 Ograničavanje privatnih skupova na 50 osoba, zabrana javih okupljanja i restrikcije na privatna okupljanja
                         [167,  0.6,   0.54,  0.013 ], #14-08-2020 Zatvaranje kafića i restorana
                         [242,  0.675, 0.64,  0.01  ], #28-10-2020 Obavezne maske svugdje
                         [273,  0.725, 0.69,  0.018 ], #28-11-2020 Poziva se na rad od doma, zatvaranje kafića, restorana i teretana
                         [300,  0.73,  0.717, 0.03  ], #25-12-2020 Lockdown 2
                         [350,  0.73,  0.69,  0.018 ]])

T = KOEFICIJENTI[:, 0]
BETA = KOEFICIJENTI[:, 1]
GAMMA = KOEFICIJENTI[:, 2]
SIGMA = KOEFICIJENTI[:, 3]

# Interpolacija koeficijenata
beta = interp1d(T, BETA, kind=0)
gamma = interp1d(T, GAMMA, kind=0)
sigma = interp1d(T, SIGMA, kind=0)

# Plot koeficijenata
plotK(beta, gamma, sigma)

# Računanje SIRD modela
days = 345
y0 = 1.0, 1/pop, 0.0, 0.0
t = np.linspace(0, days-1, days)
REZ = odeint(SIRD, y0, t)
S = REZ[: ,0]
I = REZ[: ,1]
R = REZ[: ,2]
D = REZ[: ,3]

# Plot SIRD modela uz prave podatke
plotter(t, S * pop, I * pop, R * pop, D * pop)

# Računanje RMSE
errS = mean_squared_error(susceptible[:345], S*pop, squared=False)
errI = mean_squared_error(infected[:345], I*pop, squared=False)
errR = mean_squared_error(recovered[:345], R*pop, squared=False)
errD = mean_squared_error(dead[:345], D*pop, squared=False)

print("RMSE susceptible: ", errS)
print("RMSE infected: ", errI)
print("RMSE recovered: ", errR)
print("RMSE dead: ", errD)

# True - sve na jednom grafu, False - pojedinačni grafovi
plotPred(True)












































































