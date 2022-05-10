#%%
from calendar import c
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

mi = pd.read_csv('cleveland_mi.csv')

# Assume K2O is perfectly imcompatible
mi['F'] = min(mi.K2O) / mi.K2O
mi.sort_values('F',inplace=True)

melt_fraction = 1

model = np.polyfit(mi.F, mi.SiO2, 3)
predict = np.poly1d(model)
plt.scatter(mi.F,mi.SiO2,c='red')
plt.plot(mi.F,predict(mi.F),c='k')
plt.xlabel('F')
plt.ylabel('SiO2')
plt.show()
print('SiO2',list(model))

model = np.polyfit(mi.F, mi.Al2O3, 3)
predict = np.poly1d(model)
plt.scatter(mi.F,mi.Al2O3,c='red')
plt.plot(mi.F,predict(mi.F),c='k')
plt.xlabel('F')
plt.ylabel('Al2O3')
plt.show()
print('Al2O3',list(model))

model = np.polyfit(mi.F, mi.FeO, 3)
predict = np.poly1d(model)
plt.scatter(mi.F,mi.FeO,c='red')
plt.plot(mi.F,predict(mi.F),c='k')
plt.xlabel('F')
plt.ylabel('FeO')
plt.show()
print('FeO',list(model))

model = np.polyfit(mi.F, mi.MgO, 3)
predict = np.poly1d(model)
plt.scatter(mi.F,mi.MgO,c='red')
plt.plot(mi.F,predict(mi.F),c='k')
plt.xlabel('F')
plt.ylabel('MgO')
plt.show()
print('MgO',list(model))

model = np.polyfit(mi.F, mi.CaO, 3)
predict = np.poly1d(model)
plt.scatter(mi.F,mi.CaO,c='red')
plt.plot(mi.F,predict(mi.F),c='k')
plt.xlabel('F')
plt.ylabel('CaO')
plt.show()
print('CaO',list(model))

model = np.polyfit(mi.F, mi.Na2O, 3)
predict = np.poly1d(model)
plt.scatter(mi.F,mi.Na2O,c='red')
plt.plot(mi.F,predict(mi.F),c='k')
plt.xlabel('F')
plt.ylabel('Na2O')
plt.show()
print('Na2O',list(model))

model = np.polyfit(mi.F, mi.K2O, 3)
predict = np.poly1d(model)
plt.scatter(mi.F,mi.K2O,c='red')
plt.plot(mi.F,predict(mi.F),c='k')
plt.xlabel('F')
plt.ylabel('K2O')
plt.show()
print('K2O',list(model))

model = np.polyfit(mi.H2O, mi.K2O, 1)
predict = np.poly1d(model)
plt.scatter(mi.H2O,mi.K2O,c='red')
plt.plot(mi.H2O,predict(mi.H2O),c='k')
plt.xlabel('F')
plt.ylabel('K2O')
plt.show()
print('K2O - linear',list(model))

model = np.polyfit(mi.F, mi.P2O5, 3)
predict = np.poly1d(model)
plt.scatter(mi.F,mi.P2O5,c='red')
plt.plot(mi.F,predict(mi.F),c='k')
plt.xlabel('F')
plt.ylabel('P2O5')
plt.show()
print('P2O5',list(model))

model = np.polyfit(mi.F, mi.MnO, 3)
predict = np.poly1d(model)
plt.scatter(mi.F,mi.MnO,c='red')
plt.plot(mi.F,predict(mi.F),c='k')
plt.xlabel('F')
plt.ylabel('MnO')
plt.show()
print('MnO',list(model))

