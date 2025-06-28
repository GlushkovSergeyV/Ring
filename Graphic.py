import numpy as np
from RashetRing import calculate_all
from Classes import StarterData, StarterPowers

#Ввод исходных данных
data = StarterData(
    radius=1000,       # [мм]
    thick=1,           # [мм]
    E=72000,           # [МПа] модуль Юнга
    e02=270,           # параметр материала
    mu=0.3,            # коэффициент Пуассона
    eB=440,            # предел прочности
    epc=190,           # деформация
    tay=245,           # напряжение [МПа]
    tochnost=5         # точность расчета [градусов]
)
powers = StarterPowers(
    RadialPowers=[-30000, 15000],       # [Н] радиальные силы
    RadialAngle=[90, 340],              # [град] углы приложения
    KasatPowers=[25000, 20000, -30000], # [Н] касательные силы
    KasatAngle=[45, 180, 270],          # [град] углы приложения
    MomentPowers=[9500000, 7500000],    # [Н·мм] моменты
    MomentAngle=[180, 340],             # [град] углы приложения
    k=1                                 # коэффициент прочности
)
A = calculate_all(data, powers)
#np.savetxt ('matrix.txt', A, fmt='%.2f', delimiter='|')