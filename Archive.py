#Ввод библиотек
import math
import numpy as np
#--------------------------------------------------------------------------------------------------------------------------------
#Ввод классов с начальными данными
#--------------------------------------------------------------------------------------------------------------------------------
class StarterData:
    def __init__(self, radius:float, thick:float, E:float, e02:float, mu:float, eB:float,epc:float, tay:float, tochnost:int):
        self.radius = radius
        self.thick = thick
        self.E = E
        self.e02 = e02
        self.mu = mu
        self.eB = eB
        self.epc = epc
        self.tay = tay
        self.G = E/(2*(1+mu))
        self.tochnost = tochnost
    def __str__(self):
        return (f"StarterData(radius={self.radius}, thick={self.thick}, "
                f"E={self.E}, et={self.e02}), G={self.G}")
class StarterPowers:
    def __init__(self, RadialPowers: list, RadialAngle: list, KasatPowers: list, KasatAngle: list, MomentPowers: list, MomentAngle: list, k: float):
        self.RadialPowers = [x * k for x in RadialPowers]
        self.RadialAngle = [x * 1 for x in RadialAngle]
        self.KasatPowers = [x * k for x in KasatPowers]
        self.KasatAngle = [x * 1 for x in KasatAngle]
        self.MomentPowers = [x * k for x in MomentPowers]
        self.MomentAngle = [x * 1 for x in MomentAngle]
    def __str__(self):
        return (f"StarterPowers(RadialPowers={self.RadialPowers}, {self.RadialAngle}, {self.KasatPowers}, {self.KasatAngle}, {self.MomentPowers}, {self.MomentAngle})")
#--------------------------------------------------------------------------------------------------------------------------------
#Ввод начальных данных для примерного расчета
#--------------------------------------------------------------------------------------------------------------------------------
"""
Var0Data = StarterData(1000, 1, 72000, 270, 0.3, 440, 190, 245, 5)
Var0Powers = StarterPowers(
    RadialPowers = [-30000, 15000],
    RadialAngle = [90, 340],
    KasatPowers = [25000, 20000, -30000],
    KasatAngle = [45, 180, 270],
    MomentPowers = [9500000, 7500000],
    MomentAngle = [180, 340],
    k=1
)
"""
r = Var0Data.radius
#--------------------------------------------------------------------------------------------------------------------------------
#Создание матриц
#--------------------------------------------------------------------------------------------------------------------------------

#Матрица углов
fi = [0]
tochnost = Var0Data.tochnost
j=0
i=0
while j<360:
    fi.append(fi[i] + tochnost)
    j=j+tochnost
    i=i+1
fi = np.array(fi)

#Матрица сил

rows = fi.shape
rows_int = int(rows[0])
zero_rows = np.zeros((3, rows_int))
F = np.vstack((fi, zero_rows))

#Присваивание сил в матрице

i=0
j=0
#--------------------------------------------------------------------------------------------------------------------------------
L_P=len(Var0Powers.RadialPowers)
L_T=len(Var0Powers.KasatPowers)
L_H=len(Var0Powers.MomentPowers)
#--------------------------------------------------------------------------------------------------------------------------------
while i<rows_int:
    if (Var0Powers.RadialAngle[j]) == (F[0, i]):
        F[1,i] = int(Var0Powers.RadialPowers[j])
        i=0
        if j<L_P-1:
            j=j+1
        else:
            break
    else:
        i=i+1
#--------------------------------------------------------------------------------------------------------------------------------
i=0
j=0
#--------------------------------------------------------------------------------------------------------------------------------
while i<rows_int:
    if (Var0Powers.KasatAngle[j]) == (F[0, i]):
        F[2,i] = int(Var0Powers.KasatPowers[j])
        i=0
        if j<L_T-1:
            j=j+1
        else:
            break
    else:
        i=i+1
#--------------------------------------------------------------------------------------------------------------------------------
i=0
j=0
#--------------------------------------------------------------------------------------------------------------------------------
while i<rows_int:
    if (Var0Powers.MomentAngle[j]) == (F[0, i]):
        F[3,i] = int(Var0Powers.MomentPowers[j])
        i=0
        if j<L_H-1:
            j=j+1
        else:
            break
    else:
        i=i+1
#--------------------------------------------------------------------------------------------------------------------------------
#Предварительный расчет
#--------------------------------------------------------------------------------------------------------------------------------
def psi(omega,fi):
    if fi < omega:
        return fi - omega + np.pi
    else:
        return fi - omega - np.pi
#--------------------------------------------------------------------------------------------------------------------------------
def qp(omega, fi):
    return -1/(np.pi*r)*np.sin(psi(omega,fi))
def Np(omega, fi):
    return 1/(2*np.pi)*(psi(omega,fi)*np.sin(psi(omega,fi))-3/2*np.cos(psi(omega,fi)))
def Qp(omega, fi):
    return -1/(np.pi*2)*(psi(omega,fi)*np.cos(psi(omega,fi))+1/2*np.sin(psi(omega,fi)))
def Mp(omega, fi):
    return r*(1/(2*np.pi)*(1-psi(omega,fi)*np.sin(psi(omega,fi))-1/2*np.cos(psi(omega,fi))))
def qt(omega, fi):
    return 1/(2*np.pi*r)*(2*np.cos(psi(omega, fi))-1)
def Nt(omega, fi):
    return -(1/(2*np.pi)*(1/2*np.sin(psi(omega, fi))+psi(omega, fi)*np.cos(psi(omega,fi))))
def Qt(omega, fi):
    return 1/(2*np.pi)*(1-1/2*np.cos(psi(omega, fi))-psi(omega, fi)*np.sin(psi(omega,fi)))
def Mt(omega, fi):
    return r*(1/(2*np.pi))*(psi(omega,fi)*(1+np.cos(psi(omega,fi)))-3/2*np.sin(psi(omega,fi)))
def qh(omega, fi):
    return -1/(2*np.pi*r*r)
def Nh(omega, fi):
    return 1/(np.pi*r*r)*np.sin(psi(omega,fi))
def Qh(omega, fi):
    return 1/(2*np.pi*r)*(1-2*np.cos(psi(omega,fi)))
def Mh(omega, fi):
    return 1/(2*np.pi)*(psi(omega, fi)-2*np.sin(psi(omega,fi)))
#--------------------------------------------------------------------------------------------------------------------------------
def abs_ring(fi, F):
    rows_phi = fi.shape[0]
    cols_F = F.shape[1]
    Res = np.zeros((rows_phi, 5))
    fi0 = np.zeros((rows_phi,2))

    for i in range(rows_phi):
        fi0[i] = np.array([fi[i], 0])

    for i in range(rows_phi):
        Res[i] = np.array([fi0[i, 0], 0, 0, 0, 0])
        for j in range(cols_F):
            P = F[1, j]
            T = F[2, j]
            H = F[3, j]
            omega = F[0, j] * np.pi / 180

            Res[i] += np.array([0, P*Mp(omega, fi0[i, 0]*np.pi/180), P*Np(omega, fi0[i, 0]*np.pi/180), P * Qp(omega, fi0[i, 0]*np.pi/180), P * qp(omega, fi0[i, 0]*np.pi/180)])
            Res[i] += np.array([0, T*Mt(omega, fi0[i, 0]*np.pi/180), T*Nt(omega, fi0[i, 0]*np.pi/180), T * Qt(omega, fi0[i, 0]*np.pi/180), T * qt(omega, fi0[i, 0]*np.pi/180)])
            Res[i] += np.array([0, H*Mh(omega, fi0[i, 0]*np.pi/180), H*Nh(omega, fi0[i, 0]*np.pi/180), H * Qh(omega, fi0[i, 0]*np.pi/180), H * qh(omega, fi0[i, 0]*np.pi/180)])

    return Res
#--------------------------------------------------------------------------------------------------------------------------------

# Вызов функции
AR = abs_ring(fi, F)
#print(AR)