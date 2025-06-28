#Ввод библиотек
import math
import numpy as np
import pandas as pd

from Ring.Classes import StarterData


#--------------------------------------------------------------------------------------------------------------------------------
#Ввод начальных данных для примерного расчета
def calculate_all(Var0Data, Var0Powers):
    """
    
    :param Var0Data: 
    :param Var0Powers: 
    :return: 
    """
    r = Var0Data.radius
    E = Var0Data.E
    e02 = Var0Data.e02
    delta = Var0Data.thick
    l = 2000
    tay = Var0Data.tay
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
    AR = abs_ring(fi, F)
    np.savetxt('matrix1.txt', AR, fmt='%.2f', delimiter='|')
# --------------------------------------------------------------------------------------------------------------------------------
    def get_minimal_inertia_section(F_):
        M_column = F_[:,1]
        M_max_abs = np.abs(M_column).max()
        Q_column = F_[:,3]
        Q_max_abs = np.abs(Q_column).max()
        alpha = (Q_max_abs/(4.8*E))**(1/3)
        h = (3*M_max_abs/(2*alpha*e02))**(3/7)
        F1min = M_max_abs/(2*h*e02)
        Imin = F1min*h*h
        return Imin
    Imin = get_minimal_inertia_section(AR)
    lambda_Imin = Imin*l/(delta*(r**4))
# --------------------------------------------------------------------------------------------------------------------------------
    m30 = 30
    G = Var0Data.G
    def Cm(m):
        return ((-1)**m)/(1 + lambda_Imin*m**2 * (m*m-1)**2 * (E/(2*G) + m*m*l*l/(6*r*r) ))

    def qt1(i,k):
        qp_ = qp(i,k)
        qt_ = 0
        for m in range(2, m30):
            qt_ += m * Cm(m) * np.sin(m*psi(i,k))
        return qp_ + 1/(np.pi*r)*qt_
    def km1(i,k):
        km_ = Mp(i,k)
        km1_ = 0
        for m in range(2, m30):
            km1_ += (Cm(m)/ (m*m-1) ) * np.cos(m*psi(i,k))
        return km_ - r/(np.pi)*km1_
    def kQ1(i, k):
        kQ_ = Qp(i, k)
        kQ1_ = 0
        for m in range(2, m30):
            kQ1_ += ((m*Cm(m)) / (m * m - 1)) * np.sin(m * psi(i, k))
        return kQ_ + 1/np.pi * kQ1_
    def kN1(i, k):
        kN_ = Np(i, k)
        kN1_ = 0
        for m in range(2, m30):
            kN1_ += ((m*m * Cm(m)) / (m * m - 1)) * np.cos(m * psi(i, k))
        return kN_ + 1 / np.pi * kN1_
    # -------------------------------------------------
    def qt2(i,k):
        qT_ = qt(i,k)
        qT1_ = 0
        for m in range(2, m30):
            qT1_ += Cm(m) * np.cos(m*psi(i,k))
        return qT_ - 1/(np.pi*r)*qT1_
    def km2(i,k):
        Mt_ = Mt(i,k)
        Mt1_ = 0
        for m in range(2, m30):
            Mt1_ += (Cm(m)/ (m*(m*m-1)) ) * np.sin(m*psi(i,k))
        return Mt_ - r/(np.pi)*Mt1_
    def kQ2(i, k):
        kQ_ = Qt(i, k)
        kQ1_ = 0
        for m in range(2, m30):
            kQ1_ += ((Cm(m)) / (m * m - 1)) * np.cos(m * psi(i, k))
        return kQ_ - 1/np.pi * kQ1_
    def kN2(i, k):
        kN_ = Nt(i, k)
        kN1_ = 0
        for m in range(2, m30):
            kN1_ += ((m * Cm(m)) / (m * m - 1)) * np.sin(m * psi(i, k))
        return kN_ + 1 / np.pi * kN1_
    # -------------------------------------------------
    def qt3(i,k):
        qH_ = qh(i,k)
        qH1_ = 0
        for m in range(2, m30):
            qH1_ += (m*m - 1) * Cm(m) * np.cos(m*psi(i,k))
        return qH_ + 1/(np.pi*r*r)*qH1_
    def km3(i,k):
        Mh_ = Mh(i,k)
        Mh1_ = 0
        for m in range(2, m30):
            Mh1_ += (Cm(m)/m) * np.sin(m*psi(i,k))
        return Mh_ + 1/(np.pi)*Mh1_
    def kQ3(i, k):
        kQ_ = Qh(i, k)
        kQ1_ = 0
        for m in range(2, m30):
            kQ1_ += Cm(m) * np.cos(m * psi(i, k))
        return kQ_ + 1/(np.pi*r) * kQ1_
    def kN3(i, k):
        kN_ = Nh(i, k)
        kN1_ = 0
        for m in range(2, m30):
            kN1_ += m * Cm(m) * np.sin(m * psi(i, k))
        return kN_ - 1 / (np.pi*r*r) * kN1_
    # -------------------------------------------------
    def ring(fi, F):
        rows_phi = fi.shape[0]
        cols_F = F.shape[1]
        Res = np.zeros((rows_phi, 5))
        fi0 = np.zeros((rows_phi, 2))

        for i in range(rows_phi):
            fi0[i] = np.array([fi[i], 0])

        for i in range(rows_phi):
            Res[i] = np.array([fi0[i, 0], 0, 0, 0, 0])
            for j in range(cols_F):
                P = F[1, j]
                T = F[2, j]
                H = F[3, j]
                omega = F[0, j] * np.pi / 180

                Res[i] += np.array([0, P * km1(omega, fi0[i, 0] * np.pi / 180), P * kN1(omega, fi0[i, 0] * np.pi / 180),
                                    P * kQ1(omega, fi0[i, 0] * np.pi / 180), P * qt1(omega, fi0[i, 0] * np.pi / 180)])
                Res[i] += np.array([0, T * km2(omega, fi0[i, 0] * np.pi / 180), T * kN2(omega, fi0[i, 0] * np.pi / 180),
                                    T * kQ2(omega, fi0[i, 0] * np.pi / 180), T * qt2(omega, fi0[i, 0] * np.pi / 180)])
                Res[i] += np.array([0, H * km3(omega, fi0[i, 0] * np.pi / 180), H * kN3(omega, fi0[i, 0] * np.pi / 180),
                                    H * kQ3(omega, fi0[i, 0] * np.pi / 180), H * qt3(omega, fi0[i, 0] * np.pi / 180)])

        return Res
    # -------------------------------------------------
    ARR = ring(fi,F)
    np.savetxt('matrix2.txt', ARR, fmt='%.2f', delimiter='|')
    # ----------------------------------------------------------------------------------------
    df = pd.read_excel('P:\Рабочий стол\Рабочий стол для синхронизации\YandexDisk\Рабочий стол для синхронизации\Модели_Файлы_Расчеты\Python\Ring\Справочник.xlsx')
    df.columns = ['F_mm2', 'H_mm', 'S_mm', 'R_mm', 'r_mm', 'Ix_Iy_mm4', 'x0_mm', 'n_prof', 'prof_code']
    df = df.astype({'F_mm2': float, 'H_mm': float, 'S_mm': float, 'R_mm': float, 'r_mm': float, 'Ix_Iy_mm4': float, 'x0_mm': float, 'n_prof': int, 'prof_code': str})
    sort = df.to_dict('records')
    # -------------------------------------------------
    sort_delta = np.array([
    1, 1.2, 1.5, 1.6, 1.8, 1.9, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5, 5.5, 6, 6.5, 7])
    sort_h = np.array([
    50, 56, 63, 71, 80, 90, 100, 110, 125, 140, 160, 180, 200, 220, 250, 280, 320])
    sort_d = np.array([
    1.4, 1.6, 2, 2.6, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 18, 20])
    # -------------------------------------------------
    def get_sort_delta(delta):
        for k in range(len(sort_delta)):
            if sort_delta[k] > delta:
                delta_ = sort_delta[k]
                return delta_
            else: continue
    def get_sort_h(h):
        for k in range(len(sort_h)):
            if sort_h[k] > h:
                h_ = sort_h[k]
                return h_
            else: continue

    def get_sort_d(d):
        for k in range(len(sort_d)):
            if sort_d[k] > d:
                d_ = sort_d[k]
                return d_
            else: continue

    # -------------------------------------------------
    t1 = 25
    t2 = 25

    epc = Var0Data.epc
    eB = Var0Data.eB

    def podbor(fi, F_):
        M_column = F_[:, 1]
        M_max_abs = np.abs(M_column).max()
        Q_column = F_[:, 3]
        Q_max_abs = np.abs(Q_column).max()
        qt_column = F_[:, 4]
        qt_max_abs = np.abs(qt_column).max()

        alpha = (abs(Q_max_abs)/(4.8*E))**(0.3333)
        h = ((3*M_max_abs)/(2*alpha*e02))**(0.428571)
        Fmin = M_max_abs / (2*h*e02)

        P1 = Q_max_abs * t1/(2*h)
        P2 = qt_max_abs * t2/(2)

        Var = 0

        n_rows = len(sort)
        n_cols = len(sort[0]) if n_rows > 0 else 0
        Res = np.zeros((n_rows, n_cols), dtype=object)
        for i in range(n_rows):
            if sort[i]['F_mm2'] < Fmin: continue
            F = sort[i]['F_mm2']
            H = sort[i]['H_mm']
            S = sort[i]['S_mm']
            Ix = sort[i]['Ix_Iy_mm4']
            y0 = sort[i]['x0_mm']
            Name = sort[i]['prof_code']
            ecr0 = 0.385*E*((S/H)**2)

            if ecr0 > epc:
                ez = e02 + ((e02-epc)**2)/(0.002*E)
                s = (1 + (4*ez*ecr0) / ((ez-epc)**2))**(0.5)
                ecr = ez*(s-1)/(s+1)
            else: ecr = ecr0

            hcur = M_max_abs/(2*F*ecr)

            if hcur < (2*H): continue

            hcur = get_sort_h(hcur)
            deltacur = alpha * (hcur**(0.3333))
            deltacur = get_sort_delta(deltacur)
            d1cur = 4*P1*(1.73205) / (np.pi * eB)
            d1cur = get_sort_d(d1cur)
            d2cur = 4*P2*(1.73205) / (np.pi * eB)
            d2cur = get_sort_d(d2cur)

            Res[Var][1] = F
            Res[Var][2] = Name
            Res[Var][3] = hcur.item()
            Res[Var][4] = deltacur.item()
            Res[Var][5] = d1cur.item()
            Res[Var][6] = d2cur.item()
            Res[Var][7] = Ix
            Res[Var][8] = y0

            Var = Var + 1
        return Res
    Podbor = podbor(fi, ARR)
    np.savetxt('matrixPodbor.txt', Podbor, fmt='%s', delimiter='|')
    # ----------------------------------------------------------------------------------------

    def ver_cal(F_):
        M_column = F_[:, 1]
        M_max_abs = np.abs(M_column).max()
        Q_column = F_[:, 3]
        Q_max_abs = np.abs(Q_column).max()
        qt_column = F_[:, 4]
        qt_max_abs = np.abs(qt_column).max()
        N1_column = F_[:, 2]

        Name_ = Podbor[:, 8]

        M2_min_ring = (M_column).min()
        M2_max_ring = (M_column).max()

        N2 = np.where(M_column[:] == M_column.max())[0]               #M_column.where(M_column.max)
        print(N2)
        k = abs(N2)
        N2 = N1_column[k]
        N3 = N1_column[0]
        alpha = (abs(Q_max_abs)/(4.8*E))**(0.3333)

        n_rows = len(Podbor)
        n_cols = len(Podbor[0]) if n_rows > 0 else 0
        Res = np.zeros((n_rows, n_cols), dtype=object)
        Var = 0
        for i in range(16):
            F1 = Podbor[i][1]
            I1 = Podbor[i][7]
            y0 = Podbor[i][8]
            Name = Podbor[i][2]
            h = Podbor[i][3]
            delta0 = Podbor[i][4]
            d1 = Podbor[i][5]
            d2 = Podbor[i][6]

            #name1 = np.where(sort[:, 8] == Podbor[i, 1])[0]                             Name_.index(Podbor[i][1])
            nametarget = Podbor[i][2]
            print(nametarget)
            name1 = next(
                (idx for idx, row in enumerate(sort)
                 if row.get("prof_code") == nametarget # Используем ключ словаря
                 ), None)
            g = abs(name1)
            H = sort[g]['H_mm']
            S = sort[g]['S_mm']
            ecr0 = 0.385 * E * ((S / H) ** 2)
            ez = e02 + ((e02 - epc) ** 2) / (0.002 * E)
            s = (1 + (4 * ez * ecr0) / ((ez - epc) ** 2)) ** (0.5)
            ecr = ez * (s - 1) / (s + 1)
            h1 = h - 2*((H-2*s)/2 + s)
            b0 = 30*delta + delta0 + H + S
            F = h*delta0 + 4*F1 + b0 * delta
            Sx = b0*delta*(h/2 + delta/2) - 2*d2 * (delta + S) * (h/2 + delta - (delta + S)/2)
            c = Sx/F

            Ix = delta0*(h**3)/12
            Ix = Ix + 4*(I1 + F1*(h/2 - y0)**2)
            Ix = Ix + b0*(delta**3)/12
            Ix = Ix + b0*delta0*(h/2 + delta/2)**2
            Ix = Ix - 2 * ((delta0 + 2*S)*(d1**3)/12 + d1*(delta0+2*S) * ((h1/2)**2) - 2)
            Ix = Ix - 2* (d2*((delta+S)**3)/12 + d2*(delta+S) * (h/2 - (S- delta)/2)**2)
            I = Ix - F*c**2
            eA1 = M2_max_ring/I * (h/2-c) + N2/F
            eB1 = -M2_max_ring/I * (h/2+c) + N2/F
            eA2 = M2_min_ring/I * (h/2-c) + N3/F
            eB2 = -M2_min_ring/I * (h/2+c) + N3/F

            eрастmax = max(eA1, eB1, eA2, eB2)
            eсжmax = min(eA1, eB1, eA2, eB2)

            nраст = 0.8*eB/eрастmax
            nсжим = abs(ecr/eсжmax)

            S0 =  1/2 * delta0 * (h/2-c) * (h/2-c) + b0*delta*(h/2-c+delta/2) + 2*F1 * (h/2 - c - y0) - d1*(delta0 + 2*S)*(h1/2 - c) - 2*d2*(delta+S) * (h/2 - c -(s-delta)/2)
            taymax = Q_max_abs*S0 / (I*delta0)
            taykr0 = (4.8 * E * (delta0**2))/(h1**2)
            taypc = epc/1.73205
            tay0 = ez/1.73205
            taykr = tay0*taykr0/ (taykr0 + tay0 + taypc)
            nуст = taykr/taymax

            S1 = b0 * delta * (h/2 - c + delta/2) + 2*F1 * (h/2 - c -y0) - d1*s*S*(h1/2 - c) - 2*d2*(delta + S) * (h/2 - c -(s-delta)/2)
            P1 = Q_max_abs * S1 * t1 / I
            P1ср = np.pi*d1**2 *tay / 4
            nср = 2*P1ср / P1
            P2 = qt_max_abs * t2/2
            P2ср = np.pi*d2**2 *tay / 4
            nзп = P2ср / P2

            ecm1 = P1/(d1*delta0)
            ecm2 = 1.3*eB
            ncm = ecm2/ecm1
            ecm = P2/(d2*delta)
            nпрсм = ecm/ncm

            Res[Var][0] = F1
            Res[Var][1] = Name
            Res[Var][2] = nраст
            Res[Var][3] = nсжим
            Res[Var][4] = nуст
            Res[Var][5] = nср
            Res[Var][6] = nзп
            Res[Var][7] = ncm
            Res[Var][8] = nпрсм
            Var = Var + 1
        return Res
    ARRR = ver_cal(ARR)
    np.savetxt('matrix3.txt', ARRR, fmt='%s', delimiter='|')
