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