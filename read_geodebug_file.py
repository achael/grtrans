import numpy as np

class geodebug:
    def init(self):
        self.file='geodebug.out'
    
    def __init__(self,**kwargs):
        self.init()
        self.__dict__.update(kwargs)

    def read(self):
        f = open(self.file,'r')
        data = f.readlines()
        f.close()

        self.header = np.array(data[0].split(),dtype=float)
        self.lam = np.array(data[1].split(),dtype=float)
        self.time = np.array(data[2].split(),dtype=float)
        self.r = np.array(data[3].split(),dtype=float)
        self.th = np.array(data[4].split(),dtype=float)
        self.phi = np.array(data[5].split(),dtype=float)
        self.kt = np.array(data[6].split(),dtype=float)
        self.kr = np.array(data[7].split(),dtype=float)
        self.kth = np.array(data[8].split(),dtype=float)
        self.kphi = np.array(data[9].split(),dtype=float)
        self.tpr = np.array(data[10].split(),dtype=float)
        self.tpm = np.array(data[11].split(),dtype=float)

        self.rho = np.array(data[12].split(),dtype=float)
        self.p = np.array(data[13].split(),dtype=float)
        self.b = np.array(data[14].split(),dtype=float)
        self.ut = np.array(data[15].split(),dtype=float)
        self.ur = np.array(data[16].split(),dtype=float)
        self.uth = np.array(data[17].split(),dtype=float)
        self.uphi = np.array(data[18].split(),dtype=float)
        self.bt = np.array(data[19].split(),dtype=float)
        self.br = np.array(data[20].split(),dtype=float)
        self.bth = np.array(data[21].split(),dtype=float)
        self.bphi = np.array(data[22].split(),dtype=float)
        self.n = np.array(data[23].split(),dtype=float)
        
        self.t = np.array(data[24].split(),dtype=float)
        self.b = np.array(data[25].split(),dtype=float)
        self.g = np.array(data[26].split(),dtype=float)
        self.incang = np.array(data[27].split(),dtype=float)
        self.ji = np.array(data[28].split(),dtype=float)
        self.ki = np.array(data[29].split(),dtype=float)
        self.jq = np.array(data[30].split(),dtype=float)
        self.jv = np.array(data[31].split(),dtype=float)
        self.kq = np.array(data[32].split(),dtype=float)
        self.kv = np.array(data[33].split(),dtype=float)
        self.rhoq = np.array(data[34].split(),dtype=float)
        self.rhov = np.array(data[35].split(),dtype=float)
        self.ii = np.array(data[36].split(),dtype=float)
        self.iq = np.array(data[37].split(),dtype=float)
        self.iu = np.array(data[38].split(),dtype=float)
        self.iv = np.array(data[39].split(),dtype=float)
        self.taui = np.array(data[40].split(),dtype=float)
        self.ju = np.array(data[41].split(),dtype=float)
        self.ku = np.array(data[42].split(),dtype=float)
        self.rhou = np.array(data[43].split(),dtype=float)
        self.s2xi = np.array(data[44].split(),dtype=float)
        self.c2xi = np.array(data[45].split(),dtype=float)
        self.ang = np.array(data[46].split(),dtype=float)
        self.g2 = np.array(data[47].split(),dtype=float)
        #s2psi = np.array(data[47].split(),dtype=float)
        #c2psi = np.array(data[48].split(),dtype=float)
        self.s2xi2 = np.array(data[48].split(),dtype=float)
        self.c2xi2 = np.array(data[49].split(),dtype=float)
        self.aahat = np.array(data[50].split(),dtype=float)
        self.aat = np.array(data[51].split(),dtype=float)
        self.aar = np.array(data[52].split(),dtype=float)
        self.aath = np.array(data[53].split(),dtype=float)
        self.aaph = np.array(data[54].split(),dtype=float)
        self.kht = np.array(data[55].split(),dtype=float)
        self.khr = np.array(data[56].split(),dtype=float)
        self.khth = np.array(data[57].split(),dtype=float)
        self.khph = np.array(data[58].split(),dtype=float)
        self.bht = np.array(data[59].split(),dtype=float)
        self.bhr = np.array(data[60].split(),dtype=float)
        self.bhth = np.array(data[61].split(),dtype=float)
        self.bhph = np.array(data[62].split(),dtype=float)
