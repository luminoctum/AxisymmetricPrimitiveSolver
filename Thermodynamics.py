#! /usr/bin/env python2.7
from numpy import *
from matplotlib.mlab import amap
from scipy.optimize import newton, root
from scipy.interpolate import interp1d
from scipy.integrate import quad

'''@package
Because python and c uses row-major indexing (rightmost index varies the fastest)
species are alway the right most index when needed
function Argument shollow the sequence: temperature, pressure, mixing ratio, species
Cheng Li 2014.5.26'''

# Condensates
class Condensates:
    def __init__(self, mmr = 1.E9):
        self.mmr = mmr # maximum mixing ratio
    def __repr__(self): return self.name

class Water(Condensates):
    name = 'Water(H2O)'
    mu, Rvap = 18.E-3, 461.67 # 8.31 / mu
    al, bl = 11.079, -2261.1
    Lv0, dLvdt = 2.5E6, -2.3E3 # Latent heat at 0 degree Celsius, Kirchhoff's equation, Emanual (4.4.4)
    cliq, cpvap = 4188., 2028. # those two value are from wiki water (data_page)
    Lv = 2.38E6 # convenient value for latent heat

class Ammonia(Condensates):
    name = 'Ammonia(NH3)'
    mu, Rvap = 17.E-3, 488.82
    al, bl = 10.201, -1248
    Lv0, dLvdt = 1.4E6, 0.
    cliq, cpvap = 4753., 2062.
    Lv = Lv0

class Helium:
    name = 'Helium(He)'
    mu, Rvap = 4.E-3, 2077.5
    cpvap = 5193.75 # 2.5 * R

class Hydrogen:
    name = 'Hydrogen(H2)'
    mu, Rvap = 2.E-3, 4155.
    cpvap = 14.32E3 # 3.5 * R = 14542.5, from www.engineeringtoolbox.com/specific-heat-capacity-gases-d_159.html

# Moist functions
def SatVaporP(T, sp):
    if isinstance(sp, list) or isinstance(sp, ndarray):
        result = amap(lambda s: 10. ** (s.al + s.bl / T), sp)
        return result.T
    else:
        return 10. ** (sp.al + sp.bl / T)

def MixingR(e, p):
    if sum(e) >= p:
        if isinstance(e, list) or isinstance(e, ndarray):
            return array([inf] * len(e))
        else:
            return inf
    else:
        return e / (p - sum(e))

def SatMixingRScalar(T, p, sp):
    svp = SatVaporP(T, sp)
    mixr = svp / (p - sum(svp))
    mmr = array([s.mmr for s in sp])
    i = 0
    while (any(mixr < 0) or any(mixr > mmr)):
        mixr[i] = sp[i].mmr
        i += 1
        if i == len(sp): break
        mixr[i:] = svp[i:] * (1. + sum(mixr[:i])) / (p - sum(svp[i:]))
    return mixr

def SatMixingR(T, P, sp):
    if not isinstance(T, ndarray):
        return SatMixingRScalar(T, P, sp)
    else:
        dims = T.shape
        n = size(T)
        result = amap(lambda x, y: SatMixingRScalar(x, y, sp), T.reshape(n), P.reshape(n))
        return result.reshape(dims)

def LatentH(sp, temp = 273.15):
    if isinstance(sp, list) or isinstance(sp, ndarray):
        return amap(lambda s: s.Lv, sp)
        #return amap(lambda s: s.Lv0 - s.dLvdt * (temp - 273.15), sp)
    else:
        return sp.Lv

def IsSat(temp, ptol, mixr, sp):
    return less(SatVaporP(temp, sp), mixr / (1 + sum(mixr)) * ptol)

class CloudLayer:
    def __init__(self, temp, ptol, mixr, sp):
        self.temp = temp
        self.ptol = ptol
        self.mixr = mixr
        self.species = sp
        self.name = '%-20s: p = %12.2f Pa, T = %8.2f K' % (sp, ptol, temp)
    def __repr__(self): return self.name

class CondensationSolution:
    def make(self):
        self.nclouds = len(self.cloud)
        self.pmin = min(self.sample_p)
        self.pmax = max(self.sample_p)
    def __repr__(self):
        result = 'Cloud Layers:\n'
        for i in range(self.nclouds): result += str(self.cloud[i]) + '\n'
        result += 'Maximum pressure: %10.2f\n' % self.pmax
        result += 'Minimum pressure: %10.2f\n' % self.pmin
        result += 'Temperature interpolation function: funcT<%10.2f, %10.2f>\n' % (self.pmin, self.pmax)
        result += 'Mixing ratio interpolation function: funcX<%10.2f, %10.2f>\n' % (self.pmin, self.pmax)
        result += 'Species: '
        for i in range(len(self.species)): result += str(self.species[i]) + ', '
        return result

# Atmosphere
class Atmosphere:
    def __init__(self, Tref = 134.8, Pref = 1.E5, grav = 10.44, mu = 2.198E-3, cp = 11456.):
        self.Tref = Tref
        self.Pref = Pref
        self.grav = grav
        self.mu = mu
        self.cpdry = cp
        self.Rdry = 8.31 / mu
        self.Hscale = self.Rdry * Tref / grav

    def eps(self, sp):
        if isinstance(sp, list) or isinstance(sp, ndarray): return amap(lambda s: s.mu / self.mu, sp)
        else: return sp.mu / self.mu

    def PotentialT(self, T, p, species = []):
        # Emanuel 4.5.11
        result = T * (self.Pref / (p - sum(SatVaporP(T, species), axis = 0)))**(self.Rdry / self.cpdry)
        #result = T * (self.Pref / p)**(self.Rdry / self.cpdry)
        if isinstance(species, list) or isinstance(species, ndarray):
            if len(species) == 0: return result
            else: return result * prod(exp(self.eps(species) * LatentH(species) * MixingR(SatVaporP(T, species), p) / (self.cpdry * T)))
        else: return result * exp(self.eps(species) * LatentH(species) * MixingR(SatVaporP(T, species), p) / (self.cpdry * T))

    def VirtualTScalar(self, T, x, s):
        return T * (1. + sum(x)) / (1. + sum(self.eps(s) * x))

    def VirtualT(self, T, X, s):
        if not isinstance(T, ndarray):
            result = self.VirtualTScalar(T, X, s)
            return result
        else:
            # ndarray might be zero shape and looks exactly like a float!
            dims = T.shape
            n = size(T)
            result = amap(lambda x, y: self.VirtualT(x, y, s), T.reshape(n), X.reshape(n, 2))
            return result.reshape(dims)

    def SatVirtualT(self, T, p, s):
        return self.VirtualT(T, SatMixingR(T, p, s), s)

    def Adiabats(self, pt, pmin, pmax, species = [], num = 200, endpoint = True):
        plev = logspace(log10(pmin), log10(pmax), num = num, endpoint = endpoint)
        result = amap(lambda p: 
                newton(lambda v: self.PotentialT(v, p, species) - pt, 
                    pt * (p / self.Pref)**(self.Rdry / self.cpdry)),
                plev)
        return plev, result

    def EqCondensation(self, T, P, X, sp, pmin, pmax, num = 200):
        sp = array(sp)
        XS = array([s.mmr for s in sp])
        mixr = amap(min, X, XS)
        sat = IsSat(T, P, mixr, sp)
        sat0 = sat.copy()
        mixr = logical_not(sat) * mixr + sat * XS
        # mixing ratio is the "real mixing ratio" when it is unsaturated. 
        # But it is the "maximum mixing ratio" if it is saturated.
        plev = logspace(log10(pmin), log10(pmax), num = num)
        # logspace doesn't behave too well in python, need to force boundary condition
        plev[0], plev[-1] = pmin, pmax
        ilev = where(P <= plev)[0][0]

        # scan upward
        i, i0, temp, ptol0 = ilev, ilev, T, P
        sample_p, sample_t, sample_x, cloud = [], [], [[] for j in range(len(sp))], []
        pt = self.PotentialT(T, P, sp[sat])
        while (i >= 0 and sum(logical_not(sat)) > 0):
            temp = newton(lambda v: self.PotentialT(v, plev[i], sp[sat]) - pt, temp)
            newsat = IsSat(temp, plev[i], mixr, sp)
            if all(equal(sat, newsat)): 
                i -= 1
                continue
            j = where(sat != newsat)[0][0]
            sol = root(lambda v: self._eq1_(v, sp, sat, pt, mixr, j), (temp, plev[i], plev[i]))
            temp, ptol, pdry = tuple(sol.x)
            #print "Find Cloud Layer: T = %8.2f, P = %8.2f, Pdry = %8.2f, Species = " % (temp, ptol, pdry), sp[j]
            xcloud = array([mixr[ii] if not sat[ii] else SatVaporP(temp, sp[ii]) / pdry for ii in range(len(sp))])
            cloud = [CloudLayer(temp, ptol, xcloud, sp[j])] + cloud
            _sample_p, _sample_t = self.Adiabats(pt, ptol, ptol0, species = sp[sat], num = 3 * (i0 - i), endpoint = False)
            sample_p = hstack((_sample_p, sample_p))
            sample_t = hstack((_sample_t, sample_t))
            pdry = (_sample_p - sum(SatVaporP(_sample_t, sp[sat]), axis = -1)) / (1. + sum(mixr[logical_not(sat)]))
            _sample_x = zeros((len(sp), len(pdry)))
            for k in range(len(sp)):
                if sat[k]: _sample_x[k, :] = SatVaporP(_sample_t, sp[k]) / pdry
                else: _sample_x[k, :] = array([mixr[k]] * len(pdry))
            sample_x = hstack((_sample_x, sample_x))
            sat[j] = True;
            pt = self.PotentialT(temp, ptol, sp[sat])
            i0, ptol0 = i, ptol
        if ptol0 > pmin:
            _sample_p, _sample_t = self.Adiabats(pt, pmin, ptol0, species = sp, num = 3 * i0, endpoint = False)
            sample_p = hstack((_sample_p, sample_p))
            sample_t = hstack((_sample_t, sample_t))
            pdry = (_sample_p - sum(SatVaporP(_sample_t, sp[sat]), axis = -1)) / (1. + sum(mixr[logical_not(sat)]))
            _sample_x = zeros((len(sp), len(pdry)))
            for k in range(len(sp)):
                if sat[k]: _sample_x[k, :] = SatVaporP(_sample_t, sp[k]) / pdry
                else: _sample_x[k, :] = array([mixr[k]] * len(pdry))
            sample_x = hstack((_sample_x, sample_x))

        # scan downward
        i, i0, temp, ptol0 = ilev, ilev, T, P
        sat = sat0.copy()
        id = len(sample_p) - 1
        pt = self.PotentialT(T, P, sp[sat])
        while (i < num and sum(sat) > 0):
            temp = newton(lambda v: self.PotentialT(v, plev[i], sp[sat]) - pt, temp)
            newsat = IsSat(temp, plev[i], mixr, sp)
            if all(equal(sat, newsat)): 
                i += 1
                continue
            j = where(sat != newsat)[0][0]
            sol = root(lambda v: self._eq1_(v, sp, sat, pt, mixr, j), (temp, plev[i], plev[i]), tol = 1.E-6)
            temp, ptol, pdry = tuple(sol.x)
            #print "Find Cloud Layer: T = %8.2f, P = %8.2f, Pdry = %8.2f, Species = " % (temp, ptol, pdry), sp[j]
            xcloud = array([mixr[ii] if not sat[ii] else SatVaporP(temp, sp[ii]) / pdry for ii in range(len(sp))])
            cloud = cloud + [CloudLayer(temp, ptol, xcloud, sp[j])]
            _sample_p, _sample_t = self.Adiabats(pt, ptol0, ptol, species = sp[sat], num = 3 * (i - i0), endpoint = False)
            sample_p = hstack((sample_p, _sample_p))
            sample_t = hstack((sample_t, _sample_t))
            pdry = (_sample_p - sum(SatVaporP(_sample_t, sp[sat]), axis = -1)) / (1. + sum(mixr[logical_not(sat)]))
            _sample_x = zeros((len(sp), len(pdry)))
            for k in range(len(sp)):
                if sat[k]: _sample_x[k, :] = SatVaporP(_sample_t, sp[k]) / pdry
                else: _sample_x[k, :] = array([mixr[k]] * len(pdry))
            sample_x = hstack((sample_x, _sample_x))
            sat[j] = False;
            pt = self.PotentialT(temp, ptol, sp[sat])
            i0, ptol0 = i, ptol
        if ptol0 < pmax:
            _sample_p, _sample_t = self.Adiabats(pt, ptol0, pmax, species = sp[sat], num = 3 * (num - i0))
            sample_p = hstack((sample_p, _sample_p))
            sample_t = hstack((sample_t, _sample_t))
            pdry = (_sample_p - sum(SatVaporP(_sample_t, sp[sat]), axis = -1)) / (1. + sum(mixr[logical_not(sat)]))
            _sample_x = zeros((len(sp), len(_sample_p)))
            for k in range(len(sp)):
                if sat[k]: _sample_x[k, :] = SatVaporP(_sample_t, sp[k]) / pdry
                else: _sample_x[k, :] = array([mixr[k]] * len(pdry))
            sample_x = hstack((sample_x, _sample_x))

        # interpolation function and cell averaging
        result = CondensationSolution()
        result.cloud = cloud
        result.species = sp
        result.funcT = interp1d(sample_p, sample_t, kind = 'linear')
        result.funcX = interp1d(sample_p, sample_x, kind = 'linear')
        def cellT(pres):
            presw = zeros(len(pres) + 1)
            presw[0] = pres[0] * pres[0] / sqrt(pres[0] * pres[1])
            presw[-1] = pres[-1] * pres[-1] / sqrt(pres[-1] * pres[-2])
            for i in range(1, len(pres)):
                presw[i] = sqrt(pres[i] * pres[i - 1])
            return array([quad(result.funcT, presw[i], presw[i + 1])[0] / (presw[i + 1] - presw[i]) 
                for i in range(len(pres))])
        def cellX(pres):
            presw = zeros(len(pres) + 1)
            presw[0] = pres[0] * pres[0] / sqrt(pres[0] * pres[1])
            presw[-1] = pres[-1] * pres[-1] / sqrt(pres[-1] * pres[-2])
            for i in range(1, len(pres)):
                presw[i] = sqrt(pres[i] * pres[i - 1])
            cellx = zeros((len(sp), len(pres)))
            for j in range(len(sp)):
                funcx = interp1d(sample_p, sample_x[j, :], kind = 'linear')
                for i in range(len(pres)):
                    cellx[j, i] = quad(funcx, presw[i], presw[i + 1])[0] / (presw[i + 1] - presw[i]) 
            return cellx
        result.sample_p = sample_p
        result.sample_t = sample_t
        result.sample_x = sample_x
        result.cellT = cellT
        result.cellX = cellX
        result.make()
        return result
    def _eq1_(self, var, sp, sat, pt, mixr, j):
        tc, pc, pd = var
        return (self.PotentialT(tc, pc, sp[sat]) - pt,
                SatVaporP(tc, sp[j]) - mixr[j] * pd,
                pd * (1. + sum(mixr[logical_not(sat)])) + sum(SatVaporP(tc, sp[sat])) - pc
                )
    def _eq2_(self, mixr, sp, vt, ptol):
        tc = vt * (1. + sum(self.eps(sp) * mixr)) / (1. + sum(mixr))
        return SatVaporP(tc, sp) - mixr * ptol / (1. + sum(mixr))
    def __repr__(self):
        st = '%-30s: %12.3E kg/mol\n' % ('Molecular weight of dry air', self.mu)
        st += '%-30s: %12s m/(s^2)\n' % ('Gravity', self.grav)
        st += '%-30s: %12.3f J/(kg K)\n' % ('Specific heat capacity', self.cpdry)
        st += '%-30s: %12.3E pa\n' % ('Reference pressure', self.Pref)
        st += '%-30s: %12.3f K\n' % ('Reference temperature', self.Tref)
        st += '%-30s: %12.3f m' % ('Density scale height', self.Hscale)
        return st

class Saturn(Atmosphere):
    name = 'Saturn'
    xHe = 0.11
    H2, He = Hydrogen(), Helium()
    mu = (H2.mu + He.mu * xHe) / (1. + xHe)
    # scale cp by 0.8 to match voyer observation
    cp = 0.8 * (H2.cpvap * H2.mu + H2.cpvap * He.mu * xHe) / (H2.mu + He.mu * xHe)
    def __init__(self): Atmosphere.__init__(self, mu = self.mu, cp = self.cp)


