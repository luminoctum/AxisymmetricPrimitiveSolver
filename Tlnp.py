#! /usr/bin/env python2.7
#import matplotlib
#matplotlib.use('Agg')
from scipy.optimize import fsolve
from pycli.data.SaturnData import *
from pycli.thermodynamics.Thermodynamics import *

def find_warm_t0(atm, cold_t0, sp):
    cold_sol = atm.EqCondensation(cold_t0, atm.Pref, array([1., 1.]), sp, atm.Pref, 30.E5, num = 200.)
    def _eq_(warm_t0):
        warm_t0 = squeeze(warm_t0)
        warm_sol = atm.EqCondensation(warm_t0, atm.Pref, array([1., 1.]), sp, atm.Pref, 30.E5, num = 50.)
        vt0 = atm.VirtualT(warm_sol.cloud[1].temp, warm_sol.cloud[1].mixr, sp)
        vt1 = atm.VirtualT(cold_sol.funcT(warm_sol.cloud[1].ptol), cold_sol.funcX(warm_sol.cloud[1].ptol), sp)
        return vt1 - vt0
    return fsolve(_eq_, cold_t0 + 10.)[0]

# T-lnp diagram
class Thermodiagram:
    def __init__(self, 
            tmin = 130, tmax = 360, nt = 200, 
            pmin = 1.E5, pmax = 30.E5, np = 200, 
            atm = None, sp = None):
        self.pmin, self.pmax, self.np = pmin, pmax, np
        self.tmin, self.tmax, self.nt = tmin, tmax, nt
        self.ptol = logspace(log10(pmax), log10(pmin), np)
        self.temp = linspace(tmin, tmax, nt)
        pticks = array(range(1, 9) + range(10, 50, 5))
        self.pticks = (hstack([arange(0.1, 1, 0.2), pticks]))
        self.tticks = arange(100, 400, 20)
        self.atm = atm
        self.sp = sp
    def draw(self):
        side_plot = False
        figure(figsize = (14, 10))
        if side_plot:
            axmain = axes([0.1, 0.1, 0.63, 0.8])
            axside = axes([0.1+0.63+0.02, 0.1, 0.15, 0.8])
            axside.minorticks_on()
            axside.tick_params(width = 2, length = 8)
            axside.tick_params(width = 2, which = 'minor', length = 4)
            axside.set_yscale('log')
            axside.set_yticks(self.pticks * 1E5)
            axside.yaxis.set_major_formatter(NullFormatter())
            axside.set_ylim([self.pmax, self.pmin])
        else:
            axmain = axes()
        axmain.minorticks_on()
        axmain.tick_params(width = 2, length = 8)
        axmain.tick_params(width = 2, which = 'minor', length = 4)
        axmain.set_yscale('log')
        axmain.set_yticks(self.pticks * 1E5)
        axmain.set_yticklabels(self.pticks)
        axmain.set_xticks(self.tticks)
        axmain.set_ylabel('Pressure (bar)')
        axmain.set_xlabel('Temperature (K)')
        axmain.set_xlim([self.tmin, self.tmax])
        axmain.set_ylim([self.pmax, self.pmin])

        # main plot
        cold_t0 = self.atm.Tref
        print 'Finding warm moist adiabats...'
        warm_t0 = find_warm_t0(self.atm, cold_t0, self.sp)
        print 'warm  = ', warm_t0
        #warm_t0 = 139.690562

        print 'Calculating moist adiabats...'
        warm_sol = self.atm.EqCondensation(warm_t0, self.atm.Pref, array([1,1]), self.sp, 0.9 * self.pmin, 1.1 * self.pmax)
        warm_t = warm_sol.cellT(self.ptol[::-1])[::-1]
        cold_sol = self.atm.EqCondensation(cold_t0, self.atm.Pref, array([1,1]), self.sp, 0.9 * self.pmin, 1.1 * self.pmax)
        cold_t = cold_sol.cellT(self.ptol[::-1])[::-1]
        print 'Voyager observation...'
        data = SaturnData()
        data.import_rss()
        axmain.plot(data.rss['lindal85'].temp, data.rss['lindal85'].pres, 'k', linewidth = 4)
        print 'Cloud base...'
        axmain.plot([self.tmin, self.tmax], [warm_sol.cloud[0].ptol, warm_sol.cloud[0].ptol], 'g--', linewidth = 3)
        axmain.plot([self.tmin, self.tmax], [warm_sol.cloud[1].ptol, warm_sol.cloud[1].ptol], 'g-', linewidth = 3)
        print 'Moist aidabats...'
        axmain.plot(warm_t, self.ptol, 'r-')
        i = where(self.ptol < warm_sol.cloud[1].ptol)[0][0]
        axmain.plot(cold_t[i:], self.ptol[i:], 'r--')
        print 'Constant mixing ratio...'
        mixr = zeros((len(self.sp), self.nt, self.np))
        vt = zeros((self.nt, self.np))
        for i in range(self.nt):
            for j in range(self.np):
                mixr[:, i, j] = SatMixingR(self.temp[i], self.ptol[j], [self.sp[0], Water()])
                vt[i, j] = self.atm.VirtualT(self.temp[i], mixr[:, i, j], [self.sp[0], Water()])
        h = axmain.contour(self.temp, self.ptol, 100. * mixr[1, :, :].T,
                array([0.1, 0.2, 0.4, 0.7, 1.0, 100. * self.sp[1].mmr, 2.]),
                colors = 'm', linestyles = '--', linewidths = 1)
        h.collections[5].set_linestyle('-')
        h.collections[5].set_linewidth('2')
        axmain.clabel(h, inline = 1, fontsize = 15, fmt = '%3.1f')
        print 'Virtual temperature...'
        vt0 = self.atm.VirtualT(warm_sol.cloud[1].temp, warm_sol.cloud[1].mixr, self.sp)
        h = axmain.contour(self.temp, self.ptol, vt.T, 
                arange(vt0 - 40., vt0 + 20, 10),
                colors = 'c', linestyles = '--', linewidths = 2)
        axmain.clabel(h, inline = 1, fontsize = 15, fmt = '%3.0f')
        if not side_plot: 
            savefig('Tlnp-%3.1f.pdf' % (100.*self.sp[1].mmr,), bbox_inches = 'tight')
            print 'Done!'
            return
        # side plot
        #axside.set_xticks([130., 140., 150])
        #axside.set_xlim([125., 150.])
        #cold_t0, warm_t0 = diag.atmos.Tref, 140.85
        #axside.plot(warm_ept, diag.ptol, 'r-')
        print 'Done!'

if __name__ == '__main__':
    H2O = Water(mmr = 1.1E-2)
    NH3 = Ammonia(mmr = 4.E-4)
    diagram = Thermodiagram(
            tmin = 240, tmax = 360, nt = 200,
            #tmin = 100, tmax = 360, nt = 200,
            pmin = 0.7E5, pmax = 30.E5, np = 200, 
            atm = Saturn(), sp = [NH3, H2O])
    diagram.draw()
