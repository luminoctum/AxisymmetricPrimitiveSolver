#! /usr/bin/env python2.7
from InitBase import *
from collections import OrderedDict
from pycli.ode.gridop import *
from pycli.thermodynamics.Thermodynamics import *
from pycli.data.SaturnData import *
from Tlnp import find_warm_t0

class InitPrimitive(InitBase):
    def set_variables(self):
        cold_t0 = self.atm.Tref
        warm_t0 = find_warm_t0(self.atm, cold_t0, self.sp)
        print warm_t0
        #warm_t0 = 139.69

        pmin, pmax = min(self.ptol), max(self.ptol)
        cold_sol = self.atm.EqCondensation(cold_t0, self.atm.Pref, array([1,1]), self.sp, 0.9 * pmin, 1.1 * pmax)
        cold_t = cold_sol.funcT(self.ptol)
        warm_sol = self.atm.EqCondensation(warm_t0, self.atm.Pref, array([1,1]), self.sp, 0.9 * pmin, 1.1 * pmax)
        warm_t = warm_sol.funcT(self.ptol)
        # above 30.E3 pa, cold_t is replaced by voyager observation
        data = SaturnData()
        data.import_rss(pmin = 0.5E3, pmax = 50.E3) # actual data needed is from 1E3 to 30E3
        voyager = interp1d(data.rss['lindal85'].pres[::-1], data.rss['lindal85'].temp[::-1], kind = 'linear')
        i = where(self.ptol < 30.E3)[0][0]
        cold_t[i:] = voyager(self.ptol[i:])
        # below warm_t cloud base, cold_t is replaced by warm_t
        i = where(self.ptol < warm_sol.cloud[1].ptol)[0][0]
        cold_t[:i] = warm_t[:i]
        
        cold_tc = self.atm.PotentialT(cold_t, self.ptol, [])
        warm_tc = self.atm.PotentialT(warm_t, self.ptol, [])
        cold_tv = self.atm.SatVirtualT(cold_t, self.ptol, self.sp)
        warm_tv = self.atm.SatVirtualT(warm_t, self.ptol, self.sp)

        Ptol = array([self.ptol for x in self.xaxis])
        Rdist = array([self.xaxis for y in self.yaxis]).T
        Cold_t = array([cold_t for x in self.xaxis])
        Cold_tv = array([cold_tv for x in self.xaxis])
        T = array([cold_t + (warm_t - cold_t).clip(0, 1.E6) * exp(- x**2 / (2. * self.sigma**2)) for x in self.xaxis])
        #Tc = T * (self.atmos.Pref / Ptol)**(self.atmos.Rgas/self.atmos.cpvap)
        Tc = self.atm.PotentialT(T, Ptol, [])
        Tv = self.atm.SatVirtualT(T, Ptol, self.sp)
        Eta = SatMixingR(T, Ptol, self.sp)
        
        Pdry = Ptol / (1. + sum(Eta, axis = -1))
        svp = SatVaporP(T, self.sp)
        RelativeH = Eta * dstack((Pdry, Pdry)) / svp
        Phi = zeros((self.nx, self.ny))
        for i in range(self.ny):
            Phi[:, i] = trapz(self.atm.grav * (Tv[:, :i + 1] - Cold_tv[:, :i + 1]) / self.atm.Tref, self.yaxis[:i + 1], axis = 1)
        r0 = 1.E6
        Mass = array([x/r0 * exp(- self.yaxis / self.atm.Hscale) for x in self.xaxis])
        MassX = array([x/r0 * exp(- self.yaxis / self.atm.Hscale) for x in self.xaxisb])
        # eliminate sigularity
        MassX[0, :] = 1.E-6
        MassY = array([x/r0 * exp(- self.yaxisb / self.atm.Hscale) for x in self.xaxis])
        T_ov_tc = array([ exp(- self.atm.grav * self.yaxis / (self.atm.cpdry * self.atm.Tref)) for x in self.xaxis])
        # para-hydrogen fraction
        data = genfromtxt('/home/cli/lib/pycli/data/paraH2Fraction.txt', skip_header = 1)
        fpara = interp1d(data[:, 0], data[:, 1], 'linear')

        self.var['radial'] = self.xaxis
        self.var['radialh'] = self.xaxisb
        self.var['zplev'] = self.yaxis
        self.var['zplevh'] = self.yaxisb
        self.var['rdist'] = Rdist
        self.var['ptol'] = Ptol
        self.var['pdry'] = Pdry
        self.var['tc'] = Tc
        self.var['phi'] = Phi
        self.var['tv0'] = Cold_tv
        self.var['mass'] = Mass
        self.var['massx'] = MassX
        self.var['massy'] = MassY
        self.var['t_ov_tc'] = T_ov_tc
        self.var['xNH3'] = Eta[:, :, 0]
        self.var['xH2O'] = Eta[:, :, 1]
        self.var['temp'] = T
        self.var['tempv'] = Tv
        self.var['hNH3'] = RelativeH[:, :, 0]
        self.var['hH2O'] = RelativeH[:, :, 1]
        self.var['paraH2'] = fpara(T)

if __name__ == '__main__':
    H2O = Water(mmr = 1.1E-2)
    NH3 = Ammonia(mmr = 4.E-4)
    pmin, pmax = 1.E3, 30E5
    nlev = 128
    nx, ny = 2*nlev, nlev
    sigma = 2.E5
    saturn = Saturn()

    xlen = 10.E6
    xaxisb = linspace(0, xlen, nx + 1)
    xaxis = tohalf(xaxisb)
    paxis = logspace(log10(pmax), log10(pmin), ny)
    yaxis = saturn.Hscale * log(saturn.Pref / paxis)
    yaxisb = tohalf(yaxis, ext = 'both')
    ylen = yaxisb[-1] - yaxisb[0]

    control = OrderedDict([
            ('num_grids_in_x', nx),
            ('num_grids_in_y', ny),
            ('length_x', xlen),
            ('length_y', ylen),
            ('grid_size_x', xlen / nx),
            ('grid_size_y', ylen / ny),
            ('time_start', 0.),
            ('time_step', 10.),
            ('time_end', 180000.),
            ('steps_per_frame', 20),
            ('restart', 0),
    ])
    attrlist = OrderedDict([
            ('heating_width', sigma),
            ('water_mixing_ratio', H2O.mmr),
            ('ammonia_mixing_ratio', NH3.mmr),
            ('reference_temperature', saturn.Tref),
            ('coriolis_parameter', 2.128E-4),
            ('specific_heat_capacity', saturn.cpdry),
            ('mean_molecular_weight', saturn.mu),
            ('gravity', saturn.grav),
    ])
    varlist = [
            ['time', 'di', 'time', 's'],
            ['radial', 'dk', 'radial distance', 'm'],
            ['radialh', 'dK', 'staggered radial distance', 'm'],
            ['zplev', 'dj', 'distance in log pressure coordinate', 'm'],
            ['zplevh', 'dJ', 'staggered distance in log pressure coordinate', 'm'],
            ['rdist', 'jk', 'broadcasted radial distance', 'm'],
            ['ptol', 'jk', 'broadcasted pressure', 'pa'],
            ['pdry', 'ijk', 'dry air pressure', 'pa'],
            ['uwind', 'ijk', 'zonal wind', 'm/s'],
            ['vwind', 'ijk', 'meridional wind', 'm/s'],
            ['wwind', 'ijk', 'vertical wind', 'm / s'],
            ['tc', 'ijk', 'potential temperature', 'K'],
            ['phi', 'ijk', 'geopotential height anomaly', 'm^2 / s^2'],
            ['tv0', 'jk', 'basic virtual temperature', 'K'],
            ['temp', 'ijk', 'temperature', 'K'],
            ['tempv', 'ijk', 'virtual temperature', 'K'],
            ['mass', 'jk', 'coupled mass variable = r/r0 * exp(-z / H0)', ''],
            ['massx', 'jK', 'X staggered mass variable = r/r0 * exp(-z / H0)', ''],
            ['massy', 'Jk', 'Y staggered mass variable = r/r0 * exp(-z / H0)', ''],
            ['t_ov_tc', 'jk', 'T over theta = exp(- gz / cpT0)', ''],
            ['xH2O', 'ijk', 'H2O mixing ratio', 'mol/mol'],
            ['xNH3', 'ijk', 'NH3 mixing ratio', 'mol/mol'],
            ['svpH2O', 'ijk', 'H2O saturation vapor pressure', 'pa'],
            ['svpNH3', 'ijk', 'NH3 saturation vapor pressure', 'pa'],
            ['hH2O', 'ijk', 'H2O relative humidity', ''],
            ['hNH3', 'ijk', 'NH3 relative humidity', ''],
            ['qH2O', 'ijk', 'liquid H2O mixing ratio', 'g/g'],
            ['qNH3', 'ijk', 'liquid NH3 mixing ratio', 'g/g'],
            ['paraH2', 'ijk', 'para-H2 fraction', ''],
    ]

    model = InitPrimitive(control, attrlist, varlist, 
            output = 'dy256.s%3.1f.x%3.1f.nc' % (1E-5*sigma, 100.*H2O.mmr))
    model.xaxis = xaxis
    model.xaxisb = xaxisb
    model.yaxis = yaxis
    model.yaxisb = yaxisb
    model.ptol = paxis
    model.atm = saturn
    model.sp = [NH3, H2O]
    model.sigma = sigma
    model.initialize()
