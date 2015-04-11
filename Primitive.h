#ifndef PRIMITIVE
#define PRIMITIVE
#include <vector>
#include "Variable.h"
#include "Stencil.h"
#include "MicroPhysics.h"

/**
 * This class is the base class for the class that actually does the computation.
 * It defines many useful common features such as :
 * stencil operation : cdx, cdy, fdx, fdy, finty, binty, esum, eavg, avg, laplace
 * constant test arrays : ones, onesx, onesy
 * domain range arrays : cij, sicj, cisj
 * mass variable pointer : mass_pt
 * step counter : step_count
 */
class ModelBase{
protected:
    /// runtime information
    bool unset_info;
    std::vector<std::string> info_header_fmt;
    std::vector<std::string> info_header;
    std::vector<std::string> info_value_fmt;
    std::vector<double> info_value;
    void printInfoHeader(){
        for (int i = 0; i < info_header.size(); i++)
            printf(info_header_fmt[i].c_str(), info_header[i].c_str());
        std::cout << std::endl;
    }
    void printInfoValue(){
        for (int i = 0; i < info_header.size(); i++)
            printf(info_value_fmt[i].c_str(), info_value[i]);
        std::cout << std::endl;
    }

    /// stencil operations
    Stencil<ElementSum> esum;
    Stencil<ElementAverage> eavg;
    Stencil<Average> avg;
    Stencil<FiniteDifference<D1, O2, DimX> > cdx;
    Stencil<FiniteDifference<D1, O2, DimY> > cdy;
    Stencil<FiniteDifference<D1, O1, DimX> > fdx;
    Stencil<FiniteDifference<D1, O1, DimY> > fdy;
    Integrate<Forward, Trapz, DimY> finty;
    Integrate<Backward, Staggered, DimY> binty;
    Stencil<Laplace> laplace;

    // help variables
    Array<2, double, ConstantFunction> ones;
    Array<2, Vector<2, double>, ConstantFunction> onesx, onesy;
    Interval<2> cij, sicj, cisj, sij;
    Array<2, _AtomicType_>*& mass_pt;

    // initializer and destructor
    ModelBase() : step_count(VariableBase::step_count), mass_pt(VariableBase::mass)
    {
        cpu_tic = clock();
        cpu_time = 0;
        unset_info = true;
        cij     = setups::cij;
        sicj    = setups::sicj;
        cisj    = setups::cisj;
        sij     = setups::sij;
        FiniteDifferenceBase::dx = setups::dx;
        FiniteDifferenceBase::dy = setups::dy;
        ones.initialize(cij); ones.engine().setConstant(1.);
        onesx.initialize(sicj); onesx.engine().setConstant(Vector<2>(1., 1.));
        onesy.initialize(cisj); onesy.engine().setConstant(Vector<2>(1., 1.));
    }
    ~ModelBase(){
        int hour = floor(cpu_time / 3600.);
        int min = floor((cpu_time - 3600. * hour) / 60);
        int sec = floor((cpu_time - 3600. * hour - 60 * min));
        printf("Total time elapsed: %d h %d m %d s\n", hour, min, sec);
        std::cout << "********** Simulation completed **********" << std::endl;
    }
public:
    /// store the total steps, this is a reference to the static varialbe in 'class Variable'
    long& step_count;
    /// store the model time and cpu time
    double model_time, cpu_time;
    clock_t cpu_tic, cpu_toc;
};


/** This class is actually solving the primitive equation in log pressure coordinate.
 * It is derived from 'class ModelBase' and 'class VariableList'.
 * 'class ModelBase' does the common initialization processes and
 * 'class VariableLise does the common operations of variables such as advection, write into ncfile, etc.
 */
class Primitive : public ModelBase, public VariableList{
    typedef _AtomicType_ T;
protected:
    /// these variables calculate the actual flux from two wall flux (Godunov method)
    Stencil<_FluxCalculator_<T, DimX> > fux;
    Stencil<_FluxCalculator_<T, DimY> > fuy;

    /// these variables are constants in the model, loaded from nc file
    double grav, T0, cp, f, mu; Vector<2, T> eps;
    Array<2, T> mass, rdist, t_ov_tc, tv0, ptol, massx, massy;

    /// microphysics variables
    Water H2O; Ammonia NH3;
    CondensateList<2> species;

    /// viscosity
    Array<2, T> viscosity;
public:
    /** 
     * definition of variables : Variable<type, boundary_condition, view_type>
     * where type could be : float, double, Vector<n, x>, etc
     * boundary_condition could be : Mirror, Fixed, ConstExtrap, etc
     * view_type could be : SnapshotView, MassView, AverageView, MassAverageView, etc
     */
    Variable<T, Boundary<Reflective, ConstExtrap, Absorbing<T>, Absorbing<T, 10, 20> >
    //Variable<T, Boundary<Reflective, ConstExtrap, Absorbing<T>, ConstExtrap>
        , MassView> uwind;
    Variable<T, Boundary<Reflective, ConstExtrap, Absorbing<T>, Absorbing<T, 10, 20> > 
    //Variable<T, Boundary<Reflective, ConstExtrap, Absorbing<T>, ConstExtrap> 
        > vwind;
    Variable<T, Boundary<Mirror, Fixed<T>, ConstExtrap, Fixed<T> >
        > theta;
    Variable<T, Boundary<Mirror, Fixed<T>, ConstExtrap, Fixed<T> >
        > paraH2;
    Variable<Vector<2, T>, Boundary<Mirror, Fixed<Vector<2, T> >, ConstExtrap, ConstExtrap>
        > mixr;
    Variable<Vector<2, T>, Dependent
        , AverageView> svp, rh, liq;
    Variable<T, Dependent, 
        MassView> wwind;
    Variable<T, Dependent
        > temp, tempv, pdry;
    Variable<T, Boundary<Reflective, Fixed<T>, ConstExtrap, LinearExtrap>
        > phi;

    /**
     * Initilization process.
     * 1) retrieve constants used in the model from nc file
     * 2) define wind variable for calculation of advection
     * 3) define mass variable for output
     * 4) define varialbe list for output
     * 5) define prognostic variable
     * 6) read variable from ncfile
     * 7) apply boundary condition
     */
    Primitive() : uwind("uwind"), vwind("vwind"), theta("tc"), paraH2("paraH2"), mixr({"xH2O", "xNH3"}), 
        wwind("wwind"), phi("phi"), temp("temp"), tempv("tempv"), pdry("pdry"),
        svp({"svpH2O", "svpNH3"}), rh({"hH2O", "hNH3"}), liq({"qH2O", "qNH3"}),
        H2O(1.3E-2), NH3(4E-4) ,species(H2O, NH3)
    {
        mu      = setups::ncattr["mean_molecular_weight"];
        grav    = setups::ncattr["gravity"];
        T0      = setups::ncattr["reference_temperature"];
        cp      = setups::ncattr["specific_heat_capacity"];
        f       = setups::ncattr["coriolis_parameter"];
        eps     = species.mu / mu;

        fux.function().xwind = &uwind.wallx;
        fuy.function().ywind = &wwind.wally;

        _ncInitialize7_(mass, rdist, t_ov_tc, tv0, ptol, massx, massy);

        mass_pt = &mass;
        var_list = {
            &uwind, &vwind, &wwind, &theta, &paraH2, &mixr, &phi, &temp, &tempv, &pdry, &svp, &rh, &liq
        };
        prog_list = {
            &uwind, &vwind, &theta, &paraH2, &mixr
        };
        view2cell();
        fixBoundary();

        viscosity.initialize(cij);
        viscosity = 0;
        double d = 20.;
        for (int i = 0; i < setups::nx; i++)
            for (int j = 0; j < setups::ny; j++)
                viscosity(i, j) =  _viscosity_ / d
                    * (_max_(d - i, 0) + _max_(i - setups::nx + d, 0) + _max_(d - j, 0) + _max_(j - setups::ny + d, 0));
        //std::cout << viscosity << std::endl;
        //exit(0);

        std::cout << "********** Simulation start **********" << std::endl;
    };
    /**
     * update diagnostic variables
     */
    void updateDiagnostics(){
        wwind.cell(cij) = eavg(wwind.wally, cij);
        temp.cell(cij) = t_ov_tc * theta.cell(cij);
        // stencil output is zero based, need to specify cij explicitly
        tempv.cell(cij) = temp.cell(cij) * (1. + esum(mixr.cell, cij))/(1. + esum(eps * mixr.cell, cij));
        finty(grav / T0 * (tempv.cell(cij) - tv0), phi.cell, cij);
    }
    /**
     * update body forces
     */
    void updateBodyforce(){
        uwind.tendency += - fdx(avg(phi.cell(sicj)) * massx, cij) + mass * phi.cell(cij) / rdist
            + mass * vwind.cell(cij) * (f + vwind.cell(cij) / rdist);
        vwind.tendency += - uwind.cell(cij) / mass * (f + vwind.cell(cij) / rdist);
        theta.tendency += theta.cell(cij) / (cp * temp.cell(cij)) * esum(species.Lv * eps * liq.cell(cij)) / setups::dt;
    }
    /**
     * update microphysics
     */
    void updateMicrophysics(){
        // saturation vapor pressure
        svp.cell.comp(0) = H2O.svp_from_t(temp.cell);
        svp.cell.comp(1) = NH3.svp_from_t(temp.cell);
        // relative humidity: h = x * pdry / svp
        rh.cell = mixr.cell * pdry.cell / svp.cell;
        // condensing liquid
        for (int s = 0; s < 2; s++){
            // in case at some point cfl is voilated, but instability will occur
            //rh.cell.comp(s) = where(rh.cell.comp(s) < 0., 0., rh.cell.comp(s));
            liq.cell.comp(s) = where(rh.cell.comp(s) > 1., (rh.cell.comp(s) - 1.) * svp.cell.comp(s) / pdry.cell, 0.);
            rh.cell.comp(s) = where(liq.cell.comp(s) > 0., 1., rh.cell.comp(s));
        }
        pdry.cell(cij) = ptol / (1. + esum(mixr.cell, cij));
        // mixing ratio: x = h * svp / pdry
        mixr.cell(cij) = rh.cell(cij) * svp.cell(cij) / pdry.cell(cij);
    }
    /**
     * deal with advective terms
     */
    void advection(){
        wallConstruct();
        /* uwind at the left boundary must be zero */
        uwind.wallx(0, AllDomain<1>()).comp(0) = 0;
        uwind.wallx(-1, AllDomain<1>()).comp(1) = 0;
    
        binty(- fdx(fux(onesx)), wwind.wally, cij);
        uwind.tendency -= fdx(fux(uwind.wallx) / massx) + fdy(fuy(uwind.wally) / massy);
        vwind.tendency -= (fdx(fux(vwind.wallx)) + fdy(fuy(vwind.wally))) / mass;
        theta.tendency -= (fdx(fux(theta.wallx)) + fdy(fuy(theta.wally))) / mass;
        paraH2.tendency -= (fdx(fux(paraH2.wallx)) + fdy(fuy(paraH2.wally))) / mass;
        mixr.tendency  -= (fdx(fux(mixr.wallx)) + fdy(fuy(mixr.wally))) / mass;
    }
    /**
     * viscous dissipation of wind velocity
     */
    inline void dissipation(){
        uwind.tendency += _viscosity_ / setups::dt * laplace(uwind.cell, cij);
        vwind.tendency += _viscosity_ / setups::dt * laplace(vwind.cell, cij);
        //uwind.tendency += viscosity / setups::dt * laplace(uwind.cell, cij);
        //vwind.tendency += viscosity / setups::dt * laplace(vwind.cell, cij);
        //uwind.tendency += _viscosity_ / setups::dt * laplace(laplace(uwind.cell, sij));
        //vwind.tendency += _viscosity_ / setups::dt * laplace(laplace(vwind.cell, sij));
    }
    /**
     * print out run time information and write variables into ncfile
     */
    inline void observe(){
        cpu_toc = clock();
        double _time = (double)(cpu_toc - cpu_tic) / CLOCKS_PER_SEC;
        if (unset_info){
            info_header_fmt = {"%-6s", "%-12s", "%-8s", "%-8s", "%-8s", "%-12s"};
            info_header     = {"No.", "time", "cpu", "uwind", "wwind", "xNH3"};
            info_value_fmt  = {"%-6.0f", "%-12.1f", "%-8.1f", "%-8.1f", "%-8.1f", "%-12.2E"}; 
            unset_info = false;
            printInfoHeader();
        }
        info_value = {
            (double)step_count,
            model_time,
            _time,
            max(abs(uwind.cell(cij)) / mass),
            max(abs(wwind.cell(cij)) / mass),
            max(mixr.cell(cij).comp(1))
        };
        printInfoValue();
        ncwrite(model_time);
        cpu_time += _time;
        cpu_tic = clock();
    }
    /**
     * time marching
     */
    inline void do_step(){
        advection();
        updateBodyforce();
        dissipation();
        updateMicrophysics();
        updateTendency();
        updateDiagnostics();
        fixBoundary();
    }
};

#endif
