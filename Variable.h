#ifndef VARIABLE
#define VARIABLE
#include "Pooma/Arrays.h"
#include "Boundary.h"
#include "setups.h"
#include "NumericalScheme.h"

#define _maxorder_ 5
enum{ SnapshotView, MassView, AverageView, MassAverageView };

class VariableBase{
public:
    VariableBase(){}
    virtual ~VariableBase(){}
    virtual void printInfo(int level = 0) = 0;
    virtual void fixBoundary() = 0;
    virtual const char* getName(int s = 0) = 0;
    virtual void updateTendency() = 0;
    virtual int getSize() = 0;
    virtual void wallConstruct() = 0;
    virtual void updateView() = 0;
    virtual void view2cell() = 0;
    virtual void archive(int) = 0;
    virtual void add_archive(std::vector<int>, std::vector<double>) = 0;
    virtual Array<2, _AtomicType_>& getView(int s = 0) = 0;
    static long step_count;
    static Array<2, _AtomicType_> *mass;
    _Constructor_<_SpatialOrder_> constructor;
};
long VariableBase::step_count = 0;
Array<2, _AtomicType_> *VariableBase::mass = 0;

template<class T = double, class B = Dependent, int ViewType = SnapshotView>
class Variable : public VariableBase, public B
{
    typedef T Element_t;
    typedef B Boundary_t;
public:
    Variable(){}
    Variable(std::string _name) : name(_name) {
        cij = setups::cij;
        sicj = setups::sicj;
        cisj = setups::cisj;
        cxy = setups::cxy;

        cell.initialize(cxy); cell = 0;
        for (int i = 0; i < _TemporalOrder_; i++){
            cell_t[i].initialize(cxy); cell_t[i] = 0;
        }
        tendency.initialize(cij); tendency = 0;
        view.initialize(cij); view = 0;
        wallx.initialize(sicj); wallx = 0;
        wally.initialize(cisj); wally = 0;

        view = setups::ncvar[name];
    }

    void printInfo(int level = 0){
        std::cout << "Variable: " << name << std::endl;
        std::cout << "Domain: " << cij << " , " << cxy << std::endl;
        std::cout << "Boundary: " << B() << std::endl;
        if (level > 0){
            std::cout << "Maximum Value: " << max(cell) << ", Minimum Value: " << min(cell) << std::endl;
        }
        Interval<2> pij; pij = cij;
        while (pij[0].length() > 10 && pij[1].length() > 10){
            pij[0] = Interval<1>(0.75 * pij[0].first() + 0.25 * pij[0].last(),
                    0.25 * pij[0].first() + 0.75 * pij[0].last());
            pij[1] = Interval<1>(0.75 * pij[1].first() + 0.25 * pij[1].last(),
                    0.25 * pij[1].first() + 0.75 * pij[1].last());
        }
        if (level > 1){
            std::cout << "Sample Cell Value: " << pij << std::endl;
            std::cout << cell(pij);
        }
        if (level > 2){
            std::cout << "Sample Wall Value: " << pij << std::endl;
            std::cout << wallx(pij);
        }
        std::cout << std::endl;
    }

    void fixBoundary(){ this->fix(cell, cij); }

    void updateTendency(){ cell(cij) += tendency * setups::dt; tendency = 0; };

    void wallConstruct(){ constructor(cell, wallx, wally, cij); }
    
    void updateView(){
        switch (ViewType) {
            case SnapshotView:
                if (step_count % setups::frame == 0) view = cell(cij); 
                break;
            case MassView:
                if (step_count % setups::frame == 0) view = cell(cij) / mass->read(cij); 
                break;
            case AverageView:
                view += cell(cij);
                if (step_count % setups::frame == 0) view = view / setups::frame;
                break;
            case MassAverageView:
                view += cell(cij) / mass->read(cij); 
                if (step_count % setups::frame == 0) view = view / setups::frame;
                break;
        }
    }

    void view2cell(){
        if (ViewType == MassView || ViewType == MassAverageView) cell(cij) = view * mass->read(cij);
        else cell(cij) = view;
    }

    const char* getName(int) { return name.c_str(); }

    Array<2, T>& getView(int){ return view; }

    int getSize(){ return 1; }

    void archive(int t){ cell_t[t] = cell; }

    void add_archive(std::vector<int> index, std::vector<double> weight){
        if (index[0] == -1) cell = cell * weight[0];
        else cell = cell_t[index[0]] * weight[0];
        for (int i = 1; i < index.size(); i++) cell += cell_t[index[i]] * weight[i];
    }

    std::string name;
    Interval<2> cxy, cij, sicj, cisj;
    Array<2, Element_t> cell_t[_maxorder_], cell, tendency;
    Array<2, Vector<2, Element_t> > wallx, wally;
    Array<2, Element_t> view;
};

// default template arguments may not be used in partial specializations
template<int S, class T, class B, int ViewType>
class Variable<Vector<S, T>, B, ViewType> : public VariableBase, public B
{
    typedef Vector<S, T> Element_t;
    typedef B Boundary_t;
public:
    Variable(){for (int i = 0; i < S; i++) name[i] = new char[1];}
    Variable(std::initializer_list<std::string> _name){
        int j = 0;
        for (auto i = _name.begin(); i != _name.end(); i++, j++) { name[j] = *i; }
        cij = setups::cij;
        sicj = setups::sicj;
        cisj = setups::cisj;
        cxy = setups::cxy;

        cell.initialize(cxy); cell = 0;
        for (int i = 0; i < _TemporalOrder_; i++){
            cell_t[i].initialize(cxy); cell_t[i] = 0;
        }
        tendency.initialize(cij); tendency = 0;
        for (int i = 0; i < S; i++) {
            view[i].initialize(cij); view[i] = 0;
        }
        wallx.initialize(sicj); wallx = 0;
        wally.initialize(cisj); wally = 0;

        for (int i = 0; i < S; i++) view[i] = setups::ncvar[name[i]];
    }

    void printInfo(int level = 0){
        std::cout << "Variable: (";
        for (int i = 0; i < S - 1; i++) std::cout << name[i] << ", ";
        std::cout << name[S -1] << ")" << std::endl;
        std::cout << "Domain: " << cij << " , " << cxy << std::endl;
        std::cout << "Boundary: " << B() << std::endl;
        if (level > 0){
            Element_t vmax, vmin;
            for (int i = 0; i < S; i++){
                vmax(i) = max(cell.comp(i));
                vmin(i) = min(cell.comp(i));
            }
            std::cout << "Maximum Value: " << vmax << ", Minimum Value: " << vmin << std::endl;
        }
        Interval<2> pij; pij = cij;
        while (pij[0].length() > 10 && pij[1].length() > 10){
            pij[0] = Interval<1>(0.75 * pij[0].first() + 0.25 * pij[0].last(),
                    0.25 * pij[0].first() + 0.75 * pij[0].last());
            pij[1] = Interval<1>(0.75 * pij[1].first() + 0.25 * pij[1].last(),
                    0.25 * pij[1].first() + 0.75 * pij[1].last());
        }
        if (level > 1){
            std::cout << "Sample Cell Value: " << pij << std::endl;
            std::cout << cell(pij);
            std::cout << "Sample WallY Value: " << pij << std::endl;
            std::cout << wally(pij);
        }
        if (level > 2){
            std::cout << "Sample WallX Value: " << pij << std::endl;
            std::cout << wallx(pij);
            std::cout << "Sample WallY Value: " << pij << std::endl;
            std::cout << wally(pij);
        }
        std::cout << std::endl << std::endl;
    }

    void fixBoundary(){ this->fix(cell, cij); }

    void updateTendency(){ cell(cij) += tendency * setups::dt; tendency = 0;};

    void wallConstruct(){ constructor(cell, wallx, wally, cij); }

    void updateView(){
        switch (ViewType){
            case SnapshotView:
                if (step_count % setups::frame == 0) 
                    for (int i = 0; i < S; i++) view[i] = cell(cij).comp(i);
                break;
            case MassView:
                if (step_count % setups::frame == 0) 
                    for (int i = 0; i < S; i++) view[i] = cell(cij).comp(i) / mass->read(cij); 
                break;
            case AverageView:
                for (int i = 0; i < S; i++) view[i] += cell(cij).comp(i);
                if (step_count % setups::frame == 0)
                    for (int i = 0; i < S; i++) view[i] = view[i] / setups::frame;
                break;
            case MassAverageView:
                for (int i = 0; i < S; i++) view[i] += cell(cij).comp(i) / mass->read(cij); 
                if (step_count % setups::frame == 0) 
                    for (int i = 0; i < S; i++) view[i] = view[i] / setups::frame;
                break;
        }
        return;
    }

    void view2cell(){
        if (ViewType == MassView || ViewType == MassAverageView) 
            for (int i = 0; i < S; i++) cell(cij).comp(i) = view[i] * mass->read(cij);
        else 
            for (int i = 0; i < S; i++) cell(cij).comp(i) = view[i];
    }

    const char* getName(int s) { return name[s].c_str(); }

    Array<2, T>& getView(int s){ return view[s]; }

    int getSize(){ return S; };

    void archive(int t){ cell_t[t] = cell; }

    void add_archive(std::vector<int> index, std::vector<double> weight){
        if (index[0] == -1) cell = cell * weight[0];
        else cell = cell_t[index[0]] * weight[0];
        for (int i = 1; i < index.size(); i++) cell += cell_t[index[i]] * weight[i];
    }

    std::string name[S];
    Interval<2> cxy, cij, sicj, cisj;
    Array<2, Element_t> cell_t[_maxorder_], cell, tendency;
    Array<2, Vector<2, Element_t> > wallx, wally;
    Array<2, T> view[S];
};

class VariableList{
public:
    std::vector<VariableBase*> prog_list, var_list;
    void printInfo(int level = 0){
        std::cout << "******  variable infomation *****" << std::endl;
        for (int i = 0; i < var_list.size(); i++) var_list[i]->printInfo(level);
        std::cout << "*********************************" << std::endl;
    }
    void fixBoundary(){
        for (int i = 0; i < var_list.size(); i++) var_list[i]->fixBoundary();
    }
    void updateTendency(){
        for (int i = 0; i < prog_list.size(); i++) prog_list[i]->updateTendency();
    }
    void wallConstruct(){
        for (int i = 0; i < prog_list.size(); i++) prog_list[i]->wallConstruct();
    }
    void updateView(){
        for (int i = 0; i < var_list.size(); i++) var_list[i]->updateView();
    }
    void view2cell(){
        for (int i = 0; i < var_list.size(); i++) var_list[i]->view2cell();
    }
    void archive(int t){
        for (int i = 0; i < prog_list.size(); i++) prog_list[i]->archive(t);
    }
    void add_archive(std::vector<int> index, std::vector<double> weight){
        for (int i = 0; i < prog_list.size(); i++) prog_list[i]->add_archive(index, weight);
    }
    void ncwrite(double time){
        setups::current++;
        NcFile dataFile(setups::ncfile.c_str(), NcFile::Write);
        for (int i = 0; i < var_list.size(); i++){
            for (int s = 0; s < var_list[i]->getSize(); s++){
                dataFile.get_var(var_list[i]->getName(s))->put_rec(&var_list[i]->getView(s)(0, 0), setups::current);
                var_list[i]->getView(s) = 0;
            }
        }
        dataFile.get_var("time")->put_rec(&time, setups::current);
    }
};
#undef _maxorder_

#endif

