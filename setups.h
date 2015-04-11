#ifndef SETUPS
#define SETUPS
#include "Pooma/Arrays.h"
#include "netcdf.hh"
#include "utils.h"
#include <map>
#include "configure.h"

namespace setups{
    int nx, ny, glayer, frame;
    long current;
    double xlen, ylen, dx, dy, start, end, dt;
    std::map<std::string, Array<2, double> > ncvar;
    std::map<std::string, double> ncattr;
    Interval<2> cij, sij, sicj, cisj, cxy;
    std::string ncfile;

    void initialize(char *fname, double time_start, double time_end, double time_step, 
            double time_frame = 1, int restart = 0){
        ncfile = (std::string)fname;
        NcFile dataFile(fname, NcFile::ReadOnly);
        if (!dataFile.is_valid()){ ASSERT_FILE_NOT_FOUND(fname); }

        // read time varialbe
        Array<2, double> buffer;
        NcVar *data = dataFile.get_var("time");
        long *edges = data->edges();
        buffer.initialize(1, edges[0]);
        data->get(&buffer(0, 0), edges[0]);
        ncvar["time"].initialize(buffer.domain());
        ncvar["time"] = buffer;

        // set current time
        switch (restart) {
            long total;
            case 0:
                current = 0;
                start   = time_start;
                break;
            case -1:
                current = dataFile.get_dim("time")->size() - 1;
                start   = ncvar["time"](0, current);
                break;
            default: 
                total = dataFile.get_dim("time")->size();
                while (current++ < total && ncvar["time"](0, current) < time_start);
                current--;
                start   = ncvar["time"](0, current);
        }

        // read other variables
        for (int i = 0; i < dataFile.num_vars(); i++){
            data = dataFile.get_var(i);
            edges = data->edges();
            data->set_cur(current, 0, 0);
            switch (data->num_dims()){
                case 1:
                    buffer.initialize(1, edges[0]);
                    data->get(&buffer(0, 0), edges[0]);
                    break;
                case 2:
                    buffer.initialize(edges[1], edges[0]);
                    data->get(&buffer(0, 0), edges[0], edges[1]);
                    break;
                case 3:
                    buffer.initialize(edges[2], edges[1]);
                    data->get(&buffer(0, 0), 1, edges[1], edges[2]);
                    break;
            }
            ncvar[data->name()].initialize(buffer.domain());
            ncvar[data->name()] = buffer;
        }

        end     = time_end;
        dt      = time_step;
        frame   = time_frame;
        nx      = dataFile.get_att("num_grids_in_x")->as_int(0);
        ny      = dataFile.get_att("num_grids_in_y")->as_int(0);
        xlen    = dataFile.get_att("length_x")->as_double(0);
        ylen    = dataFile.get_att("length_y")->as_double(0);
        dx      = dataFile.get_att("grid_size_x")->as_double(0);
        dy      = dataFile.get_att("grid_size_y")->as_double(0);

        cij     = Interval<2>(Interval<1>(0, nx - 1), Interval<1>(0, ny - 1));
        sij[0]  = Interval<1>(cij[0].first() - 1, cij[0].last() + 1);
        sij[1]  = Interval<1>(cij[1].first() - 1, cij[1].last() + 1);
        sicj    = Interval<2>(sij[0], cij[1]);
        cisj    = Interval<2>(cij[0], sij[1]);
        glayer  = _SpatialOrder_;
        cxy     = Interval<2>(Interval<1>(- glayer, nx - 1 + glayer), Interval<1>(- glayer, ny - 1 + glayer));

        // read other global attributes
        for (int i = 0; i < dataFile.num_atts(); i++)
            ncattr[dataFile.get_att(i)->name()] = dataFile.get_att(i)->as_double(0);
    }
}

#endif
