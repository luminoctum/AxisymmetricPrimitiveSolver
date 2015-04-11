#include <iostream>
#include "getopt.h"
#include "stdio.h"
#include "math.h"
#include "NumericalScheme.h"
#include "Primitive.h"
#include "jansson.h"
using namespace std;

static int newline_offset(const char *text){
    const char *newline = strchr(text, '\n');
    if (!newline) return strlen(text);
    else return (int)(newline - text);
}

int main(int argc, char **argv){
    FILE *fp; long fsize = 0; char *control;
    double time_start, time_end, time_step;
    int time_frame, restart;
    json_t *root; json_error_t error;
    char *input_name = "dynamics.nc";
    char *output_name = "dynamics.nc";
    char *control_name = "control.json";

    int op = 0;
    while ((op = getopt(argc, argv, "i:o:c:")) != -1)
        switch (op){
            case 'i':
                input_name = optarg;
                break;
            case 'o':
                output_name = optarg;
                break;
        }
    printf("Model runs in %s\n", input_name);
    Pooma::initialize(argc, argv);

    if (!(fp = fopen("control.json", "r"))) perror("error opening control file\n");
    while ( fgetc(fp) != EOF ) fsize++;
    rewind(fp);
    control = (char*) malloc(fsize * sizeof(char));
    fread(control, sizeof(char), fsize, fp);
    fclose(fp);

    root = json_loads(control, 0, &error);
    if (!root) {
        fprintf(stderr, "error: on line %d: %s\n", error.line, error.text);
        return 1;
    }
    delete control;
    if (!json_is_array(root)){
        fprintf(stderr, "error: root is not an array\n");
        json_decref(root);
        return 1;
    }
    for (int i = 0; i < json_array_size(root); i++){
        json_t *data, *j_time_start, *j_time_end, *j_time_step, *j_time_frame, *j_restart;
        data = json_array_get(root, i);
        if (!json_is_object(data)){
            fprintf(stderr, "error: field %d is not an valid field\n", i + 1);
            json_decref(root);
            return 1;
        }

        j_time_start = json_object_get(data, "time_start");
        if (!json_is_number(j_time_start)){
            fprintf(stderr, "error: time_start is not a number", i + 1);
            json_decref(root);
            return 1;
        }
        time_start = json_real_value(j_time_start);

        j_time_end = json_object_get(data, "time_end");
        if (!json_is_number(j_time_end)){
            fprintf(stderr, "error: time_end is not a number", i + 1);
            json_decref(root);
            return 1;
        }
        time_end = json_real_value(j_time_end);

        j_time_step = json_object_get(data, "time_step");
        if (!json_is_number(j_time_step)){
            fprintf(stderr, "error: time_step is not a number", i + 1);
            json_decref(root);
            return 1;
        }
        time_step = json_real_value(j_time_step);

        j_time_frame = json_object_get(data, "steps_per_frame");
        if (!json_is_integer(j_time_frame)){
            fprintf(stderr, "error: steps_per_frame is not a integer", i + 1);
            json_decref(root);
            return 1;
        }
        time_frame = json_integer_value(j_time_frame);

        j_restart = json_object_get(data, "restart");
        if (!json_is_integer(j_restart)){
            fprintf(stderr, "error: restart is not a integer", i + 1);
            json_decref(root);
            return 1;
        }
        restart = json_integer_value(j_restart);
    }

    setups::initialize(input_name, time_start, time_end, time_step, time_frame, restart);
    _System_ system;
    _TimeMarching_<_TemporalOrder_> solve;
    solve(system, setups::start, setups::end, setups::dt, setups::frame);

    Pooma::finalize();
}
