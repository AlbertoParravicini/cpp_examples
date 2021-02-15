#pragma once

#include <getopt.h>
#include <cstdlib>

//////////////////////////////
//////////////////////////////

#define DEBUG true
#define NUM_ITER 1
#define DEFAULT_SIZE 10
#define RESET false

//////////////////////////////
//////////////////////////////

struct Options {
    // Testing options;
    bool debug = DEBUG;
    uint num_iter = NUM_ITER;
    uint N = DEFAULT_SIZE;
    bool reset = RESET;

    //////////////////////////////
    //////////////////////////////

    Options(int argc, char *argv[]) {

        int opt;
        static struct option long_options[] = {{"debug", no_argument, 0, 'd'},
                                               {"num_iter", required_argument, 0, 'i'},
                                               {"N", required_argument, 0, 'n'},
                                               {"reset", no_argument, 0, 'r'},
                                               {0, 0, 0, 0}};
        // getopt_long stores the option index here;
        int option_index = 0;

        while ((opt = getopt_long(argc, argv, "di:n:r", long_options, &option_index)) != EOF) {
            switch (opt) {
                case 'd':
                    debug = true;
                    break;
                case 'i':
                    num_iter = atoi(optarg);
                    break;
                case 'n':
                    N = atoi(optarg);
                    break;
                case 'r':
                    reset = true;
                    break;
                default:
                    break;
            }
        }
    }
};
