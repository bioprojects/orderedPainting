#ifdef WIN32
#include "getopt.h"
#else
#include <unistd.h>
#endif

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <sstream>
#include <iostream>
#include <fstream>

#include <map>
#include <string>
#include <vector>

#include <regex.h> 

//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_rng.h>

using namespace std;
//using namespace boost;

static const char * help=
        "\
        Usage: gzip -dc file.gz | sp [OPTIONS] \n\
            \n\
            Options:\n\
                -p outPrefix_ \n\
                \n";


// ######################################################################
// main 
// ######################################################################
int main(int argc, char **argv)
{

    // getopt
    int c;
    bool verbose=false;
    char * outPrefix=NULL;

    if (argc==1) {printf("%s",help);exit(0);}
    while ((c = getopt (argc, argv, "p:v")) != -1)
        switch (c)
    {
        case('p'):outPrefix=optarg;break;
        case '?':
            if (isprint (optopt))
                fprintf (stderr, "Unknown option `-%c'.\n", optopt);
            else
                fprintf (stderr,
                "Unknown option character `\\x%x'.\n",
                optopt);
            return 1;
        default:
            abort ();
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // variables
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    time_t timer;
    char *stamp;

    FILE *fh_out1;
    FILE *fh_out2;
    FILE *fh_out3;
    FILE *fh_out4;
    FILE *fh_out5;
    FILE *fh_out6;
    FILE *fh_out7;
    FILE *fh_out8;
    FILE *fh_out9;

    static const char *reg1 = "^(1)";
    static const char *reg2 = "^(2)";
    static const char *reg3 = "^(3)";
    static const char *reg4 = "^(4)";
    static const char *reg5 = "^(5)";
    static const char *reg6 = "^(6)";
    static const char *reg7 = "^(7)";
    static const char *reg8 = "^(8)";
    static const char *reg9 = "^(9)";

    regex_t regst1;
    regex_t regst2;
    regex_t regst3;
    regex_t regst4;
    regex_t regst5;
    regex_t regst6;
    regex_t regst7;
    regex_t regst8;
    regex_t regst9;
    
    regmatch_t match[1];

    char fname[512];
    
    string line;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // prepare regexpr
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (regcomp(&regst1, reg1, REG_EXTENDED)) {
        return 1;
    }
    if (regcomp(&regst2, reg2, REG_EXTENDED)) {
        return 1;
    }
    if (regcomp(&regst3, reg3, REG_EXTENDED)) {
        return 1;
    }
    if (regcomp(&regst4, reg4, REG_EXTENDED)) {
        return 1;
    }
    if (regcomp(&regst5, reg5, REG_EXTENDED)) {
        return 1;
    }
    if (regcomp(&regst6, reg6, REG_EXTENDED)) {
        return 1;
    }
    if (regcomp(&regst7, reg7, REG_EXTENDED)) {
        return 1;
    }
    if (regcomp(&regst8, reg8, REG_EXTENDED)) {
        return 1;
    }
    if (regcomp(&regst9, reg9, REG_EXTENDED)) {
        return 1;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // fopen
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    sprintf( fname, "%s%02d", outPrefix, 1);
    ofstream ofs1(fname);

    sprintf( fname, "%s%02d", outPrefix, 2); 
    ofstream ofs2(fname);

    sprintf( fname, "%s%02d", outPrefix, 3); 
    ofstream ofs3(fname);

    sprintf( fname, "%s%02d", outPrefix, 4); 
    ofstream ofs4(fname);

    sprintf( fname, "%s%02d", outPrefix, 5); 
    ofstream ofs5(fname);

    sprintf( fname, "%s%02d", outPrefix, 6); 
    ofstream ofs6(fname);

    sprintf( fname, "%s%02d", outPrefix, 7); 
    ofstream ofs7(fname);

    sprintf( fname, "%s%02d", outPrefix, 8); 
    ofstream ofs8(fname);

    sprintf( fname, "%s%02d", outPrefix, 9); 
    ofstream ofs9(fname);


    // ########################################################################################

    //
    // calculate an average matrix by processing each site in gz_sortfile_catOrderings_copyprob
    //
    timer = time(NULL); stamp = ctime(&timer); stamp[strlen(stamp)-1] = '\0';
    printf("%s: start\n", stamp);

    while (getline(cin, line)) {

        if (!regexec(&regst1, line.c_str(), 1, match, 0)) {
            ofs1 << line << endl;
        } else if (!regexec(&regst2, line.c_str(), 1, match, 0)) {
            ofs2 << line << endl;
        } else if (!regexec(&regst3, line.c_str(), 1, match, 0)) {
            ofs3 << line << endl;
        } else if (!regexec(&regst4, line.c_str(), 1, match, 0)) {
            ofs4 << line << endl;
        } else if (!regexec(&regst5, line.c_str(), 1, match, 0)) {
            ofs5 << line << endl;
        } else if (!regexec(&regst6, line.c_str(), 1, match, 0)) {
            ofs6 << line << endl;
        } else if (!regexec(&regst7, line.c_str(), 1, match, 0)) {
            ofs7 << line << endl;
        } else if (!regexec(&regst8, line.c_str(), 1, match, 0)) {
            ofs8 << line << endl;
        } else if (!regexec(&regst9, line.c_str(), 1, match, 0)) {
            ofs9 << line << endl;
        } else {
           // do nothing for "pos" or "HAP"
        }

    } // while read line

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // finish
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    timer = time(NULL); stamp = ctime(&timer); stamp[strlen(stamp)-1] = '\0';
    printf("%s: finished \n", stamp);

    return 0;


}
