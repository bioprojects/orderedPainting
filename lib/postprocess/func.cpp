#include "func.h"

using namespace std;

bool has_key_int2string(map<int, string> &m, int n){
    if (m.count(n) == 0){
        return false;
    }
    else{
        return true;
    }
}

bool has_key_string2int(map<string, int> &m, string n){
    if (m.count(n) == 0){
        return false;
    }
    else{
        return true;
    }
}

bool has_pairkey_int2double(map<pair<int, int>, double> &m, pair<int, int> n){
    if (m.count(n) == 0){
        return false;
    }
    else{
        return true;
    }
}

bool has_pairkey_int2string(map<pair<int, int>, string> &m, pair<int, int> n){
    if (m.count(n) == 0){
        return false;
    }
    else{
        return true;
    }
}

bool has_pairkey_intstring2int(map<pair<int, string>, int> &m, pair<int, string> n){
    if (m.count(n) == 0){
        return false;
    }
    else{
        return true;
    }
}

double square(double x)
{
    return x*x;
}
