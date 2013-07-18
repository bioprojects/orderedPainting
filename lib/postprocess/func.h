#ifndef FUNC_H
#define FUNC_H

#include <iostream>
#include <map>
#include <string>
#include <stdio.h>

using namespace std;

bool has_key_int2string(map<int, string> &m, int n);
bool has_key_string2int(map<string, int> &m, string n);
bool has_pairkey_int2double(map<pair<int, int>, double> &m, pair<int, int> n);
bool has_pairkey_int2string(map<pair<int, int>, string> &m, pair<int, int> n);
bool has_pairkey_intstring2int(map<pair<int, string>, int> &m, pair<int, string> n);
double square(double x);

#endif
