#ifndef FUNC_H
#define FUNC_H

#include <map>
#include <stdio.h>

#include <cerrno>

#include <iostream>
#include <string>
#include <sstream>

#include <stdlib.h>
#include "zlib.h"

using namespace std;

bool has_key_int2string(map<int, string> &m, int n);
bool has_key_string2int(map<string, int> &m, string n);
bool has_pairkey_int2double(map<pair<int, int>, double> &m, pair<int, int> n);
bool has_pairkey_int2string(map<pair<int, int>, string> &m, pair<int, int> n);
bool has_pairkey_intstring2int(map<pair<int, string>, int> &m, pair<int, string> n);
double square(double x);

// The following functions are kindly provided to the public domain by Dr. Timothee Flutre
// https://openwetware.org/wiki/User:Timothee_Flutre/Notebook/Postdoc/2012/09/12
void openFile (
  const string & pathToFile,
  gzFile & fileStream,
  const char * mode);

void closeFile (
  const string & pathToFile,
  gzFile & fileStream);

int getline (
  gzFile & fileStream,
  string & line);

#endif
