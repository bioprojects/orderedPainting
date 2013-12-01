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

double square(double x) {
    return x*x;
}

void openFile (
  const string & pathToFile,
  gzFile & fileStream,
  const char * mode)
{
  fileStream = gzopen (pathToFile.c_str(), mode);
  if (fileStream == NULL)
  {
    cerr << "ERROR: can't open file " << pathToFile
         << " with mode " << *mode
         << " (errno=" << errno << ")" << endl;
    exit (1);
  }
}

void closeFile (
  const string & pathToFile,
  gzFile & fileStream)
{
  int ret = gzclose (fileStream);
  if (ret != Z_OK)
  {
    cerr << "ERROR: can't close the file " << pathToFile
         << ", gzclose() returned " << ret << endl;
    exit (1);
  }
}

int getline (
  gzFile & fileStream,
  string & line)
{
  int res = 1, c;
  line.clear ();
  while (true)
  {
    c = gzgetc (fileStream);
    if (c == -1) // eof or error
    {
      res = 0;
      break;
    }
    else if (c == 10) // 10 is ASCII code for '\n'
      break;
    else
      line.push_back (c);
  }
  return res;
}