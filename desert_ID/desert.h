#ifndef DESERT_H_
#define DESERT_H_
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

struct methRegion{
  vector<int> chr;
  vector<float> meth;
  vector<int> location;
};

class pmdCall{
  public:
    pmdCall();
    ~pmdCall();

    void addInfo(string line);

    void addMeth(vector<string> region);

    void methDropoff(vector<string> region);

  private:
    int chr, start, end, length;
    float avgMeth, stdDevMeth, methDropoffHead, methDropoffTail; //normal/tumor .7-.8 meth value
    methRegion methVals;

};

#endif
