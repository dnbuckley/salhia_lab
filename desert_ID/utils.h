#ifndef DESERT_H_
#define DESERT_H_
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>



using namespace std;

struct methRegion{
  int chr;
  vector<float> methVals;
  vector<int> location;
  int size;
};

class pmdCall{
  public:
    pmdCall();
    ~pmdCall();
    void addInfo(string line);
    void addMeth(vector<float> methVals, int start, int end);
    void methDropoff(vector<string> region);
    int returnStart(vector<int> locations);
    int returnEnd(vector<int> locations);
    int chrom();
    void printPMD();
    string printTSV();
    string printBED();

    double score, zScore;
    int chr, start, end, length, startMismatch, endMismatch, numCpGs;
    float avgMeth, stdDevMeth, methDropoffHead, methDropoffTail; //normal/tumor .7-.8 meth value
    vector<float> methVals;
    string label, file;
    char strand;

};

struct pmdFile{
  vector<pmdCall> pmdCalls;
  string filename;
};

struct overlapCall{
  int overlaps, start, end, chr;
  vector<int> startMismatch;
  vector<int> endMismatch;
  vector<float> percentOverlap;
  vector<string> files;
  vector<pmdCall> calls;
};

void assignZ(vector<pmdCall> &pmdCalls);
void displayStats(vector<pmdCall> pmdCalls);
vector<pmdCall> filterPMDsZscore(vector<pmdCall> pmdCalls, float cutoff);
vector<pmdCall> filterPMDsRawScore(vector<pmdCall> pmdCalls, float cutoff);
vector <pmdCall> filterLength(vector<pmdCall> pmdCalls, int length);
vector<overlapCall> callOverlaps(pmdFile pmdTemplate, vector<pmdFile> pmdFiles,
                                  float minOverlap);
vector<overlapCall> callOverlaps(vector<pmdFile> pmdFiles, float minOverlap);
vector<pmdCall> findNewPmds(pmdFile pmdTemplate, vector<pmdFile> pmdFiles);
int matchExists(vector<overlapCall> overlapCalls, int start, int end, int chr);
bool containedBy(pmdCall refPmd, pmdCall tempPmd);
bool doesOverlap(pmdCall refPmd, pmdCall tempPmd, float minOverlap);
float overlapPercent(pmdCall refPmd, pmdCall tempPmd);


#endif
