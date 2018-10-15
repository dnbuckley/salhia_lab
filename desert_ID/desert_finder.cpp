#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "utils.h"

using namespace std;

int chrToInt(string chrom){
  string chromNumStr = chrom.substr(3, chrom.length() - 1);
  return std::stoi(chromNumStr);
}

vector<string> tsvLines(vector<pmdCall> pmdCalls){
  vector<string> lines;
  ostringstream os;
  os << "chr" << '\t' << "start" << '\t' << "end" << '\t' << "length" << '\t'
  << "avgMeth" << '\t' <<" stdDevMeth" << '\t' << "numCpGs" << '\t'
  << "raw score" << '\t' << "Z score" <<  endl;
  lines.push_back(os.str());
  for (size_t i = 0; i < pmdCalls.size(); i++){
    lines.push_back(pmdCalls[i].printTSV());
  }
  return lines;
}

vector<string> bedLines(vector<pmdCall> pmdCalls){
  vector<string> lines;
  for(size_t i = 0; i < pmdCalls.size(); i++)
    lines.push_back(pmdCalls[i].printBED());
  return lines;
}
// first = pmd, second = meth file, third = outfile
int main(int argc, char* argv[]){

  ifstream pmdIn(argv[1]);
  ifstream methIn(argv[2]);
  ofstream out(argv[3]);
  vector<pmdCall> pmdCalls;
  //vector<string> methVals;
  string bit, line, chrCheck;
  vector<methRegion> wholeGenome;
  int counter = 0;

  while(!pmdIn.eof()){ //read in pmds

    pmdCall temp;
    getline(pmdIn, line);
    //cout << line << endl;
    if (line == "" || line == "\n") continue;
    temp.addInfo(line);
    pmdCalls.push_back(temp);

  }

  cerr << "reading in methylation values" << endl;
  string prevChrom = "chr1";

  vector<float> methVals;
  vector<int> locations;
  int size = 0;
  cerr << "reading chr1" << endl;
  while(!methIn.eof()){ //read in meth values (whole file)

    //if (counter%1000000 == 0) cerr << "on line " << counter << endl;
    getline(methIn, line);
    //cout <<  line << endl;
    stringstream ss(line);
    ss >> bit;
    if (line == "" || line == "\n" || bit == "chrX"
    || bit == "chrY" || bit == "chrM"){
      methRegion chrBuild;
      chrBuild.chr = chrToInt(prevChrom);
      chrBuild.methVals = methVals;
      chrBuild.location = locations;
      chrBuild.size = size;
      wholeGenome.push_back(chrBuild);
      break;
    }
    if (bit != prevChrom){
      methRegion chrBuild;
      chrBuild.chr = chrToInt(prevChrom);
      chrBuild.methVals = methVals;
      chrBuild.location = locations;
      chrBuild.size = size;
      wholeGenome.push_back(chrBuild);
      methVals.clear();
      locations.clear();
      size = 0;
      cerr << "reading " << bit << endl;
    }
    prevChrom = bit;
    ss >> bit;
    locations.push_back(stoi(bit));
    ss >> bit;
    ss >> bit;
    methVals.push_back(stof(bit));
    counter++;
    size++;

  }

  //cerr << "num of vals in wholeGenome: " << wholeGenome.meth.size() << endl;
  cerr << "number of pmds = " << pmdCalls.size() << endl;
  cerr << "number of chrs in wholeGenome = " << wholeGenome.size() << endl;
  int start, end;

  for (size_t i = 0; i < pmdCalls.size(); i++){
    for(size_t j = 0; j < wholeGenome.size(); j++){
      if(pmdCalls[i].chrom() == wholeGenome[j].chr){
        start = pmdCalls[i].returnStart(wholeGenome[j].location);
        end = pmdCalls[i].returnEnd(wholeGenome[j].location);
        pmdCalls[i].addMeth(wholeGenome[j].methVals, start, end);
      }
    }
  }
  assignZ(pmdCalls);
  displayStats(pmdCalls);
  float cutoff;
  if (argc == 4){
    cerr << "what length cutoff would you like to set?" << endl;
    cin >> cutoff;
  }
  if (argc == 5) cutoff = float(stof(argv[4]));
  pmdCalls = filterPMDsRawScore(pmdCalls, cutoff);
  pmdCalls = filterLength(pmdCalls, 200000);
  cout << "cutoff = " << cutoff << endl;
  string choice;
  //cerr << "bed or tsv?" << endl;
  //cin >> choice;
  choice = "bed"; //hard coded for now
  if (choice == "tsv"){
    vector<string> lines = tsvLines(pmdCalls);
    for (size_t i = 0; i < lines.size(); i++) out << lines[i];
  }

  if (choice == "bed"){
    vector<string> lines = bedLines(pmdCalls);
    for (size_t i = 0; i < lines.size(); i++) out << lines[i];
  }

  else cerr << "invalid choice" << endl;

  return 0;



/*
  if (start < 0){
    pmdCalls[i].printPMD();
    cout << "ERROR mismatch = " << wholeGenome[j].location[-start]-pmdCalls[i].start << endl;
    cout << "+- = " << wholeGenome[j].location[-start+1]-pmdCalls[i].start << ", "
    << wholeGenome[j].location[-start-1]-pmdCalls[i].start << endl;
    cout << "pmdCHR = " << pmdCalls[i].chrom() << ", methCHR = " << wholeGenome[j].chr << endl << endl;
    //cout << "locations"

  }
*/


}
