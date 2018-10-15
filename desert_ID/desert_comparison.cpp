#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "utils.h"


vector<string> printOverlapBed(vector<overlapCall> overlapCalls){
  vector<string> lines;
  for (size_t i = 0; i < overlapCalls.size(); i++){
    ostringstream os;
    os << "chr" <<  overlapCalls[i].chr << '\t' << overlapCalls[i].start << '\t'
    << overlapCalls[i].end << '\t' << overlapCalls[i].overlaps << endl;
    lines.push_back(os.str());
  }
  return lines;
}

vector<string> printOverlapTSV(vector<overlapCall> overlapCalls){
  vector<string> lines;
  ostringstream os;
  os << "chr" << '\t' << "start" << '\t' << "end" << '\t' << "length" << '\t'
  << "overlaps called" << '\t' <<"mean overlap" << '\t' << "files" <<endl;
  lines.push_back(os.str());
  for (size_t i = 0; i < overlapCalls.size(); i++){
    float sum = 0;
    ostringstream fileNames;
    if (overlapCalls[i].percentOverlap.size() != overlapCalls[i].files.size()){
      cerr << "problem with overlapcall, unequal vectors.  " << "percentOverlap.size() = "
      << overlapCalls[i].percentOverlap.size() << "\tfiles.size() = "
      << overlapCalls[i].files.size() <<  endl;  //sanity check
      for (size_t j = 0; j < overlapCalls[i].percentOverlap.size(); j++){
        cerr << "percentOverlap = " << overlapCalls[i].percentOverlap[j] << endl;
      }
      for (size_t j = 0; j < overlapCalls[i].files.size(); j++){
        cerr << "percentOverlap = " << overlapCalls[i].files[j] << endl;
      }
      cerr << "mismatch in lengths" << endl;
    }
    for (size_t j = 0; j < overlapCalls[i].percentOverlap.size(); j++){
      sum += overlapCalls[i].percentOverlap[j];
      fileNames << overlapCalls[i].files[j] << " || ";
    }
    float avgOverlap = float(sum/overlapCalls[i].percentOverlap.size());
    ostringstream os;
    os << "chr" << overlapCalls[i].chr << '\t' << overlapCalls[i].start << '\t'
    << overlapCalls[i].end << '\t' << overlapCalls[i].end - overlapCalls[i].start
    << '\t' << overlapCalls[i].overlaps << '\t' << avgOverlap << '\t'
    << fileNames.str() << endl;
    lines.push_back(os.str());
  }
  return lines;
}

// arg 1 = template, arg 2 = outfile, arg 3-n sample files
// or arg 1 = outfile, arg 2-n = compareison files
int main (int argc, char* argv[]){

  char choice;
  vector<string> tsvLines;
  vector<string> bedLines;
  vector<pmdFile> pmdFiles;
  string line;
  int totalPMDs = 0;

  cerr << "compare all? y/n" << endl;
  cin >> choice;

  if (choice == 'n'){
    ifstream templateFile(argv[1]);
    ofstream out(argv[2]);
    pmdFile pmdTemplate;
    cerr << "template file = " << argv[1] << endl;

    while (!templateFile.eof()){
      getline(templateFile, line);
      if (line == "") break;
      pmdCall temp;
      temp.addInfo(line);
      pmdTemplate.pmdCalls.push_back(temp);
    }
    for (int i = 3; i < argc; i++){
      cerr << "reading " << argv[i] << endl;
      ifstream in(argv[i]);
      vector<pmdCall> tempCalls;
      while (!in.eof()){
        getline(in, line);
        if (line == "") break;
        pmdCall temp;
        temp.addInfo(line);
        temp.file = argv[i];
        tempCalls.push_back(temp);
        totalPMDs++;
      }
      pmdFile tempFile;
      tempFile.pmdCalls = tempCalls;
      tempFile.filename = argv[i]; // this may cause problems
      pmdFiles.push_back(tempFile);
    }
    // float minOverlap;
    // cerr << "min overlap? (percent)" << endl;
    // cin >> minOverlap;
    // minOverlap = minOverlap/100;
    //vector<overlapCall> overlapCalls;
    //overlapCalls = callOverlaps(pmdTemplate, pmdFiles, minOverlap);
    //tsvLines = printOverlapTSV(overlapCalls);
    //bedLines = printOverlapBed(overlapCalls);

    vector<pmdCall> newPMDs = findNewPmds(pmdTemplate, pmdFiles);

    cerr << "totalPMDs = " << totalPMDs << endl;
    cerr << "new pmds = " << newPMDs.size() << endl;

    for (size_t i = 0; i < newPMDs.size(); i++)  //uncomment for tsv
      out << "chr" << newPMDs[i].chr << '\t' << newPMDs[i].start << '\t'
      << newPMDs[i].end << '\t' << newPMDs[i].file << '\n';
    //
    // for (size_t i = 0; i < bedLines.size(); i++)  //uncomment for bed
    //   out << bedLines[i];
  }

  else if (choice == 'y'){
    ofstream out(argv[1]);
    for (int i = 2; i < argc; i++){
      cout << "reading " << argv[i] << endl;
      ifstream in(argv[i]);
      vector<pmdCall> tempCalls;
      while (!in.eof()){
        getline(in, line);
        if (line == "") break;
        pmdCall temp;
        temp.addInfo(line);
        tempCalls.push_back(temp);
      }
      pmdFile tempFile;
      tempFile.pmdCalls = tempCalls;
      tempFile.filename = argv[i]; // this may cause problems
      pmdFiles.push_back(tempFile);
    }
    float minOverlap;
    cerr << "min overlap? (percent)" << endl;
    cin >> minOverlap;
    minOverlap = minOverlap/100;
    vector<overlapCall> overlapCalls;
    overlapCalls = callOverlaps(pmdFiles, minOverlap);
    tsvLines = printOverlapTSV(overlapCalls);
    bedLines = printOverlapBed(overlapCalls);

    for (size_t i = 0; i < tsvLines.size(); i++)  //uncomment for tsv
      out << tsvLines[i];
    //
    // for (size_t i = 0; i < bedLines.size(); i++)  //uncomment for bed
    //   out << bedLines[i];
  }


  return 0;
}
