#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

int main(int argc, char* argv[]){

  ifstream in(argv[1]);
  ofstream out(argv[2]);
  int cutoff = stoi(argv[3]);
  int start, end;
  string line, bit;
  while(!in.eof()){
    getline(in, line);
    stringstream ss(line);
    ss >> bit;
    ss >> bit;
    start = stoi(bit);
    ss >> bit;
    end = stoi(bit);
    if ((end - start) > cutoff) out << line << endl;
  }
  return 0;
}
