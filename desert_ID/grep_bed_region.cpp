#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

using namespace std;

int chrToInt(string chrom){
  string chromNumStr = chrom.substr(3, chrom.length() - 1);
  return stoi(chromNumStr);
}

int main(int argc, char* argv[]){

  ifstream in(argv[1]);
  ofstream out(argv[2]);
  string line, bit;
  int chr, start, end, chrx, startx, endx;

  cerr << "chr, start, end?" << endl;
  cin >> chr;
  cin >> start;
  cin >> end;

  while (!in.eof()){
    getline(in, line);
    stringstream ss(line);
    ss >> bit;
    chrx = chrToInt(bit);
    if (chrx != chr) continue;
    ss >> bit;
    startx = stoi(bit);
    ss >> bit;
    endx = stoi(bit);
    if (startx >= start && endx <= end)
      out << line << '\n';
  }
  return 0;
}
