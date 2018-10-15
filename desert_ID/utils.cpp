#include "utils.h"

pmdCall::pmdCall(){}
pmdCall::~pmdCall(){}

void pmdCall::addInfo(string line){

  string bit;
  stringstream ss(line);
  ss >> bit;
  chr = stoi(bit.substr(3, bit.length()-1));
  ss >> bit;
  start = stoi(bit);
  ss >> bit;
  end = stoi(bit);
  length = end - start;
  return;

}

void pmdCall::addMeth(vector<float> meth, int start, int end){
  long double sum = 0;
  bool biscuit = false;
  int n = 0;
  for (int i = start; i <= end; i++){
    methVals.push_back(meth[i]);
    sum += meth[i];
    n++;
  }
  avgMeth = float(sum/n);
  if (avgMeth < 1){
    biscuit = true;
    avgMeth = avgMeth*100;
  }
  numCpGs = n;
  sum = 0;
  for (size_t i = 0; i < methVals.size(); i++){
    if (biscuit) sum += pow((methVals[i]*100)-avgMeth, 2);
    else sum += pow(methVals[i]-avgMeth, 2);
  }
  stdDevMeth = sqrt(sum/(n-1));

  //add score here??

  //score = length/(avgMeth*stdDevMeth); //first gen, cutoff ~1000-2000
  //score = sqrt(length)/(avgMeth*stdDevMeth); //second gen, cutoff ~2-4
  score = 1/(avgMeth*stdDevMeth); //third gen, cutoff
  //score = length/(stdDevMeth);
  //score = log(length)/(avgMeth*stdDevMeth); //NO


  return;
}

void pmdCall::methDropoff(vector<string> regio){
  return;
}

int pmdCall::returnStart(vector<int> locations){
  int m;
  int hi = locations.size() - 1;
  int lo = 0;
  while (lo <= hi){
    m = (lo+hi)/2;
    if(start == locations[m]){
      startMismatch = 0;
      return m;
    }
    if(start < locations[m]) hi = m-1;
    else lo = m+1;
  }
  if (abs(locations[m+1]-start) < abs(locations[m]-start)) m++;
  if (abs(locations[m-1]-start) < abs(locations[m]-start)) m--;
  startMismatch = start - locations[m];
  return m;
}

int pmdCall::returnEnd(vector<int> locations){
  int m;
  int hi = locations.size() - 1;
  int lo = 0;
  while (lo <= hi){
    m = (lo+hi)/2;
    if(end == locations[m]){
      endMismatch = 0;
      return m;
    }
    if(end < locations[m]) hi = m-1;
    else lo = m+1;
  }
  if (abs(locations[m+1]-end) < abs(locations[m]-end)) m++;
  if (abs(locations[m-1]-end) < abs(locations[m]-end)) m--;
  endMismatch = end - locations[m];
  return m;
}

int pmdCall::chrom(){return chr;}

void pmdCall::printPMD(){
  cout << "chr = " << chr << endl
  << "start = " << start << endl
  << "end = " << end << endl
  << "length = " << length << endl;
  return;
}


string pmdCall::printTSV(){
  ostringstream os;
  os << chr << '\t' << start << '\t' << end << '\t' << length << '\t'
  << avgMeth << '\t' << stdDevMeth << '\t' << numCpGs << '\t' << score
  << '\t' << zScore << endl;
  return os.str();
}

string pmdCall::printBED(){
  ostringstream os;
  os << "chr" <<  chr << '\t' << start << '\t' << end << '\t' << zScore << endl;
  return os.str();
}

void assignZ(vector<pmdCall> &pmdCalls){

  long double sum = 0;
  for (size_t i = 0; i < pmdCalls.size(); i++) sum += pmdCalls[i].score;
  double meanScore = sum/pmdCalls.size();
  double var = 0;
  for (size_t i = 0; i < pmdCalls.size(); i++){
    var += pow(pmdCalls[i].score - meanScore, 2);
  }
  double stdDev = sqrt(var/(pmdCalls.size()-1));

  for (size_t i = 0; i < pmdCalls.size(); i++){
    pmdCalls[i].zScore = (pmdCalls[i].score-meanScore)/stdDev;
  }
  return;
}

vector <pmdCall> filterPMDsZscore(vector<pmdCall> pmdCalls, float cutoff){
  vector<pmdCall> filteredCalls;
  for (size_t i = 0; i < pmdCalls.size(); i++)
    if (pmdCalls[i].zScore >= cutoff)
      filteredCalls.push_back(pmdCalls[i]);
  return filteredCalls;
}

vector<pmdCall> filterPMDsRawScore(vector<pmdCall> pmdCalls, float cutoff){
  vector<pmdCall> filteredCalls;
  if (cutoff == -1){

  }
  else{
    for (size_t i = 0; i < pmdCalls.size(); i++)
      if (pmdCalls[i].score >= cutoff)
        filteredCalls.push_back(pmdCalls[i]);
  }
  return filteredCalls;
}

vector <pmdCall> filterLength(vector<pmdCall> pmdCalls, int length){
  vector<pmdCall> filteredCalls;
  for (size_t i = 0; i < pmdCalls.size(); i++)
    if (pmdCalls[i].length >= length)
      filteredCalls.push_back(pmdCalls[i]);
  return filteredCalls;
}


void displayStats(vector<pmdCall> pmdCalls){

  vector<int> lengths;
  vector<float> meths;
  vector<float> stdDevs;
  vector<double> scores;
  long double totalLen = 0;
  long double totalMeth = 0;
  long double totalSD = 0;
  long double totalScore = 0;
  int n = 0;
  float avgLen, avgMeth, avgSD, stdDevLen, stdDevMeth, stdDevSD, varLen;
  float varMeth, varSD;
  float avgScore, varScore, scoreSD;

  for (size_t i = 0; i < pmdCalls.size(); i++){
    lengths.push_back(pmdCalls[i].length);
    meths.push_back(pmdCalls[i].avgMeth);
    stdDevs.push_back(pmdCalls[i].stdDevMeth);
    scores.push_back(pmdCalls[i].score);
    totalLen += pmdCalls[i].length;
    totalMeth += pmdCalls[i].avgMeth;
    totalSD += pmdCalls[i].stdDevMeth;
    totalScore += pmdCalls[i].score;
    n++;
  }
  avgLen = totalLen/n;
  avgMeth = totalMeth/n;
  avgSD = totalSD/n;
  avgScore = totalScore/n;

  sort(lengths.begin(), lengths.end());
  sort(meths.begin(), meths.end());
  sort(stdDevs.begin(), stdDevs.end());
  sort(scores.begin(), scores.end());

  for (size_t i = 0; i < pmdCalls.size(); i++){
    varLen += pow(lengths[i]-avgLen, 2);
    varMeth += pow(meths[i]-avgMeth, 2);
    varSD += pow(stdDevs[i]-avgSD, 2);
    varScore += pow(scores[i]-avgScore, 2);
  }

  stdDevLen = sqrt(varLen/(n-1));
  stdDevMeth = sqrt(varMeth/(n-1));
  stdDevSD = sqrt(varSD/(n-1));
  scoreSD = sqrt(varScore/(n-1));

  cerr << "number of pmds = " << pmdCalls.size() << endl;

  cerr << "mean length (kb) = " << avgLen/1000 << ", " << "SD = "
        << stdDevLen/1000 << endl;
  cerr << "Q1 = " << floor(lengths[floor(lengths.size()*.25)]/1000) << endl;
  cerr << "Q2 = " << floor(lengths[floor(lengths.size()*.5)]/1000) << endl;
  cerr << "Q3 = " << floor(lengths[floor(lengths.size()*.75)]/1000) << endl;

  cerr << "mean meth = " << avgMeth << ", " << "SD =" << stdDevMeth << endl;
  cerr << "Q1 = " << meths[floor(meths.size()*.25)] << endl;
  cerr << "Q2 = " << meths[floor(meths.size()*.5)] << endl;
  cerr << "Q3 = " << meths[floor(meths.size()*.75)] << endl;

  cerr << "mean SD = " << avgSD << ", " << "SD = " << stdDevSD << endl;
  cerr << "Q1 = " << stdDevs[floor(stdDevs.size()*.25)] << endl;
  cerr << "Q2 = " << stdDevs[floor(stdDevs.size()*.5)] << endl;
  cerr << "Q3 = " << stdDevs[floor(stdDevs.size()*.75)] << endl;

  cerr << "mean score = " << avgScore << ", " << "SD = " << scoreSD << endl;
  cerr << "Q1 = " << scores[floor(scores.size()*.25)] << endl;
  cerr << "Q2 = " << scores[floor(scores.size()*.5)] << endl;
  cerr << "Q3 = " << scores[floor(scores.size()*.75)] << endl;

}

int matchExists(vector<overlapCall> overlapCalls, int start, int end, int chr){
  for (size_t i = 0; i < overlapCalls.size(); i++){
    if (overlapCalls[i].start == start && overlapCalls[i].end == end
        && overlapCalls[i].chr == chr) return i;
  }
  return -1;
}

bool containedBy(pmdCall refPmd, pmdCall tempPmd){ // may want to take out second contains
  return ((refPmd.start <= tempPmd.start && refPmd.end >= tempPmd.end) ||
          (tempPmd.start <= refPmd.start && tempPmd.end >= refPmd.end));
}

bool doesOverlap(pmdCall refPmd, pmdCall tempPmd, float minOverlap){ // this might be fucked
  float overlap = 0;
  float refLen = refPmd.end - refPmd.start;
  if (refPmd.start <= tempPmd.start && refPmd.end <= tempPmd.end &&
      refPmd.end >= tempPmd.start){
    overlap = float((refLen - (tempPmd.start - refPmd.start))/refLen);
    // cout << "overlap1 = " << overlap << "\treflen = " << refLen
    //  << "\tminOverlap = " << minOverlap << "\t tempPmd.start = " << tempPmd.end
    //  << "\trefPmd.start = "<< refPmd.start <<  endl;
  }
  if (refPmd.start >= tempPmd.start && refPmd.end >= tempPmd.end &&
          refPmd.start <= tempPmd.end){
    overlap = float((refLen - (refPmd.end - tempPmd.end))/refLen);
    // cout << "overlap2 = " << overlap << "\treflen = " << refLen
    //  << "\tminOverlap = " << minOverlap << "\t refPmd.end = " << refPmd.end
    //  << "\ttempPmd.end = "<< tempPmd.end << endl;
  }

  //else return false;
  return (overlap > minOverlap);
}

float overlapPercent(pmdCall refPmd, pmdCall tempPmd){
  float refLen = refPmd.end - refPmd.start;
  if (refPmd.start <= tempPmd.start && refPmd.end <= tempPmd.end)
    return float((refLen - (tempPmd.start - refPmd.start))/refLen);
  if (refPmd.start >= tempPmd.start && refPmd.end >= tempPmd.end){
    //cout << "percent = " << float((refLen - (refPmd.end - tempPmd.end))/refLen) << endl;
    return float((refLen - (refPmd.end - tempPmd.end))/refLen);
  }
  else{
    cerr << "error in overlapPercent" << endl;
    return -1;
  }
}

vector<overlapCall> callOverlaps(pmdFile pmdTemplate, vector<pmdFile> pmdFiles,
                                  float minOverlap){ //may need to restructure like below
  vector<overlapCall> overlapCalls;
  for (size_t i = 0; i < pmdTemplate.pmdCalls.size(); i++){
    pmdCall refPmd = pmdTemplate.pmdCalls[i];
    for (size_t j = 0; j < pmdFiles.size(); j++){
      for (size_t m = 0; m < pmdFiles[j].pmdCalls.size(); m++){
        pmdCall tempPmd = pmdFiles[j].pmdCalls[m];
        if (refPmd.chr == tempPmd.chr){
          if (containedBy(refPmd, tempPmd) ||
              doesOverlap(refPmd, tempPmd, minOverlap)){
            int idx = matchExists(overlapCalls, refPmd.start, refPmd.end, refPmd.chr);
            if (idx == -1){
              overlapCall newCall;
              newCall.overlaps = 1;
              newCall.start = refPmd.start;
              newCall.end = refPmd.end;
              newCall.chr = refPmd.chr;
              newCall.calls.push_back(tempPmd);
              newCall.files.push_back(pmdFiles[j].filename);
              newCall.startMismatch.push_back(tempPmd.start - refPmd.start);
              newCall.endMismatch.push_back(tempPmd.end - refPmd.end);
              if (containedBy(refPmd, tempPmd))
                newCall.percentOverlap.push_back(100);
              else newCall.percentOverlap.push_back(overlapPercent(refPmd, tempPmd));
              overlapCalls.push_back(newCall);
            }
            else{
              overlapCalls[idx].calls.push_back(tempPmd);
              overlapCalls[idx].files.push_back(pmdFiles[j].filename);
              overlapCalls[idx].startMismatch.push_back(tempPmd.start - refPmd.start);
              overlapCalls[idx].endMismatch.push_back(tempPmd.end - refPmd.end);
              if (containedBy(refPmd, tempPmd))
                overlapCalls[idx].percentOverlap.push_back(100);
              else overlapCalls[idx].percentOverlap.push_back(overlapPercent(refPmd, tempPmd));
              overlapCalls[idx].overlaps++;
            }
          }
        }
      }
    }
  }
  return overlapCalls;
}


vector<overlapCall> callOverlaps(vector<pmdFile> pmdFiles, float minOverlap){ //may need to restructure like below

  vector<overlapCall> overlapCalls;
  for (size_t p = 0; p < pmdFiles.size(); p++){
    pmdFile pmdTemplate = pmdFiles[p];
    for (size_t i = 0; i < pmdTemplate.pmdCalls.size(); i++){
      pmdCall refPmd = pmdTemplate.pmdCalls[i];
      for (size_t j = p+1; j < pmdFiles.size(); j++){
        for (size_t m = 0; m < pmdFiles[j].pmdCalls.size(); m++){
          pmdCall tempPmd = pmdFiles[j].pmdCalls[m];
          if (refPmd.chr == tempPmd.chr){
            if (containedBy(refPmd, tempPmd) ||
                doesOverlap(refPmd, tempPmd, minOverlap)){
              int idx = matchExists(overlapCalls, refPmd.start, refPmd.end, refPmd.chr);
              if (idx == -1){
                overlapCall newCall;
                newCall.overlaps = 1;
                newCall.start = refPmd.start;
                newCall.end = refPmd.end;
                newCall.chr = refPmd.chr;
                newCall.calls.push_back(tempPmd);
                newCall.files.push_back(pmdFiles[j].filename);
                newCall.startMismatch.push_back(tempPmd.start - refPmd.start);
                newCall.endMismatch.push_back(tempPmd.end - refPmd.end);
                if (containedBy(refPmd, tempPmd))
                  newCall.percentOverlap.push_back(100);
                else newCall.percentOverlap.push_back(overlapPercent(refPmd, tempPmd));
                overlapCalls.push_back(newCall);
              }
              else{
                overlapCalls[idx].calls.push_back(tempPmd);
                overlapCalls[idx].files.push_back(pmdFiles[j].filename);
                overlapCalls[idx].startMismatch.push_back(tempPmd.start - refPmd.start);
                overlapCalls[idx].endMismatch.push_back(tempPmd.end - refPmd.end);
                if (containedBy(refPmd, tempPmd))
                  overlapCalls[idx].percentOverlap.push_back(100);
                else overlapCalls[idx].percentOverlap.push_back(overlapPercent(refPmd, tempPmd));
                overlapCalls[idx].overlaps++;
              }
            }
          }
        }
      }
    }
  }
  return overlapCalls;
}


vector<pmdCall> findNewPmds(pmdFile pmdTemplate, vector<pmdFile> pmdFiles){ // TODO fix this

  vector<pmdCall> newPMDs;
  bool found;
  for (size_t j = 0; j < pmdFiles.size(); j++){
    for (size_t m = 0; m < pmdFiles[j].pmdCalls.size(); m++){
      pmdCall tempPmd = pmdFiles[j].pmdCalls[m];
      found = false;
      for (size_t i = 0; i < pmdTemplate.pmdCalls.size(); i++){
        pmdCall refPmd = pmdTemplate.pmdCalls[i];
        if (refPmd.chr == tempPmd.chr){
          if (containedBy(refPmd, tempPmd) || doesOverlap(refPmd, tempPmd, 0)){
            found = true;
            // cout << "ref:\t" << refPmd.chr << '\t' << refPmd.start << '\t' << refPmd.end
            // << '\t' << containedBy(refPmd, tempPmd) << '\t' << doesOverlap(refPmd, tempPmd, 0) << endl;
            // cout << "temp:\t" << tempPmd.chr << '\t' << tempPmd.start
            // << '\t' << tempPmd.end << endl;
          }
        }
      }
      if (!found) newPMDs.push_back(tempPmd);
    }
  }
  return newPMDs;
}


// if (abs(refPmd.start - tempPmd.start) < maxGap &&
//     abs(refPmd.end - tempPmd.end) < maxGap &&
//     refPmd.chr == tempPmd.chr){
