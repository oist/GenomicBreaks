#include <Rcpp.h>

#include <string.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <streambuf>

using namespace Rcpp;
//' Read a MAF file
//'
//' @param inputFileName The name of the file to read
//' @return a GenomicBreaks object
//' @export


static void err(const std::string& s) {
  throw std::runtime_error(s);
}

static std::istream& openIn(const std::string& fileName, std::ifstream& ifs) {
  if (fileName == "-") return std::cin;
  ifs.open(fileName.c_str());
  if (!ifs) err("can't open file: " + fileName);
  return ifs;
}

// [[Rcpp::export]]
Rcpp::CharacterVector readMAF (std::string inputFileName) {
  Rcpp::CharacterVector seqnames1;
  Rcpp::CharacterVector seqnames2;
  std::string linetype;
  std::string seqname;
  int startpos;
  int length;
  std::string strand;
  int seqlength;
  std::string seq;
  std::ifstream inFileStream;
  std::istream& input = openIn(inputFileName, inFileStream);
  std::string line;
  int score;
  while(getline(inFileStream, line)) {
    if (line[0] != '#') {
      // do nothing
    }
    if (line[0] == 'a') {
      // dummy score for the moment
      score = 1;
    }
    if (line[0] == 's') {
      // Assume that there are only two 's' lines in a row
      // Example line: s BK006934.2   19127 1679 +  562643 AACCAATCCAAAA...
      std::stringstream(line) >> linetype >> seqname >> startpos >> length >> strand >> seqlength >> seq;
      seqnames1.push_back(seqname);
      getline(inFileStream, line);
      std::stringstream(line) >> linetype >> seqname >> startpos >> length >> strand >> seqlength >> seq;
      seqnames2.push_back(seqname);
    }
  }
  return seqnames1;
}
