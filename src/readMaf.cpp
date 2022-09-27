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
  std::string seqname1;
  std::ifstream inFileStream;
  std::istream& input = openIn(inputFileName, inFileStream);
  std::string line;
  while(getline(inFileStream, line)) {
    if (line[0] != '#') {
      std::stringstream(line) >> seqname1;
      seqnames1.push_back(seqname1);
    }
  }
  return seqnames1;
}
