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
Rcpp::NumericVector readMAF (std::string inputFileName) {
  std::ifstream inFileStream;
  std::istream& input = openIn(inputFileName, inFileStream);
  return 1;
}
