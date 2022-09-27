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
Rcpp::List readMAF (std::string inputFileName) {
  Rcpp::CharacterVector seqnames1;
  Rcpp::CharacterVector seqnames2;
  Rcpp::IntegerVector seqlengths1;
  Rcpp::IntegerVector seqlengths2;
  Rcpp::CharacterVector strands;
  Rcpp::IntegerVector start1;
  Rcpp::IntegerVector start2;
  Rcpp::IntegerVector length1;
  Rcpp::IntegerVector length2;
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
      // Example line: s BK006934.2   19127 1679 +  562643 AACCAATCCAAAA...
      std::stringstream(line) >> linetype >> seqname >> startpos >> length >> strand >> seqlength >> seq;
      seqnames1.push_back(seqname);
      seqlengths1.push_back(seqlength);
      start1.push_back(startpos);
      length1.push_back(length);
      // Assume that there are always two 's' lines in a row
      getline(inFileStream, line);
      std::stringstream(line) >> linetype >> seqname >> startpos >> length >> strand >> seqlength >> seq;
      seqnames2.push_back(seqname);
      seqlengths2.push_back(seqlength);
      start2.push_back(startpos);
      length2.push_back(length);
      strands.push_back(strand); // Assume that the strand of seq1 is always '+'
    }
  }
  Rcpp::List outputList = Rcpp::List::create(Named("seqnames1") = seqnames1,
                                             Named("seqlengths1") = seqlengths1,
                                             Named("start1") = start1,
                                             Named("length1") = length1,
                                             Named("seqnames2") = seqnames2,
                                             Named("seqlengths2") = seqlengths2,
                                             Named("start2") = start2,
                                             Named("length2") = length2,
                                             Named("strand") = strands);
  return outputList;
}
