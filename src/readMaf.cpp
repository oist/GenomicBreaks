#include <Rcpp.h>

#include <string.h>

#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <streambuf>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace Rcpp;
//' Read a MAF file
//'
//' @param inputFileName The name of the file to read
//' @return a GenomicBreaks object
//' @export

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
  int score;

  std::ifstream file(inputFileName, std::ios_base::in | std::ios_base::binary);
  boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
  //in.push(boost::iostreams::gzip_decompressor());
  in.push(file);
  std::istream incoming(&in);

  std::string line;
  while(getline(incoming, line)) {
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
      getline(incoming, line);
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
