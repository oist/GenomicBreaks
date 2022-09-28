#include <Rcpp.h>

#include <fstream>
#include <iostream>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace Rcpp;

//' Read a MAF file
//'
//' Reads a pairwise genome alignment in MAF format.  The file can be plain
//' text or compressed with `gzip`.
//'
//' Known limitations: Does not expand shell metacharacters.  Trusts blindly
//' file extension to determine compression.  Does not perform any validation on
//' the file format.  Assumes that the score comes first in the 'a' lines.
//'
//' @param inputFileName The name of the file to read
//' @return a `GenomicBreaks` object.
//' @importFrom Rcpp evalCpp
//' @useDynLib GenomicBreaks, .registration = TRUE
//' @export
// [[Rcpp::export]]

Rcpp::List readMAF (std::string inputFileName) {
  Rcpp::CharacterVector seqnames1;
  Rcpp::CharacterVector seqnames2;
  Rcpp::IntegerVector seqlengths1;
  Rcpp::IntegerVector seqlengths2;
  Rcpp::CharacterVector strands;
  Rcpp::IntegerVector scores;
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

  std::ifstream file(inputFileName, std::ios_base::in | std::ios_base::binary);
  boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
  if(inputFileName.substr(inputFileName.find_last_of(".") + 1) == "gz") {
    in.push(boost::iostreams::gzip_decompressor());
  }
  in.push(file);
  std::istream incoming(&in);

  std::string line;
  while(getline(incoming, line)) {
    if (line[0] != '#') {
      // do nothing
    }
    if (line[0] == 'a') {
      // Example line: a score=762 mismap=1e-09
      // I assume score comes first after 'a'.
      std::string first;
      std::string second;
      std::stringstream(line) >> first >> second;
      std::string delimiter = "=";
      std::string token = second.substr(second.find(delimiter) + 1);
      scores.push_back(std::stoi(token));
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
                                             Named("strand") = strands,
                                             Named("scores") = scores);
  return outputList;
}
