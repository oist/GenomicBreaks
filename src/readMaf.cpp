#include <Rcpp.h>

#include <fstream>
#include <iostream>
#include "mcf_zstream.h"
#include <string>
#include <cctype> // for std::toupper

using namespace Rcpp;

// ChatGPT answer for replacing boost::to_upper for ASCII input.
void toUpperAscii(std::string &str) {
  for (char &c : str) {
    c = std::toupper(static_cast<unsigned char>(c));
  }
}

//' Read a MAF file
//'
//' Reads a pairwise genome alignment in MAF format.  The file can be plain
//' text or compressed with `gzip`.
//'
//' Known limitations: Does not expand shell metacharacters.  Trusts blindly
//' file extension to determine compression.  Does not perform any validation on
//' the file format.  Assumes that the score comes first in the 'a' lines.
//'
//' @family Data loading functions
//'
//' @param inputFileName The name of the file to read
//' @return a `list` object with containg coordinates of the alignments in both
//' genomes and other information such as alignment width and number of matches.
//' @importFrom Rcpp evalCpp
//' @useDynLib GenomicBreaks, .registration = TRUE
// [[Rcpp::export]]
Rcpp::List readMAF (std::string inputFileName) {
  Rcpp::CharacterVector seqnames1;
  Rcpp::CharacterVector seqnames2;
  Rcpp::IntegerVector seqlengths1;
  Rcpp::IntegerVector seqlengths2;
  Rcpp::CharacterVector strand1;
  Rcpp::CharacterVector strand2;
  Rcpp::IntegerVector scores;
  Rcpp::IntegerVector start1;
  Rcpp::IntegerVector start2;
  Rcpp::IntegerVector length1;
  Rcpp::IntegerVector length2;
  Rcpp::IntegerVector matches; // Number of matches (case-insensitive)
  Rcpp::IntegerVector aLength; // Alignment length

  std::string linetype;
  std::string seqname;
  int startpos;
  int length;
  std::string strand;
  int seqlength;
  std::string seq1;
  std::string seq2;

  mcf::izstream incoming(inputFileName.c_str());

  if (!incoming.is_open()) {
    Rcpp::Rcerr << "Failed to open file" << std::endl;
    return 1;
  }

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
      std::stringstream(line) >> linetype >> seqname >> startpos >> length >> strand >> seqlength >> seq1;
      seqnames1.push_back(seqname);
      seqlengths1.push_back(seqlength);
      start1.push_back(startpos);
      length1.push_back(length);
      strand1.push_back(strand);
      // Assume that there are always two 's' lines in a row
      getline(incoming, line);
      std::stringstream(line) >> linetype >> seqname >> startpos >> length >> strand >> seqlength >> seq2;
      seqnames2.push_back(seqname);
      seqlengths2.push_back(seqlength);
      start2.push_back(startpos);
      length2.push_back(length);
      strand2.push_back(strand);
      // Compute the number of matches in the alignment
      toUpperAscii(seq1);
      toUpperAscii(seq2);
      Rcpp::CharacterVector seq1Vec(seq1.begin(),seq1.end());
      Rcpp::CharacterVector seq2Vec(seq2.begin(),seq2.end());
      auto matchesVec(seq1Vec == seq2Vec);
      int match(sum(matchesVec));
      matches.push_back(match);
      aLength.push_back(seq1.length());
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
                                             Named("strand1") = strand1,
                                             Named("strand2") = strand2,
                                             Named("scores") = scores,
                                             Named("matches") = matches,
                                             Named("aLength") = aLength);
  return outputList;
}
