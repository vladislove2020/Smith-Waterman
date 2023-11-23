#ifndef alignment_h
#define alignment_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <list>
#include <ios>

struct Alignment {
public:
	int start1 = 0;
	int start2 = 0;
	std::list<bool> alignment;
};

Alignment Align(std::string seq1, std::string seq2,
	double gap_extend, double match, double mismatch);

Alignment Align_affine(std::string seq1, std::string seq2,
	double gap_open, double gap_extend, double match, double mismatch);

void Print_alignment(std::istream&, std::ostream&,
	double, double, double, double);

#endif