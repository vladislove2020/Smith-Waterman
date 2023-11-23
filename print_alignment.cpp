#include "SW.h"

void Print_alignment(std::istream& in, std::ostream& out,
	double gap_extend, double match, double mismatch, double gap_open = NULL) {
	//sequence1 from in (theoretically it could be modified to read from file and/or write to file)
	out << "Sequence 1:" << std::endl;
	std::string seq1;
	std::getline(in, seq1);

	//sequence2 from in
	out << "Sequence 2:" << std::endl;
	std::string seq2;
	std::getline(in, seq2);

	//create an alignment object from sequences
	Alignment A;
	if (gap_open) {
		A = Align_affine(seq1, seq2, gap_open, gap_extend, match, mismatch);
	}
	else {
		A = Align(seq1, seq2, gap_extend, match, mismatch);
	}

	//start from local alignment starting positions
	int i = A.start1;
	int j = A.start2;
	bool gap = false;

	//convert binary alignment representation into more human-readable double-string representation
	out << std::endl;
	//for (bool el : A.alignment) {
	for (auto iter = A.alignment.begin(); iter != A.alignment.end(); iter++) {
		if (!(*iter)) {
			out << seq1[i];
			i++;
			gap = false;
		}
		else {
			if (!gap) {
				gap = true;
			}
			else {
				out << "-";
			}
		}
	}

	out << std::endl;
	gap = false;

	for (auto iter = A.alignment.begin(); iter != A.alignment.end(); iter++) {
		if (*iter) {
			out << seq2[j];
			j++;
			gap = false;
		}
		else {
			if (!gap) {
				gap = true;
			}
			else {
				out << "-";
			}
		}
	}
}