#include "SW.h"

Alignment Align(std::string seq1, std::string seq2,
	double gap_extend = -2, double match = 2, double mismatch = -1) {
	//seq1 - horizontal, seq2 - vertical
	std::vector<std::vector<double>> matrix;
	matrix.resize(seq1.size() + 1, std::vector<double>(seq2.size() + 1));

	//first col and first row are zeros
	for (int i = 0; i < matrix.size(); i++) {
		matrix[i][0] = 0;
	}
	for (int j = 0; j < matrix[0].size(); j++) {
		matrix[0][j] = 0;
	}

	//fill matrix with weights, remember (i,j) for the biggest weight
	int i_max = 0;
	int j_max = 0;
	double weight_max = 0;
	for (int i = 1; i < matrix.size(); i++) {
		for (int j = 1; j < matrix[0].size(); j++) {
			std::vector<double> variants = { matrix[i - 1][j - 1] + match * (seq1[i - 1] == seq2[j - 1]) + mismatch * (seq1[i - 1] != seq2[j - 1]),
				matrix[i - 1][j] + gap_extend,
				matrix[i][j - 1] + gap_extend,
				0 };
			matrix[i][j] = *std::max_element(variants.begin(), variants.end());
			if (matrix[i][j] > weight_max) {
				i_max = i;
				j_max = j;
				weight_max = matrix[i][j];
			}
		}
	}

	//initialise Alignment object
	Alignment A;

	//backward propagation through the matrix from (i_max,j_max)
	int i = i_max;
	int j = j_max;
	while (true) {
		//if we reached start of any sequence
		if ((i == 0) || (j == 0)) {
			break;
		}

		//or if we reached non-profit position
		if (matrix[i][j] <= 0) {
			break;
		}

		//we must choose previous position with the biggest weight from possible variants
		char opt_from = '-';
		double opt_weight = -std::numeric_limits<double>::max();
		if ((matrix[i][j] == matrix[i - 1][j - 1] + match) || (matrix[i][j] == matrix[i - 1][j - 1] + mismatch)) {
			if (matrix[i - 1][j - 1] > opt_weight) {
				opt_from = 'd';
				opt_weight = matrix[i - 1][j - 1];
			}
		}
		if (matrix[i][j] == matrix[i - 1][j] + gap_extend) {
			if (matrix[i - 1][j] > opt_weight) {
				opt_from = 'v';
				opt_weight = matrix[i - 1][j];
			}
		}
		if (matrix[i][j] == matrix[i][j - 1] + gap_extend) {
			if (matrix[i][j - 1] > opt_weight) {
				opt_from = 'h';
				opt_weight = matrix[i][j - 1];
			}
		}

		//add suitable variant into A.alignment relying on the value of opt_prev
		//we use list to efficiently add elements to the front of our alignment
		switch (opt_from) {
		case 'd':
			A.alignment.push_front(true);
			A.alignment.push_front(false);
			i--;
			j--;
			break;
		case 'v':
			A.alignment.push_front(false);
			i--;
			break;
		case 'h':
			A.alignment.push_front(true);
			j--;
			break;
		}
	}

	//(start1,start2) is the first position of local alignment, and A.alignment is its binary representation
	A.start1 = i;
	A.start2 = j;

	return A;
}


Alignment Align_affine(std::string seq1, std::string seq2,
	double gap_open = -5, double gap_extend = -0.5, double match = 2, double mismatch = -1) {
	//seq1 - horizontal, seq2 - vertical
	std::vector<std::vector<double>> matrix;
	matrix.resize(seq1.size() + 1, std::vector<double>(seq2.size() + 1));
	std::vector<std::vector<double>> gap_x;
	gap_x.resize(seq1.size() + 1, std::vector<double>(seq2.size() + 1));
	std::vector<std::vector<double>> gap_y;
	gap_y.resize(seq1.size() + 1, std::vector<double>(seq2.size() + 1));

	//first col and first row are zeros
	for (int i = 0; i < matrix.size(); i++) {
		matrix[i][0] = 0;
		gap_x[i][0] = 0;
		gap_y[i][0] = 0;
	}
	for (int j = 0; j < matrix[0].size(); j++) {
		matrix[0][j] = 0;
		gap_x[0][j] = 0;
		gap_y[0][j] = 0;
	}

	//fill matrix with weights, remember (i,j) for the biggest weight in matrix
	//think that we cannot start with gap or finish with gap (because gaps are penalted)
	int i_max = 0;
	int j_max = 0;
	double weight_max = 0;
	for (int i = 1; i < matrix.size(); i++) {
		for (int j = 1; j < matrix[0].size(); j++) {
			//deletion in seq2 way
			std::vector<double> variants = { matrix[i - 1][j] + gap_open,
				gap_x[i - 1][j] + gap_extend };
			gap_x[i][j] = *std::max_element(variants.begin(), variants.end());

			//deletion in seq1 way
			variants = { matrix[i][j - 1] + gap_open,
				gap_y[i][j - 1] + gap_extend };
			gap_y[i][j] = *std::max_element(variants.begin(), variants.end());

			//match/mismatch way
			variants = { matrix[i - 1][j - 1] + match * (seq1[i - 1] == seq2[j - 1]) + mismatch * (seq1[i - 1] != seq2[j - 1]),
				gap_x[i][j],
				gap_y[i][j],
				0 };
			matrix[i][j] = *std::max_element(variants.begin(), variants.end());
			if (matrix[i][j] > weight_max) {
				i_max = i;
				j_max = j;
				weight_max = matrix[i][j];
			}
		}
	}

	//initialise Alignment object
	Alignment A;

	//backward propagation through the matrix from (i_max,j_max)
	int i = i_max;
	int j = j_max;
	while (true) {
		//if we reached start of any sequence
		if ((i == 0) || (j == 0)) {
			break;
		}

		//or if we reached non-profit position
		if (matrix[i][j] <= 0) {
			break;
		}

		//we must choose previous position with the biggest weight from possible variants
		//now we move through all three 2d vectors, so we must set a variable for the current condition (pointer is better to save memory)
		char opt_from = '-';
		double opt_weight = -std::numeric_limits<double>::max();
		std::vector<std::vector<double>>* current = &matrix;

		if (*current == matrix) {
			if ((matrix[i][j] == matrix[i - 1][j - 1] + match) || (matrix[i][j] == matrix[i - 1][j - 1] + mismatch)) {
				if (matrix[i - 1][j - 1] > opt_weight) {
					opt_from = 'd';
					opt_weight = matrix[i - 1][j - 1];
				}
			}
			if (matrix[i][j] == gap_x[i][j]) {
				if (gap_x[i][j] > opt_weight) {
					opt_from = 'v';
					opt_weight = gap_x[i][j];
					current = &gap_x;
				}
			}
			if (matrix[i][j] == gap_y[i][j]) {
				if (gap_y[i][j] > opt_weight) {
					opt_from = 'h';
					opt_weight = gap_y[i][j];
					current = &gap_y;
				}
			}

			//add suitable variant into A.alignment relying on the value of opt_prev
			//we use list to efficiently add elements to the front of our alignment
			switch (opt_from) {
			case 'd':
				A.alignment.push_front(true);
				A.alignment.push_front(false);
				i--;
				j--;
				break;
			case 'v':
				A.alignment.push_front(false);
				i--;
				break;
			case 'h':
				A.alignment.push_front(true);
				j--;
				break;
			}
		}

		else if (*current == gap_x) {
			if (gap_x[i][j] == matrix[i - 1][j] + gap_open) {
				if (matrix[i - 1][j] > opt_weight) {
					opt_from = 'm';
					opt_weight = matrix[i - 1][j];
					current = &gap_x;
				}
			}
			if (gap_x[i][j] == gap_x[i - 1][j] + gap_extend) {
				if (gap_x[i - 1][j] > opt_weight) {
					opt_from = 'x';
					opt_weight = gap_x[i - 1][j];
				}
			}

			switch (opt_from) {
			case 'm':
				A.alignment.push_front(true);
				A.alignment.push_front(false);
				i--;
				j--;
				break;
			case 'x':
				A.alignment.push_front(false);
				i--;
				break;
			}
		}

		else if (*current == gap_y) {
			if (gap_y[i][j] == matrix[i][j - 1] + gap_open) {
				if (matrix[i][j - 1] > opt_weight) {
					opt_from = 'm';
					opt_weight = matrix[i][j - 1];
					current = &matrix;
				}
			}
			if (gap_y[i][j] == gap_y[i][j - 1] + gap_extend) {
				if (gap_y[i][j - 1] > opt_weight) {
					opt_from = 'y';
					opt_weight = gap_y[i][j - 1];
				}
			}

			switch (opt_from) {
			case 'm':
				A.alignment.push_front(true);
				A.alignment.push_front(false);
				i--;
				j--;
				break;
			case 'y':
				A.alignment.push_front(true);
				j--;
				break;
			}
		}
	}

	//(start1,start2) is the first position of local alignment, and A.alignment is its binary representation
	A.start1 = i;
	A.start2 = j;

	return A;
}