#include "SW.h"

int main() {
	//collect all parameters from stdin
	std::cout << "This is Smith-Waterman algorithm implementation" << std::endl;
	std::cout << "Let's try to align smth" << std::endl;
	std::cout << "Provide alignment parameters. If none is provided, default value will be used" << std::endl << std::endl;

	std::cout << "Points per match (floating point number, default 2):" << std::endl;
	double match = 2;
	if (!(std::cin.peek() == '\n')) {
		std::cin >> match;
	}

	std::cout << "Points per mismatch (floating point number, default -1):" << std::endl;
	double mismatch = -1;
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	if (!(std::cin.peek() == '\n')) {
		std::cin >> mismatch;
	}

	std::cout << "Use affine gap penalties (y/n, default n):" << std::endl;
	char answer = 'n';
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	if (!(std::cin.peek() == '\n')) {
		std::cin >> answer;
	}

	//make alignment with either linear or affine gap penalties
	if (answer == 'y') {
		std::cout << "Points per gap extension (floating point number, default -0.5):" << std::endl;
		double gap_extend = -0.5;
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		if (!(std::cin.peek() == '\n')) {
			std::cin >> gap_extend;
		}

		std::cout << "Points per gap opening (floating point number, default -5):" << std::endl;
		double gap_open = -5;
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		if (!(std::cin.peek() == '\n')) {
			std::cin >> gap_open;
		}

		std::cout << std::endl;
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		Print_alignment(std::cin, std::cout, gap_extend, match, mismatch, gap_open);
	}
	else {
		std::cout << "Points per gap extension (floating point number, default -2):" << std::endl;
		double gap_extend = -2;
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		if (!(std::cin.peek() == '\n')) {
			std::cin >> gap_extend;
		}

		std::cout << std::endl;
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		Print_alignment(std::cin, std::cout, gap_extend, match, mismatch, NULL);
	}

	return 0;
}