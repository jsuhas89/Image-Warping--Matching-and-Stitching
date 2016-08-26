// B657 assignment 2 skeleton code

//

// Compile with: "make"

//

// See assignment handout for command line and project specifications.

//Link to the header file

#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Sift.h>
#include <SiftMatcher.cpp>
#include "Transforms.cpp"

//Use the cimg namespace to access the functions easily

using namespace cimg_library;
using namespace std;

int main(int argc, char **argv) {

	try {

		if (argc < 2) {
			cout << "Insufficent number of arguments; correct usage:" << endl;
			cout << "    a2-p1 part_id ..." << endl;
			return -1;
		}

		string part = argv[1];
		string inputFile = argv[2];

		string targetFile;
		if (argc > 3)
			targetFile = argv[3];

		if (part == "part1.1") {

			CImg<double> input_image(inputFile.c_str());
			string matchFile = argv[3];

			cout << "Part I task 1 processing started" << endl;

			//Task 1 - Two images to match
			CImg<double> match_image(matchFile.c_str());

			//matched image with lines drwan on it
			CImg<double> merge_img = SiftMatcher::get_sift_matched_result(
					input_image, match_image, false);

			merge_img.get_normalize(0, 255).save("result_sift_match.jpg");
			cout
					<< "Part I task 1 complete, result saved in result_sift_match.jpg file"
					<< endl;

		} else if (part == "part1.2") {

			CImg<double> input_image(inputFile.c_str());
			vector<string> compare_images;

			//match query image with all other images in the arguments and print in descending order
			for (int i = 3; i < argc; i++) {
				compare_images.push_back(argv[i]);
			}

			SiftMatcher::print_ordered_sift_matches(inputFile, compare_images,
					compare_images.size(), false, false);

		} else if (part == "part1.3") {

			string compare_image = "";

			if (argc > 3) {
				compare_image = argv[3];
			}

			SiftMatcher::match_directory(inputFile, false, compare_image);

		} else if (part == "part1.4") {

			string compare_image = "";

			if (argc > 3) {
				compare_image = argv[3];
			}

			SiftMatcher::match_directory(inputFile, true, compare_image);

		} else if (part == "part2") {

			CImg<double> input_image(inputFile.c_str());

			if (argc == 3) {
				double arr[3][3] = { 0.907, 0.258, -182, -0.153, 1.44, 58,
						-0.000306, 0.000731, 1 };
				projective_transform(input_image, arr).save("res.png");
			} else {
				// Q2 of part 2
				int st = 3;
				while (st < argc) {
					string matchfile = argv[st];
					string dot = ".jpg";
					int matchpos = matchfile.find(dot);
					string out = matchfile.substr(0, matchpos) + "-warped.jpg";
					cout << " Out file name is " << out << endl;
					CImg<double> match_image(matchfile.c_str());
					produce_warped(input_image, match_image).save(out.c_str());
					st++;
				}
			}
		} else

			throw std::string("unknown part!");

	} catch (const string &err) {
		cerr << "Error: " << err << endl;
	}
}

