/*
 * SiftMatcher.cpp
 *
 *  Created on: 19 Feb 2016
 *      Author: sagar
 */

#include <vector>
#include "CImg.h"
#include <Sift.h>
#include <queue>
#include<dirent.h>

using namespace std;
using std::priority_queue;

class SiftMatcher {
public:
	static const float MATCH_THRESHOLD = 0.63;

	class SiftMatch {
	public:
	  SiftMatch() {
	    weight = 0;
	    tries = 0;
	    toDelete = false;
	  }
		int src_x, src_y, dest_x, dest_y;
		double compare_val;
	  int weight;
	  int tries;
	  bool toDelete;
		vector<float> sift_vector;
		vector<int> summary_vector;
		string file_name;
		void display();
	  bool operator > (const SiftMatcher::SiftMatch  &s) const
	  {
	    if (s.tries == 0)
	      return false;
	    if (tries == 0)
	      return true;
	    return weight/tries > s.weight/s.tries;
	  }
	};



	struct SiftMatchLesserCompare {
		bool operator()(const SiftMatcher::SiftMatch &t1,
				const SiftMatcher::SiftMatch &t2) const {
			return t1.compare_val < t2.compare_val;
		}
	};

	struct SiftMatchGreaterCompare {
		bool operator()(const SiftMatcher::SiftMatch &t1,
				const SiftMatcher::SiftMatch &t2) const {
			return t1.compare_val > t2.compare_val;
		}
	};
  	static vector<SiftMatcher::SiftMatch> match_descriptors(
			const vector<SiftDescriptor> &src,
			const vector<SiftDescriptor> &dest);

	static vector<SiftMatcher::SiftMatch> match_descriptors(
			const vector<SiftDescriptor> &src,
			const vector<SiftDescriptor> &dest,
			double THRESHOLD,double mat[3][3]);

	static CImg<double> get_sift_matched_result(const CImg<double> &input_image,
			const CImg<double> &match_image, bool useQuantized);

	static vector<string> get_ordered_sift_matches(const string inputFile,
			const vector<string> query_images, bool useQuantized);

	static void print_ordered_sift_matches(const string inputFile,
			const vector<string> query_images, int print_limit, bool task_three,
			bool useQuantized);

	static void match_directory(const string inputDir, bool useQuantized,
			string compare_image);
};

class QuantizedSiftMatcher {

public:
	static const int sift_size = 128;

	static const double quant_threshold = 10.0;

	static const double sample_count = 40;

	static vector<double> get_sampled_vectors(int count);

	static vector<SiftMatcher::SiftMatch> get_summary_vector(
			const vector<SiftDescriptor> &src, vector<double> sampled_vectors);

	static vector<SiftMatcher::SiftMatch> match_descriptors(
			const vector<SiftDescriptor> &src,
			const vector<SiftDescriptor> &dest);
};

/**
 * returns vector of k*128 sampled values of uniform distribution of [0,1], equivalent to k vectors of length 128
 */
inline vector<double> QuantizedSiftMatcher::get_sampled_vectors(int k) {

	vector<double> samped_vector;

	for (int i = 0; i < QuantizedSiftMatcher::sift_size * k; i++) {
		double r = ((double) rand() / (RAND_MAX));
		samped_vector.push_back(r);
	}

	return samped_vector;
}

/**
 * calculates summary vectors for sift vectors image by dot product and applying quantization threshold
 */
inline vector<SiftMatcher::SiftMatch> QuantizedSiftMatcher::get_summary_vector(
		const vector<SiftDescriptor> &src, vector<double> sampled_vectors) {

	vector<SiftMatcher::SiftMatch> summary_vector;

	//for each descriptor in sift match
	for (int i = 0; i < src.size(); i++) {

		vector<int> summary;
		//for each sample vector
		for (int j = 0;
				j < sampled_vectors.size() / QuantizedSiftMatcher::sift_size;
				j++) {

			double dist = 0.0;
			//multiply by each vector
			for (int k = 0; k < QuantizedSiftMatcher::sift_size; k++) {
				int pos = j * QuantizedSiftMatcher::sift_size + k;
				dist += src[i].descriptor[k] * sampled_vectors[pos];
			}
			summary.push_back(dist);
		}
		SiftMatcher::SiftMatch sMatch;
		sMatch.summary_vector = summary;
		summary_vector.push_back(sMatch);
	}

	return summary_vector;
}

inline vector<SiftMatcher::SiftMatch> QuantizedSiftMatcher::match_descriptors(
		const vector<SiftDescriptor>& src, const vector<SiftDescriptor>& dest) {

	cout
			<< "Matching the sift descriptors for source and query image using quantized descriptor"
			<< endl;

	vector<SiftMatcher::SiftMatch> matches;

	//crude implementation

	//iterate over all the descriptor values of destination for each value of source to find minimum. preserve best and second best

	//Generate sampled projection vectors
	vector<double> samped_vector = QuantizedSiftMatcher::get_sampled_vectors(
			QuantizedSiftMatcher::sample_count);

	//calculate summary vector for first image
	vector<SiftMatcher::SiftMatch> src_summary_vector =
			QuantizedSiftMatcher::get_summary_vector(src, samped_vector);

	//calculate summary vector for second image
	vector<SiftMatcher::SiftMatch> dest_summary_vector =
			QuantizedSiftMatcher::get_summary_vector(dest, samped_vector);

	//compare both vectors and get best matches
	for (int i = 0; i < src_summary_vector.size(); i++) {

		vector<int> src_summary = src_summary_vector[i].summary_vector;

		//priority_queue
		priority_queue<SiftMatcher::SiftMatch, vector<SiftMatcher::SiftMatch>,
				SiftMatcher::SiftMatchGreaterCompare> que;

		for (int j = 0; j < dest_summary_vector.size(); j++) {

			vector<int> dest_summary = dest_summary_vector[j].summary_vector;

			double diff = 0.0;
			//compare it with source vector
			for (int k = 0; k < src_summary.size(); k++) {
				diff += pow(src_summary[k] - dest_summary[k],2);
			}

			//add results to priority queue
			SiftMatcher::SiftMatch match;
			match.compare_val = diff;
			match.src_x = j;
			match.src_y = i;

			que.push(match);
		}

		//take first 3 matches
		SiftMatcher::SiftMatch best_match = que.top();
		que.pop();
		SiftMatcher::SiftMatch second_match = que.top();
		que.pop();
		SiftMatcher::SiftMatch last_match = que.top();

		vector<SiftDescriptor> best_matches;
		best_matches.push_back(dest[best_match.src_x]);
		best_matches.push_back(dest[second_match.src_x]);
		best_matches.push_back(dest[last_match.src_x]);

		vector<SiftDescriptor> temp;
		temp.push_back(src[best_match.src_y]);

		vector<SiftMatcher::SiftMatch> mat = SiftMatcher::match_descriptors(
				temp, best_matches);

		if (mat.size() > 0) {
			matches.push_back(mat[0]);
		}
	}

	cout << "Match complete, total points matched : " << matches.size()
			<< " current threshold is: " << SiftMatcher::MATCH_THRESHOLD
			<< endl;

	return matches;
}

void SiftMatcher::SiftMatch::display() {
	std::cout << "src X: " << this->src_x << " src Y: " << this->src_y
			<< " dest X: " << this->dest_x << " dest Y: " << this->dest_y
			<< endl;
}

inline vector<SiftMatcher::SiftMatch> SiftMatcher::match_descriptors(const vector<SiftDescriptor>& src, const vector<SiftDescriptor>& dest) {
  return SiftMatcher::match_descriptors(src,dest,0,NULL);
}
inline vector<SiftMatcher::SiftMatch> SiftMatcher::match_descriptors(const vector<SiftDescriptor>& src, const vector<SiftDescriptor>& dest,double THRESHOLD,double mat[3][3]) {

	/*cout << "Matching the sift descriptors for source and query image" << endl;*/

	vector<SiftMatcher::SiftMatch> matches;
	THRESHOLD = THRESHOLD > 0 ? THRESHOLD : MATCH_THRESHOLD;
	//crude implementation

	//iterate over all the descriptor values of destination for each value of source to find minimum. preserve best and second best

	//for each key-point from source
	for (int i = 0; i < src.size(); i++) {

		SiftMatcher::SiftMatch best, second;
		long bestDist = std::numeric_limits<long>::max(), secondDist =
				std::numeric_limits<long>::max();

		best.src_x = src[i].col;
		best.src_y = src[i].row;
		second.src_x = src[i].col;
		second.src_y = src[i].row;

		//for each key-point in destination
		for (int j = 0; j < dest.size(); j++) {

			long dist = 0;

			for (int k = 0; k < 128; k++) {
				int diff = src[i].descriptor[k] - dest[j].descriptor[k];
				dist += pow(diff, 2);
			}

			//update best match
			if (dist < bestDist) {
				second = best;
				secondDist = bestDist;
				bestDist = dist;
				best.dest_x = dest[j].col;
				best.dest_y = dest[j].row;
			} else if (dist < secondDist) {
				//update second best match
				secondDist = dist;
				second.dest_x = dest[j].col;
				second.dest_y = dest[j].row;
			}
		}

		//if ratio of best to second best is less than threshold, accept that point
		if ((float) bestDist / (float) secondDist
				< THRESHOLD) {
			matches.push_back(best);
		}
	}
	
	
/*	Otimized sift matcher (Commented out due to performance reason)
	Optimized matching strategy for sift descriptors for obtaining better matches. Uses approximate homographic matrix to evaluate sift key-point match

 	//for each key-point from source
	for (int i = 0; i < src.size(); i++) {

		//priority_queue
		priority_queue<SiftMatcher::SiftMatch, vector<SiftMatcher::SiftMatch>,

		SiftMatcher::SiftMatchGreaterCompare> que;

		//for each key-point in destination
		for (int j = 0; j < dest.size(); j++) {
			double dist = 0;
			for (int k = 0; k < 128; k++) {
				int diff = src[i].descriptor[k] - dest[j].descriptor[k];
				dist += pow(diff, 2);
			}

			SiftMatcher::SiftMatch best;
			best.src_x = src[i].col;
			best.src_y = src[i].row;
			best.dest_x = dest[j].col;
			best.dest_y = dest[j].row;
			best.compare_val = dist;
			que.push(best);
		}

		if (mat != NULL) {
			SiftMatcher::SiftMatch best;

			bool bestSet = false;

			for (int m = 0; m < 10; m++) {
				SiftMatcher::SiftMatch match = que.top();
				que.pop();

				SiftMatcher::SiftMatch match2 = que.top();

				//if ratio of best to second best is less than threshold, accept that point
				if (match.compare_val / match2.compare_val < 0.7) {
					if (!bestSet) {
						best = match;
						bestSet = true;
					}

					if (is_matches_homographic(match.src_x, match.src_x,
					match.dest_x, match.dest_y, mat)) {
						matches.push_back(match);
						bestSet = false;
						break;
					}
				}
			}
			if (bestSet) {
				matches.push_back(best);
			}
		} else {
			SiftMatcher::SiftMatch match = que.top();
			que.pop();
			SiftMatcher::SiftMatch match2 = que.top();
			//if ratio of best to second best is less than threshold, accept that point
			if ((float) match.compare_val / (float) match2.compare_val
					< THRESHOLD) {
				matches.push_back(match);
			}
		}
	}*/

	/*cout << "Match complete, total points matched : " << matches.size()
	 << " current threshold is: " << SiftMatcher::MATCH_THRESHOLD
	 << endl;
	 */
	return matches;
}

/* Commented out for performance reason
//returns true if calculated match follows homographics obtained from rainsac algorithm, used by optimized sift matcher
bool is_matches_homographic(double srcX, double srcY, double dstX, double dstY,
		double mat[3][3]) {
	double x = srcX, y = dstX;
	double w = mat[2][0] * x + mat[2][1] * y + 1;
	double x1 = (mat[0][0] * x + mat[0][1] * y + mat[0][2]) / w;
	double y1 = (mat[1][0] * x + mat[1][1] * y + mat[1][2]) / w;

	if (abs(dstX - x1) < 30 && abs(dstY - y1) < 30) {
		return true;
	} else {
		return false;
	}
}
*/

inline CImg<double> SiftMatcher::get_sift_matched_result(
		const CImg<double>& input_image, const CImg<double>& match_image,
		bool useQuantized) {

	cout << "Creating merged image from input and query image" << endl;

	//create single image by appending two images
	CImg<double> merge_img(input_image);
	merge_img.append(match_image);

	// convert input image to grayscale
	CImg<double> input_gray = input_image.get_RGBtoHSI().get_channel(2);
	vector<SiftDescriptor> input_descriptor = Sift::compute_sift(input_gray);

	// convert match image to grayscale
	CImg<double> match_gray = match_image.get_RGBtoHSI().get_channel(2);
	vector<SiftDescriptor> match_descriptor = Sift::compute_sift(match_gray);

	//match the key-points using the above descriptors
	vector<SiftMatcher::SiftMatch> matches =
			useQuantized ?
					QuantizedSiftMatcher::match_descriptors(input_descriptor,
							match_descriptor) :
					SiftMatcher::match_descriptors(input_descriptor,
							match_descriptor);

	//now we have matches with us, lets draw line between matched points
	for (int i = 0; i < matches.size(); i++) {
		int sourceX = matches[i].src_x;
		int sourceY = matches[i].src_y;
		int destX = matches[i].dest_x + input_image.width();
		int destY = matches[i].dest_y;

		//yellow color
		const unsigned char color[] = { 255, 255, 0 };

		merge_img.draw_line(sourceX, sourceY, destX, destY, color, 1);
	}

	cout << "Drawing lines based on matched points" << endl;

	return merge_img;
}

void get_files(string dir, vector<string> &files) {
	DIR *dp = opendir(dir.c_str());
	struct dirent *directory_p;
	while ((directory_p = readdir(dp)) != NULL) {
		if (string(directory_p->d_name) != "."
				&& string(directory_p->d_name) != "..") {
			files.push_back(dir + string(directory_p->d_name));
		}
	}
	closedir(dp);
}

inline vector<string> SiftMatcher::get_ordered_sift_matches(
		const string inputFile, const vector<string> query_images,
		bool useQuantized) {

	cout << "matching query image with each from the argument the list" << endl;
	CImg<double> input_image(inputFile.c_str());

	// convert input image to grayscale
	CImg<double> input_gray = input_image.get_RGBtoHSI().get_channel(2);
	vector<SiftDescriptor> input_descriptor = Sift::compute_sift(input_gray);

	priority_queue<SiftMatcher::SiftMatch, vector<SiftMatcher::SiftMatch>,
			SiftMatcher::SiftMatchLesserCompare> que;

	for (int i = 0; i < query_images.size(); i++) {
		string matchFile = query_images[i];
		CImg<double> match_image(matchFile.c_str());

		cout << "Matching image: " << matchFile << endl;
		// convert match image to grayscale
		CImg<double> match_gray = match_image.get_RGBtoHSI().get_channel(2);
		vector<SiftDescriptor> match_descriptor = Sift::compute_sift(
				match_gray);

		//match the key-points using the above descriptors
		vector<SiftMatcher::SiftMatch> matches =
				useQuantized ?
						QuantizedSiftMatcher::match_descriptors(
								input_descriptor, match_descriptor) :
						SiftMatcher::match_descriptors(input_descriptor,
								match_descriptor);

		//add matched feature size to priority queue
		SiftMatcher::SiftMatch match;
		match.compare_val = matches.size();
		match.file_name = matchFile;
		que.push(match);
	}

	cout
			<< "Part 1, task 2. Printing matched images by size of matched features"
			<< endl;

	vector<string> ordered_matches;

	while (!que.empty()) {
		SiftMatcher::SiftMatch match = que.top();
		que.pop();
		ordered_matches.push_back(match.file_name);
	}

	return ordered_matches;
}

inline void SiftMatcher::print_ordered_sift_matches(const string inputFile,
		const vector<string> query_images, int print_limit, bool task_three,
		bool useQuantized) {

	vector<string> matchedFiles = SiftMatcher::get_ordered_sift_matches(
			inputFile, query_images, useQuantized);

	print_limit = print_limit < 0 ? matchedFiles.size() : print_limit;

	for (int i = 0; i < print_limit; i++) {
		cout << matchedFiles[i] << endl;
	}

	if (task_three) {

		int end = inputFile.find_last_of("_");
		int start = inputFile.find_last_of("//") + 1;

		string src_attraction = inputFile.substr(start, end - start);
		int correct = 0;
		//print accuracy
		for (int i = 0; i < 10; i++) {

			end = matchedFiles[i].find_last_of("_");
			start = matchedFiles[i].find_last_of("//") + 1;

			string detected_attraction = matchedFiles[i].substr(start,
					end - start);

			if (detected_attraction == src_attraction) {
				correct++;
			}
		}

		cout << "Precision of system: " << endl << "Correct detection: "
				<< correct * 10 << " \%" << endl;
	}
}

inline void SiftMatcher::match_directory(const string inputDir,
		bool useQuantized, string compare_image) {

	vector<string> images;
	get_files(inputDir, images);

	if (compare_image.length() == 0) {

		srand(time(NULL));
		int query_image_id = rand() % images.size();
		compare_image = images[query_image_id];
	}

	cout << "Matching for image " << compare_image << endl;
	SiftMatcher::print_ordered_sift_matches(compare_image, images, 10, true,
			useQuantized);
}
