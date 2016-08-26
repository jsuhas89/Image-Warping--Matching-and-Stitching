// All Part2 methods and functions
//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <set>
#include <vector>

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;
// Some COnfig stuff
bool isCountBased = false ;

bool invert_matrix(double mat[3][3]) {
  // Get determinanant
  double determinant=0;
  for (int i = 0; i < 3; i++) {
    determinant += mat[0][i]*(mat[1][(i+1)%3]*mat[2][(i+2)%3] - mat[1][(i+2)%3]*mat[2][(i+1)%3]);
  }
  if (determinant == 0)
    return false;
  // Cofactor matrix
  double temp[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      temp[i][j] = ((mat[(i+1)%3][(j+1)%3] * mat[(i+2)%3][(j+2)%3]) - (mat[(i+1)%3][(j+2)%3]*mat[(i+2)%3][(j+1)%3]))/ determinant;
    }
  }

  for (int i = 0 ; i < 3; i++) {
    for (int j = 0 ; j < 3 ; j++) {
      mat[j][i] = temp[i][j];
    }
  }
  return true;
}

// 2.1
CImg<double> projective_transform(const CImg<double> &input,double mat[3][3]) {
  CImg<double> inp(input,"xyzc");
  bool isInverted = invert_matrix(mat);
  // Inverse warping
  // cout << "Projective transorm"  << "\n";
  if (isInverted) {
    // //Use inverse warping to get results
    // cout << "MAtrix "  << "\n";
    // for (int i = 0; i < 3; i++) {
    //   for (int j = 0; j < 3; j++) {
    // 	cout<< mat[i][j] << "\t";
    //   }
    //   cout<<"\n";
    // }
    // cout << "END "  << "\n";
    // Use inverse warping to get results
    cimg_forXY(inp,x,y) {
      double xyw[3] = {x,y,1};
      double res[3];
      for (int i = 0 ; i < 3 ; i++) {
	res[i] = 0;
	for (int j = 0; j < 3 ; j++) {
	  res[i]+=xyw[j] * mat[i][j];
	}
      }
      if (res[2] != 0) {
	double x1 = res[0]/res[2];
	double y1 = res[1]/res[2];
	if (x1 > 0 && x1 < input._width && y1 > 0 &&  y1 < input._height) {
	  inp(x,y,0) = input(x1,y1,0);
	  inp(x,y,1) = input(x1,y1,1);
	  inp(x,y,2) = input(x1,y1,2);
	}
      }
    }
  } else {
    cimg_forXY(inp,x,y) {
      double xyw[3] = {x,y,1};
      double res[3];
      for (int i = 0 ; i < 3 ; i++) {
	res[i] = 0;
	for (int j = 0; j < 3 ; j++) {
	  res[i]+=xyw[j] * mat[i][j];
	}
      }
      if (res[2] != 0) {
	double x1 = res[0]/res[2];
	double y1 = res[1]/res[2];
	if (x1 > 0 && x1 < input._width && y1 > 0 &&  y1 < input._height) {
	  inp(x1,y1,0) = input(x,y,0);
	  inp(x1,y1,1) = input(x,y,1);
	  inp(x1,y1,2) = input(x,y,2);
	}
      }
    }
    return inp;
  }
}

// Gets 4 random number from 0 to max-1
set<int> get_4_random_numbers(int max) {
  set<int> index_set;
  for (int i = 0 ; i < 4 ; i++) {
    while (1) {
      int ri = rand() %  max;
      // ri should not exist in the set already
      if (index_set.find(ri) == index_set.end()) {
	index_set.insert(ri);
	break;
      } else {
	// Try another random number
	continue;
      }
    }
  }
  return index_set;
}

// Gets SIFT matches for input image and match image. Passes optional hint matrix to the matcher
vector<SiftMatcher::SiftMatch> get_matches(const CImg<double> input_image, const CImg<double> match_image,double mat[3][3]=NULL) {
  // Lets do some sift matching here and populate matrix here
  // convert input image to grayscale
  CImg<double> input_gray = input_image.get_RGBtoHSI().get_channel(2);
  vector<SiftDescriptor> input_descriptor = Sift::compute_sift(input_gray);

  // convert match image to grayscale
  CImg<double> match_gray = match_image.get_RGBtoHSI().get_channel(2);
  vector<SiftDescriptor> match_descriptor = Sift::compute_sift(match_gray);
  vector<SiftMatcher::SiftMatch> matches;
  double threshold = SiftMatcher::MATCH_THRESHOLD - 0.1;
  // Need 10 points atleast
  // match the key-points using the above descriptors 0.53
  // Number of matches should be higher than 30 for this to have any chance of workign - Lets just move the thtreshold till it reaches that
  while (1) {
    matches = SiftMatcher::match_descriptors(input_descriptor, match_descriptor,threshold,mat);
    cout<< "No. of matches are " << matches.size() << "\n";
    if (matches.size() > 10 || threshold >= 1)
      break;
    else {
      threshold += 0.05;
    }
  }
  return matches;
}

// Takes in matches and set of 4 numbers , creates the matrix and returns the solved matrix as a Cimg
CImg<double> get_solved(const vector<SiftMatcher::SiftMatch> &matches,const set<int> index_set) {
  // 2 Cimages with 4p
  CImg<double> a(1,8);
  CImg<double> b(8,8);
  int st = 0;
  for (set<int>::iterator it=index_set.begin(); it!=index_set.end(); ++it) {
    SiftMatcher::SiftMatch m = matches[*it];
    double x = m.src_x , y = m.src_y , u = m.dest_x , v = m.dest_y;
    // cout << "Match " << x << ","<<y << "   " << u <<","<< v  << "\n";
    // double u = m.src_x , v = m.src_y , x = m.dest_x , y = m.dest_y;
    b(0,st) = x;b(1,st) = y;b(2,st) = 1;b(3,st) = 0;
    b(4,st) = 0, b(5,st) = 0, b(6,st) = -u * x , b(7,st) = -u * y;
    a(0,st) = u;
    st++;
    b(0,st) = 0; b(1,st) = 0; b(2,st) = 0; b(3,st) = x;
    b(4,st) = y; b(5,st) = 1; b(6,st) = -v * x; b(7,st) = -v * y;
    a(0,st) = v;
    st++;
  }
  a.solve(b);
  return CImg<double>(a);
}

// Takes [x,y] as arr reference and modifies it to x1 and y1
void transform_xy (double arr[2],double mat[3][3]) {
  double x = arr[0] , y = arr[1];
  double w = mat[2][0] * x + mat[2][1] * y + 1;
  double x1 = (mat[0][0] * x + mat[0][1] * y + mat[0][2])/w;
  double y1 = (mat[1][0] * x + mat[1][1] * y + mat[1][2])/w;
  arr[0] = x1;
  arr[1] = y1;
}


struct MatchStats {
  int totalMatches;
  double probability;
};

// Computes statistics between the matcher and the solution
MatchStats compute_stats(const vector<SiftMatcher::SiftMatch> &matches,const CImg<double> a) {
  MatchStats stats;
  double DELTA = 10 , PROB = 0.7;
  int match_count = 0;
  for (int it = 0; it != matches.size(); ++it) {
    SiftMatcher::SiftMatch m = matches[it];
    double x = m.src_x , y = m.src_y , u = m.dest_x , v = m.dest_y;
    double w = a(0,6) * x + a(0,7) * y + 1;
    double x1 = (a(0,0) * x + a(0,1) * y + a(0,2))/w;
    double y1 = (a(0,3) * x + a(0,4) * y + a(0,5))/w;
    if (x1 >= u - DELTA  && x1 <= u + DELTA && y1 >= v - DELTA && y1 <= v + DELTA)
      match_count++;
  }
  stats.totalMatches = match_count;
  stats.probability = match_count / matches.size();
  return stats;
}

// Modifies match vector and sets final Match according to match stats of last solution
void modify_matches_on_stats(vector<SiftMatcher::SiftMatch> &matches,const MatchStats &stat,
			     const set<int> &index_set,MatchStats &maxStat,const  CImg<double> &a
			     ,CImg<double> &finalMatch) {
  int match_count = stat.totalMatches;
  double probability = stat.probability;
  double maxProb = maxStat.probability;
  double PROB_LOWER_BOUND = 0.3;
  if (isCountBased) {
    // Increase weight accordingly
    for (set<int>::iterator it=index_set.begin(); it!=index_set.end(); ++it) {
      matches[*it].weight += match_count;
      matches[*it].tries++;
    }
  } else {
    // Is Probability based
    double prob = (double)match_count / matches.size();
    if (prob >= maxProb) {
      maxProb = prob;
      finalMatch.assign(a);
    } else if (prob < PROB_LOWER_BOUND) {
      // Reject these matches - ie delete them
      for (set<int>::iterator it=index_set.begin(); it!=index_set.end(); ++it) {
	matches[*it].toDelete = true;
      }
      // Delete those 4 points
      for (int i = 0; i < matches.size() && matches.size() > 4 ; i++) {
	if (matches[i].toDelete) {
	  matches.erase(matches.begin() + i);
	} else {
	  i++;
	}
      }
    }
  }
  maxStat.probability = maxProb;
}

// Copy values into matrix
void copy_into_matrix(CImg<double> finalMatch , double mat[3][3]) {
  mat[0][0] = finalMatch(0,0);mat[0][1] = finalMatch(0,1);mat[0][2] = finalMatch(0,2);
  mat[1][0] = finalMatch(0,3);mat[1][1] = finalMatch(0,4);mat[1][2] = finalMatch(0,5);
  mat[2][0] = finalMatch(0,6);mat[2][1] = finalMatch(0,7);mat[2][2] = 1;
}

// 2.2
// Computes homography matrix and assign it to passed matrix
bool compute_homography(const CImg<double> input_image, const CImg<double> match_image, double mat[3][3]) {
  int rounds = 100;
  MatchStats maxStat;
  CImg<double> finalMatch(1,8,1,1,1);
  double maxProb = 0;
  maxStat.totalMatches = 0;
  maxStat.probability = 0;
  vector<SiftMatcher::SiftMatch> matches;
  for (int r = 0 ; r < rounds ; r++) {
    // Try changing matches twicebn n8
    if (r % 50 == 0) {
      matches = get_matches(input_image,match_image,mat);
    }
    // Use the matches and calculate teta and keep removing outliers from the vector until
    // no more can be removed or only 4 matches are left
    // Now need to get 4 probable points from these matches using RANSAC
    if (matches.size() <= 4)
      break;
    // Choose 4 random points and use its range
    set<int> index_set = get_4_random_numbers(matches.size());
    int match_count = 0;
    CImg<double> a = get_solved(matches,index_set);
    MatchStats stat = compute_stats(matches,a);
    modify_matches_on_stats(matches,stat,index_set,maxStat,a,finalMatch);
    // Copy match data into matrix
    copy_into_matrix(finalMatch,mat);
  }

  if (isCountBased) {
    std::sort(matches.begin(), matches.end(),greater<SiftMatcher::SiftMatch>());
    // Take top 4 vectors and solve and use it to make matrix
    // 2 Cimages with 4p
    set<int> numset;
    // First 4
    numset.insert(0);numset.insert(2);numset.insert(2);numset.insert(3);
    finalMatch.assign(get_solved(matches,numset));
  } else {
    cout << "Got Final Probability " << maxStat.probability << "\n";
  }
  copy_into_matrix(finalMatch,mat);

  return true;
}

// 2.3 - Produces the warped image by getting homography matrix and then transforming it
CImg<double> produce_warped(CImg<double> img1, CImg<double> img2) {
  double mat[3][3];
  bool s = compute_homography(img2,img1,mat);
  for (int i = 0 ; i < 3;i++) {
    for (int j = 0 ; j < 3;j++) {
      cout << mat[i][j]  << "\t";
    }
    cout<< "\n";
  }
  return projective_transform(img2,mat);
}

