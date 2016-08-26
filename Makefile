all: CImg.h a2.cpp
	g++ a2.cpp -o a2 -lX11 -lpthread -I. -Isiftpp -O3 siftpp/sift.cpp

clean:
	rm a2
	- rm ./a2-images/part2_images/seq2/*-warped.jpg
	- rm ./a2-images/part2_images/seq1/*-warped.jpg

run-part2: all
	- rm ./a2-images/part2_images/seq2/*-warped.jpg
	find ./a2-images/part2_images/seq2/ -type f | xargs ./a2 part2

run-part1: all
	- rm ./a2-images/part2_images/seq1/*-warped.jpg
	find ./a2-images/part2_images/seq1/ -type f | xargs ./a2 part2
