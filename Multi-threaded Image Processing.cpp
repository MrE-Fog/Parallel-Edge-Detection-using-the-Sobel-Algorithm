// readWrite-bmp.cc
//
// extracts pixel data from user-specified .bmp file
// inserts data back into new .bmp file
//
// gw

// uncomment for MSVS#
#include "stdafx.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <omp.h>
#include <stdio.h>
#include <thread>



using namespace std;

#pragma pack(1)
typedef struct {
	char id[2];
	int file_size;
	int reserved;
	int offset;
}  header_type;

#pragma pack(1)
typedef struct {
	int header_size;
	int width;
	int height;
	unsigned short int color_planes;
	unsigned short int color_depth;
	unsigned int compression;
	int image_size;
	int xresolution;
	int yresolution;
	int num_colors;
	int num_important_colors;
} information_type;


//sobel algorithm declaration
vector <vector <int>> sobel(int width, int height, vector <vector <int> >, vector <vector <int> >);
void sobel_thread(int width, int start_height, int stop_height, vector <vector <int> > newImageData, vector <vector <int> > oldData);

vector <vector <int>> x;
clock_t start, finish;

//vector <vector <int>> newData;

//vector<vector<int>> newDatar(416, vector<int>(2, 0));
vector<vector<int> > newData;



int main(int argc, char* argv[])
{
	header_type header;
	information_type information;
	string imageFileName, newImageFileName;
	unsigned char tempData[3];
	int row, col, row_bytes, padding;
	vector <vector <int> > data, newData, something;

	// prepare files
	cout << "Original imagefile? ";
	cin >> imageFileName;
	ifstream imageFile;
	imageFile.open(imageFileName.c_str(), ios::binary);
	if (!imageFile) {
		cerr << "file not found" << endl;
		exit(-1);
	}
	cout << "New imagefile name? ";
	cin >> newImageFileName;
	ofstream newImageFile;
	newImageFile.open(newImageFileName.c_str(), ios::binary);

	// read file header
	imageFile.read((char *)&header, sizeof(header_type));
	if (header.id[0] != 'B' || header.id[1] != 'M') {
		cerr << "Does not appear to be a .bmp file.  Goodbye." << endl;
		exit(-1);
	}
	// read/compute image information
	imageFile.read((char *)&information, sizeof(information_type));
	row_bytes = information.width * 3;
	padding = row_bytes % 4;
	if (padding)
		padding = 4 - padding;

	// extract image data, initialize vectors
	for (row = 0; row < information.height; row++) {
		data.push_back(vector <int>());
		for (col = 0; col < information.width; col++) {
			imageFile.read((char *)tempData, 3 * sizeof(unsigned char));
			data[row].push_back((int)tempData[0]);
		}
		if (padding)
			imageFile.read((char *)tempData, padding * sizeof(unsigned char));
	}
	cout << imageFileName << ": " << information.width << " x " << information.height << endl;

	// this loop shows how to simply recreate the original Black-and-White image
	for (row = 0; row < information.height; row++) {
		newData.push_back(vector <int>());
		for (col = 0; col < information.width; col++) {
			newData[row].push_back(data[row][col]);
		}
	}



	//divide 2d matrix for threading
	int total_rows = information.height;
	int rows = double(total_rows % 4);
	int rows_per_thread = floor(total_rows / 4);
	//Calculate start and stop rows for threads
	int row_t1 = rows_per_thread;
	int row_t2 = rows_per_thread * 2;
	int row_t3 = rows_per_thread * 3;
	int row_t4 = rows_per_thread * 4 + rows;

	//thread thread_1 (sobel_thread, information.width, 0, row_t1- 1, newData, data);
	//thread thread_2(sobel_thread, information.width, row_t1, row_t2 - 1, newData, data);
	//thread thread_3(sobel_thread, information.width, row_t2, row_t3 - 1, newData, data);
//	thread thread_4(sobel_thread, information.width, row_t3, row_t4, newData, data);

	//thread_1.join();
	//thread_2.join();
	//thread_3.join();
	//thread_4.join();

	//start clock for timing
	start = clock();

	//call sobel function and store matrix data into variable newData
	//calls perform sobel function on different rows - currently last one run is the 
	//only one that is written to the file
	//sobel_thread(information.width,0, row_t1 - 1, newData, data);
	sobel_thread(information.width, row_t1, row_t2 - 1, newData, data);
	//sobel_thread(information.width, row_t2, row_t3 - 1, newData, data);
	//sobel_thread(information.width, row_t3, row_t4, newData, data);
	
	//stop clock
	finish = clock();

	// write header to new image file
	newImageFile.write((char *)&header, sizeof(header_type));
	newImageFile.write((char *)&information, sizeof(information_type));

	// write new image data to new image file
	for (row = 0; row < information.height; row++) {
		for (col = 0; col < information.width; col++) {
			tempData[0] = (unsigned char)x[row][col];
			tempData[1] = (unsigned char)x[row][col];
			tempData[2] = (unsigned char)x[row][col];
			newImageFile.write((char *)tempData, 3 * sizeof(unsigned char));
		}
		if (padding) {
			tempData[0] = 0;
			tempData[1] = 0;
			tempData[2] = 0;
			newImageFile.write((char *)tempData, padding * sizeof(unsigned char));
		}
	}
	cout << newImageFileName << " done." << endl;
	imageFile.close();
	newImageFile.close();

	return 0;
}

vector <vector <int> > sobel(int width, int height, vector <vector <int> > newImageData, vector <vector <int> > oldData) {
	//filter
	int dx[3][3] = { { 1,0,-1 },{ 2,0,-2 },{ 1,0,-1 } };
	int SUM;

	for (int y = 0; y < height - 2; y++) {
		for (int x = 0; x < width - 2; x++) {
			if (y == 0 || y >= height - 1 || x == 0 || x >= width - 1) {
				newImageData[y][x] = oldData[y][x];
				continue;
			}

			int sumX = 0;
			int sumY = 0;

			for (int i = -1; i < 2; i++) {
				for (int j = -1; j < 2; j++) {
					sumX = sumX + dx[j + 1][i + 1] * (int)oldData[y + j][x + i];
					sumY = sumY + dx[j + 1][i + 1] * (int)oldData[y + j][x + i];
				}
			}
			/*Edge strength*/
			SUM = sqrt(pow((double)sumX, 2) + pow((double)sumY, 2));
			//threshold
			if (SUM > 255) SUM = 255;
			if (SUM < 20) SUM = 0;

			newImageData[y][x] = SUM;
		}
	}
	return newImageData;
}

void sobel_thread(int width, int start_height, int stop_height, vector <vector <int> > newImageData, vector <vector <int> > oldData) {
	//filter
	int dx[3][3] = { { 1,0,-1 },{ 2,0,-2 },{ 1,0,-1 } };
	int SUM;

	for (int y = start_height; y < stop_height - 2; y++) {
		for (int x = 0; x < width - 2; x++) {
			if (y == 0 || y >= stop_height - 1 || x == 0 || x >= width - 1) {
				newImageData[y][x] = oldData[y][x];
				continue;
			}

			int sumX = 0;
			int sumY = 0;

			for (int i = -1; i < 2; i++) {
				for (int j = -1; j < 2; j++) {
					sumX = sumX + dx[j + 1][i + 1] * (int)oldData[y + j][x + i];
					sumY = sumY + dx[j + 1][i + 1] * (int)oldData[y + j][x + i];
				}
			}
			/*Edge strength*/
			SUM = sqrt(pow((double)sumX, 2) + pow((double)sumY, 2));
			//threshold
			if (SUM > 255) SUM = 255;
			if (SUM < 20) SUM = 0;

			newImageData[y][x] = SUM;
		}
	}
	/*
	for (int row = start_height; row < stop_height; row++) {
		newData.push_back(vector <int>());
		for (int col = 0; col < width; col++) {
			newData[row].push_back(newImageData[row][col]);
		}
	}
	*/
	x = newImageData;

}




