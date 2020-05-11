
/* Data Mining implementation for Fuzzy c-Means Algorithm */
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <string>
#include <windows.h>
#include <algorithm>
#include <vector>
#include <cstdio>
#include <sstream>
#include <iomanip>
using namespace std;

#define MAX_DATA_POINTS 1000
#define MAX_CLUSTERS 100
#define MAX_DATA_DIMENSION 100

int num_data_points = 1000; // Number of data points 
int num_clusters = 5;	// Number of clusters
int num_dimensions = 100; // Number of dimension >= 1 -- Fixed parameter

int num_sparity = 100000; // Level of sparity of data points
double MostLowHighPos[MAX_DATA_DIMENSION][2];
double degreeOfMemb[MAX_DATA_POINTS][MAX_CLUSTERS];
double dataPoint[MAX_DATA_POINTS][MAX_DATA_DIMENSION];
double clusterCentre[MAX_CLUSTERS][MAX_DATA_DIMENSION];

#ifdef _WIN32
//#include <Windows.h>
#else
#include <sys/time.h>
#include <ctime>
#endif

/* Remove if already defined */
typedef long long int64; typedef unsigned long long uint64;

/* Returns the amount of milliseconds elapsed since the UNIX epoch. Works on both
* windows and linux. */

uint64 GetTimeMs64()
{
#ifdef _WIN32
	/* Windows */
	FILETIME ft;
	LARGE_INTEGER li;

	/* Get the amount of 100 nano seconds intervals elapsed since January 1, 1601 (UTC) and copy it
	* to a LARGE_INTEGER structure. */
	GetSystemTimeAsFileTime(&ft);
	li.LowPart = ft.dwLowDateTime;
	li.HighPart = ft.dwHighDateTime;

	uint64 ret = li.QuadPart;
	ret -= 116444736000000000LL; /* Convert from file time to UNIX epoch time. */
	ret /= 10000; /* From 100 nano seconds (10^-7) to 1 millisecond (10^-3) intervals */

	return ret;
#else
	/* Linux */
	struct timeval tv;

	gettimeofday(&tv, NULL);

	uint64 ret = tv.tv_usec;
	/* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
	ret /= 1000;

	/* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
	ret += (tv.tv_sec * 1000);

	return ret;
#endif
}

static int compare(const void* a, const void* b)
{
	int* arr1 = (int*)b;
	int* arr2 = (int*)a;
	int diff1 = arr1[0] - arr2[0];
	if (diff1) return diff1;
	return arr1[1] - arr2[1];
}

string ExePath() {
	char buffer[MAX_PATH];
	GetModuleFileName(NULL, buffer, MAX_PATH);
	string::size_type pos = string(buffer).find_last_of("\\/");
	return string(buffer).substr(0, pos);
}

void readDataFromFile(string file) {
	ifstream fileInput(file);

	if (fileInput.fail())
	{
		cout << "Cannot open file at " << file << endl;
		return;
	}

	while (!fileInput.eof())
	{
		char line[255];
		fileInput.getline(line, 255);
		cout << line << endl;
	}

	fileInput.close();
}

// Read data from file to vector
vector<unsigned> readDataFromFileToVector(string file) {
	std::vector<unsigned> numbers;
	ifstream inputFile(file);
	// Check if exists and then open the file.
	if (inputFile.good()) {
		// Push items into a vector
		int current_number = 0;
		while (inputFile >> current_number) {
			numbers.push_back(current_number);
		}

		// Close the file.
		inputFile.close();

		// Display the numbers read:
		cout << "The numbers are: ";
		for (unsigned count = 0; count < numbers.size(); count++) {
			cout << numbers[count] << " ";
		}

		cout << endl;
	}
	else {
		cout << "Error!";
		_exit(0);

	}
	return numbers;
}

// Create membership file
void createMembershipFile(ofstream &fileMembership) {
	int i, j;
	for (i = 0; i < num_data_points; i++) {
		fileMembership << "Point[" << i << "] ";
		for (j = 0; j < num_clusters; j++) {
			fileMembership << degreeOfMemb[i][j] << " ";
		}
		fileMembership << std::endl;
	}

	fileMembership.close();
}

ofstream fileMembership(string fileName) {
	string currentDirectory = ExePath();
	replace(currentDirectory.begin(), currentDirectory.end(), '\\', '/');

	string fullFilePath = currentDirectory + "/" + fileName;
	ofstream fileOutput;
	fileOutput.open(fullFilePath);
	if (fileOutput.fail())
	{
		std::cout << "Cannot open file at " << fullFilePath << std::endl;
		//return 0;
	}
	return fileOutput;
}
//// Write new line to plot.dat file
//void writeNewLineDataToPlotFile() {
//	string fileName;
//	string currentDirectory = ExePath();
//	replace(currentDirectory.begin(), currentDirectory.end(), '\\', '/');
//	fileName = "plot.dat";
//	string fullFilePath = currentDirectory + "/" + fileName;
//
//	string file = fullFilePath;
//	ofstream fileOutput(file, std::ios_base::app);
//
//	if (fileOutput.fail())
//	{
//		std::cout << "Cannot open file at " << file << std::endl;
//		return;
//	}
//
//	fileOutput << '\n';
//	fileOutput.close();
//}

// Get current folder of input.txt file
string fullInputFilePath()
{
	string currentDirectory = ExePath();
	replace(currentDirectory.begin(), currentDirectory.end(), '\\', '/');
	string fileName = "input.txt";
	string fullFilePath = currentDirectory + "/" + fileName;
	return fullFilePath;

}

// Get current folder of input.txt file
string fullInputFilePath(int indexOfInputFiles)
{
	string fileName;
	string Result;//string which will contain the result
	stringstream convert; // stringstream used for the conversion
	convert << indexOfInputFiles;//add the value of Number to the characters in the stream
	Result = convert.str();//set Result to the content of the stream

	string currentDirectory = ExePath();
	replace(currentDirectory.begin(), currentDirectory.end(), '\\', '/');
	//string fileName = "input.txt";
	fileName = "input" + Result + ".txt";
	string fullFilePath = currentDirectory + "/" + fileName;
	return fullFilePath;

}

// Get current folder of output.txt file
string fullOutputFilePath(int indexOfInputFiles)
{
	string currentDirectory = ExePath();
	replace(currentDirectory.begin(), currentDirectory.end(), '\\', '/');
	string fileName = "output.txt";
	string fullFilePath = currentDirectory + "/" + fileName;
	return fullFilePath;
}

// Create data set for all data points
void createDataset() {
	int i, j;
	for (i = 0; i < num_data_points; i++) {
		for (j = 0; j < num_dimensions; j++) {
			dataPoint[i][j] = 1 + (rand() % num_sparity);
		}
	}
}

// Create overview file including all information
void createDatasetInputFile(ofstream &fileInput, int numOfDataPoints, int numOfCluster, int dimension, int m, double epsilon) {
	/*fileInput << numOfDataPoints << " " << numOfCluster << " " << dimension << "\n";
	fileInput << m << " " << epsilon << "\n";*/
	fileInput << "Input.txt file contains coordinates of all of points" << "\n";
	
	int i, j;
	for (i = 0; i < numOfDataPoints; i++) {
		for (j = 0; j < dimension; j++) {
			//fileInput << 1 + (rand() % num_sparity) << " ";
			fileInput << dataPoint[i][j] << " ";
		}
		fileInput << std::endl;
	}

	fileInput.close();
}


void initializeDoM23(int fuzziness, double epsilon) {
	int i, j;
	int p, randomVal;
	double sum;

	//Open input file and get initial parameters for running
	string inputFile = fullInputFilePath();
	ifstream fin(inputFile);
	fin >> num_data_points >> num_clusters >> num_dimensions;
	fin >> fuzziness >> epsilon;
	// Each row is a dimension of MostLowHighPos array. This code return the 
	// lowest and highest point position of each dimension.
	for (i = 0; i < num_data_points; i++) {
		for (j = 0; j < num_dimensions; j++) {
			fin >> dataPoint[i][j];
			if (dataPoint[i][j] < MostLowHighPos[j][0])
				MostLowHighPos[j][0] = dataPoint[i][j];
			if (dataPoint[i][j] > MostLowHighPos[j][1])
				MostLowHighPos[j][1] = dataPoint[i][j];
		}
	}

	// Initialize degree of membership of a given data point i to cluster j.
	// Each has value from 0 to 1.
	for (i = 0; i < num_data_points; i++) {
		sum = 0.0;
		p = 100;
		for (j = 1; j < num_clusters; j++) {
			randomVal = rand() % (p + 1);
			p -= randomVal;
			degreeOfMemb[i][j] = randomVal / 100.0;
			sum += degreeOfMemb[i][j];
		}
		degreeOfMemb[i][0] = 1.0 - sum;
	}
}


void initializeDoM(int fuzziness, double epsilon) {
	int i, j;
	int p, randomVal;
	double sum;

	// Open input file and get initial parameters for running
	//string inputFile = fullInputFilePath();
	//ifstream fin(inputFile);
	//fin >> num_data_points >> num_clusters >> num_dimensions;
	//fin >> fuzziness >> epsilon;
	//// Each row is a dimension of MostLowHighPos array. This code return the 
	//// lowest and highest point position of each dimension.
	//for (i = 0; i < num_data_points; i++) {
	//	for (j = 0; j < num_dimensions; j++) {
	//		fin >> dataPoint[i][j];
	//			if (dataPoint[i][j] < MostLowHighPos[j][0])
	//				MostLowHighPos[j][0] = dataPoint[i][j];
	//			if (dataPoint[i][j] > MostLowHighPos[j][1])
	//				MostLowHighPos[j][1] = dataPoint[i][j];
	//	}
	//}

	// Initialize degree of membership of a given data point i to cluster j.
	// Each has value from 0 to 1.
	for (i = 0; i < num_data_points; i++) {
		sum = 0.0;
		p = 100;
		for (j = 1; j < num_clusters; j++) {
			randomVal = rand() % (p + 1);
			p -= randomVal;
			degreeOfMemb[i][j] = randomVal / 100.0;
			sum += degreeOfMemb[i][j];
		}
		degreeOfMemb[i][0] = 1.0 - sum;
	}
}

int calculateCentreVectors(int fuzziness) {
	int i, j, k;
	double numerator, denominator;
	double t[MAX_DATA_POINTS][MAX_CLUSTERS];

	// Calculate array of power degree of member of fuzziness
	for (i = 0; i < num_data_points; i++) {
		for (j = 0; j < num_clusters; j++) {
			t[i][j] = pow(degreeOfMemb[i][j], fuzziness);
		}
	}

	// Calculate cluster centre
	for (j = 0; j < num_clusters; j++) {
		for (k = 0; k < num_dimensions; k++) {
			numerator = 0.0;
			denominator = 0.0;
			for (i = 0; i < num_data_points; i++) {
				numerator += t[i][j] * dataPoint[i][k];
				denominator += t[i][j];
			}
			clusterCentre[j][k] = numerator / denominator;
		}
	}
	return 0;
}

// Get Norm
double getNorm(int i, int j) {
	int k;
	double total = 0.0;
	for (k = 0; k < num_dimensions; k++) {
		total += pow(dataPoint[i][k] - clusterCentre[j][k], 2);
	}
	return sqrt(total);
}

// Get new degree of membership
double getNewDoM(int fuzziness, int i, int j) {
	int k;
	double t, power;
	double total = 0.0;
	power = 2 / (fuzziness - 1);

	for (k = 0; k < num_clusters; k++) {
		t = getNorm(i, j) / getNorm(i, k);
		t = pow(t, power);
		total += t;
	}
	return (1.0 / total);
}

// Update degree of membership of all data points on each cluster
double updateDoM(int fuzziness) {
	int i, j;
	double diff;
	double max_diff = 0.0;
	double newDoM;

	// Loop for each cluster and data points then get maximum difference in degree of membership
	for (j = 0; j < num_clusters; j++) {
		for (i = 0; i < num_data_points; i++) {
			newDoM = getNewDoM(fuzziness, i, j);
			diff = newDoM - degreeOfMemb[i][j];
			diff = abs(diff);
			if (diff > max_diff)
				max_diff = diff;
			degreeOfMemb[i][j] = newDoM;
		}
	}
	return max_diff;
}

// Main function of fuzzy C-means
int fuzzyCMeans23(int fuzziness, double epsilon) {
	double max_difference;
	initializeDoM23(fuzziness, epsilon);
	do {
		calculateCentreVectors(fuzziness);
		max_difference = updateDoM(fuzziness);
	} while (max_difference > epsilon);

	return 0;
}

// Main function of fuzzy C-means
int fuzzyCMeans(int fuzziness, double epsilon, int &numOfSteps) {
	double max_difference;
	// Initialize degree of membership of a given data point i to cluster j.
	initializeDoM(fuzziness, epsilon);
	createMembershipFile(fileMembership("Initial_Mebership.txt"));
	do {
		calculateCentreVectors(fuzziness);
		max_difference = updateDoM(fuzziness);
		numOfSteps = numOfSteps + 1;
	} while (max_difference > epsilon);
	/*ofstream test = fileMembership("mebershipBeforeGetoutFuzzyCMeans.txt");
	createMembershipFile(test);*/
	return 0;
}

// Create input file for overivew all information
ofstream fileInput() {
	string currentDirectory = ExePath();
	replace(currentDirectory.begin(), currentDirectory.end(), '\\', '/');
	string fileName = "input.txt";
	string fullFilePath = currentDirectory + "/" + fileName;
	ofstream fileInput;
	fileInput.open(fullFilePath);
	if (fileInput.fail())
	{
		std::cout << "Cannot open file at " << fullFilePath << std::endl;
		//return 0;
	}
	return fileInput;
}



ofstream fileReport() {
	string currentDirectory = ExePath();
	replace(currentDirectory.begin(), currentDirectory.end(), '\\', '/');
	string fileName = "report.txt";
	string fullFilePath = currentDirectory + "/" + fileName;
	ofstream fileOutput;
	fileOutput.open(fullFilePath, std::ios_base::app);
	if (fileOutput.fail())
	{
		std::cout << "Cannot open file at " << fullFilePath << std::endl;
	}
	return fileOutput;
}




// Get cluster file path correspondence to cluster index
string clusterFilePath(int clusterIndex) {
	string fileName;
	string Result;//string which will contain the result
	stringstream convert; // stringstream used for the conversion
	convert << clusterIndex;
	Result = convert.str();//set Result to the content of the stream

	string currentDirectory = ExePath();
	replace(currentDirectory.begin(), currentDirectory.end(), '\\', '/');
	fileName = "Cluster_Num" + Result + ".txt";
	string fullFilePath = currentDirectory + "/" + fileName;
	return fullFilePath;
}

// Create cluster file
void createClusterFile(string file) {
	ofstream fileOutput(file);

	if (fileOutput.fail())
	{
		std::cout << "Cannot open file at " << file << std::endl;
		return;
	}
	fileOutput.close();
}

// Write one line to a file
void writeOneLineToFile(ofstream &fileReport, string str, bool append = false) {
	fileReport.precision(15);
	fileReport.setf(ios::fixed, ios::floatfield);
	fileReport.setf(ios::showpoint);
	fileReport << str;
	fileReport << std::endl;
	fileReport.close();
}

// Write one data point to one cluster file
void writeOneDataPointToClusterFile(const std::string& name, double d1, bool append = false) {
	std::ofstream outfile;
	if (append)
		outfile.open(name, std::ios_base::app);
	else
		outfile.open(name);
	outfile << d1 << " " ;
}

// Write one data point to one cluster file
void writeOneDataPointToClusterFileLineBreak(const std::string& name, bool append = false) {
	std::ofstream outfile;
	if (append)
		outfile.open(name, std::ios_base::app);
	else
		outfile.open(name);
	outfile << std::endl;
}


// Write d1 and d2 for 2 dimension
void filePutContents(const std::string& name, double d1, double d2, bool append = false) {
	std::ofstream outfile;
	if (append)
		outfile.open(name, std::ios_base::app);
	else
		outfile.open(name);
	outfile << d1 << " " << d2 << std::endl;
}

// Write d1, d2 and d3 for 3 dimension
void filePutContents3d(const std::string& name, double d1, double d2, double d3, bool append = false) {
	std::ofstream outfile;
	if (append)
		outfile.open(name, std::ios_base::app);
	else
		outfile.open(name);
	outfile << d1 << " " << d2 << " " << d3 << std::endl;
}

// Write new line to gnuplot.gp file
void writeNewLineDataToClusterFiles() {
	int clusterNum;
	int i, j, k;
	double highestDoM;
	for (i = 0; i < num_data_points; i++) {
		clusterNum = 0;
		highestDoM = 0.0;
		for (j = 0; j < num_clusters; j++) {
			if (degreeOfMemb[i][j] > highestDoM) {
				highestDoM = degreeOfMemb[i][j];
				clusterNum = j;
			}
		}
		/*if (num_dimensions == 2) {
			filePutContents(clusterFilePath(clusterNum), dataPoint[i][0], dataPoint[i][1], true);
		}
		if (num_dimensions == 3) {
			filePutContents3d(clusterFilePath(clusterNum), dataPoint[i][0], dataPoint[i][1], dataPoint[i][2], true);
		}*/
		for (k = 0; k < num_dimensions; k++) {
			writeOneDataPointToClusterFile(clusterFilePath(clusterNum), dataPoint[i][k], true);
		}
		writeOneDataPointToClusterFileLineBreak(clusterFilePath(clusterNum), true);
	}
}

// Remove file before running application
void removeFile(string fileName) {
	string currentDirectory = ExePath();
	replace(currentDirectory.begin(), currentDirectory.end(), '\\', '/');
	string fullFilePath = currentDirectory + "/" + fileName;

	string file = fullFilePath;
	remove(file.c_str());

	//if (remove(file.c_str()) == 0)
		//perror("Error deleting filename file");
	//else
		//puts("File successfully deleted");
}

// Remove gnuplot.gp file before running application
void removePlotFile() {
	string fileName;
	string currentDirectory = ExePath();
	replace(currentDirectory.begin(), currentDirectory.end(), '\\', '/');
	fileName = "gnuplot.gp";
	string fullFilePath = currentDirectory + "/" + fileName;

	string file = fullFilePath;

	if (remove(file.c_str()) != 0)
		perror("Error deleting gnuplot.gp file");
	else
		puts("File successfully deleted");
}

// Create gnuplot.gp file
void createGnuPlotFile() {
	int i;
	string fileName;
	string currentDirectory = ExePath();
	replace(currentDirectory.begin(), currentDirectory.end(), '\\', '/');
	fileName = "gnuplot.gp";
	string fullFilePath = currentDirectory + "/" + fileName;

	string file = fullFilePath;
	ofstream fileOutput(file, std::ios_base::app);

	if (fileOutput.fail())
	{
		std::cout << "Cannot open file at " << file << std::endl;
		return;
	}

	fileOutput << "set terminal png medium" << std::endl;
	fileOutput << "set output \"clusterPlot.png\"" << std::endl;
	fileOutput << "set title \"Implementation for Fuzzy c-Means Algorithm\"" << std::endl;
	fileOutput << "set xlabel \"x-coordinate\"" << std::endl;
	fileOutput << "set ylabel \"y-coordinate\"" << std::endl;
	fileOutput << "set xrange [" << MostLowHighPos[0][0] << " : " << MostLowHighPos[0][1] << "]" << std::endl;
	fileOutput << "set yrange [" << MostLowHighPos[1][0] << " : " << MostLowHighPos[1][1] << "]" << std::endl;
	fileOutput << "plot 'Cluster_Num0.txt' using 1:2 with points pt 7 ps 1 lc 1 notitle";
	for (i = 1; i < num_clusters; i++) {
		fileOutput << ", \\" << std::endl;
		fileOutput << "'Cluster_Num" << i << ".txt' using 1:2 with points pt 7 ps 1 lc " << i + 1 << " notitle";
	}

	fileOutput.close();
}

// Create gnuplot.gp 3d file
void createGnuPlotFile3d() {

	int i;
	string fileName;
	string currentDirectory = ExePath();
	replace(currentDirectory.begin(), currentDirectory.end(), '\\', '/');
	fileName = "gnuplot.gp";
	string fullFilePath = currentDirectory + "/" + fileName;

	string file = fullFilePath;
	ofstream fileOutput(file, std::ios_base::app);

	if (fileOutput.fail())
	{
		std::cout << "Cannot open file at " << file << std::endl;
		return;
	}

	fileOutput << "set view equal xyz" << std::endl;
	fileOutput << "set xyplane at 0" << std::endl;
	fileOutput << "set title \"Implementation for Fuzzy c-Means Algorithm\"" << std::endl;
	fileOutput << "set border 4095" << std::endl;
	fileOutput << "set xlabel \"x-coordinate\"" << std::endl;
	fileOutput << "set ylabel \"y-coordinate\"" << std::endl;
	fileOutput << "set zlabel \"z-coordinate\"" << std::endl;
	fileOutput << "set xrange [" << MostLowHighPos[0][0] << " : " << MostLowHighPos[0][1] << "]" << std::endl;
	fileOutput << "set yrange [" << MostLowHighPos[1][0] << " : " << MostLowHighPos[1][1] << "]" << std::endl;
	fileOutput << "set zrange [" << MostLowHighPos[2][0] << " : " << MostLowHighPos[2][1] << "]" << std::endl;
	fileOutput << "splot 'Cluster_Num0.txt' using 1:2:3 with points pt 7 ps 1 lc 1 notitle";
	for (i = 1; i < num_clusters; i++) {
		fileOutput << ", \\" << std::endl;
		fileOutput << "'Cluster_Num" << i << ".txt' using 1:2:3 with points pt 7 ps 1 lc " << i + 1 << " notitle";
	}

	fileOutput.close();
}

// Create report file with multiple parameters
void createReportFileStep(ofstream &fileReport, int m, double epsilon, uint64 runningTime, int numOfSteps) {
	fileReport.precision(15);
	fileReport.setf(ios::fixed, ios::floatfield);
	fileReport.setf(ios::showpoint);
	fileReport << m << " " << epsilon << " " << runningTime << " " << numOfSteps;
	fileReport << std::endl;
	fileReport.close();
}

void createReportFile(ofstream &fileReport, int m, double epsilon, uint64 runningTime) {
	fileReport.precision(15);
	fileReport.setf(ios::fixed, ios::floatfield);
	fileReport.setf(ios::showpoint);
	fileReport << m << " " << epsilon << " " << runningTime;
	fileReport << std::endl;
	fileReport.close();
}

void cleanUpFiles() {
	removeFile("input.txt");
	removeFile("membership.txt");
	removeFile("Initial_Mebership.txt");
	removeFile("report.txt");

	int maxCluster = MAX_CLUSTERS;
	string fileName;
	string Result;//string which will contain the result
	stringstream convert; // stringstream used for the conversion
	for (int i = 0; i < maxCluster; i++) {
		convert << i;
		Result = convert.str();//set Result to the content of the stream
		fileName = "Cluster_Num" + Result + ".txt";
		removeFile(fileName);
	}
}
int main()
{
	cleanUpFiles();
	cout << "\n" << "Generate random data points" << "\n";

	// Create data set for all points
	createDataset();
	/*double epsilon = 0.00005;
	int fuzziness = 300;*/
	double epsilon = 0.00005;
	int fuzziness = 2;
	int numOfSteps;
	writeOneLineToFile(fileReport(), "Fuzziness = Epsilon = RunningTime = Steps", false);
	/*for (fuzziness = 2; fuzziness < 3; fuzziness++) {
		for (epsilon = 0.0000005; epsilon <= 0.0005; epsilon = epsilon * 2) {*/
			numOfSteps = 0;
			createDatasetInputFile(fileInput(), num_data_points, num_clusters, num_dimensions, fuzziness, epsilon);
			uint64 startTime, finishTime, runningTime;
			startTime = GetTimeMs64();
			fuzzyCMeans(fuzziness, epsilon, numOfSteps);
			finishTime = GetTimeMs64();
			runningTime = finishTime - startTime;
			createReportFileStep(fileReport(), fuzziness, epsilon, runningTime, numOfSteps);

			// //Without Steps
			////createReportFile(fileReport(), fuzziness, epsilon, runningTime);

			createMembershipFile(fileMembership("mebership.txt"));
			//// Create cluster files
			//for (int i = 0; i < num_clusters; i++) {
			//	createClusterFile(clusterFilePath(i));
			//}

			//// Write data to cluster files
			//writeNewLineDataToClusterFiles();
		//}
	//}

	


	return 0;
}


