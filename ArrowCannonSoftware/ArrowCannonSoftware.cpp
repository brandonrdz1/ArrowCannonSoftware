#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <thread>
#include <mutex>
#include <limits>
#include <chrono>
#include <atomic>

#include "Entity.h"

using namespace std;

#define velToIndex 15.7056868866
#define indexToVel 0.0636712044
#define centerX -86214.5
#define cneterZ 87930.5

void printArray(double gains[], int size) {
    cout << "Gains: {";
    for (int i = 0; i < size; i++) {
        cout << gains[i] << ((i != size - 1) ? ", " : "");
    }
    cout << "}" << endl;
}

void printVector(vector<int> vector) {
    cout << "{";
    for (int i = 0; i < vector.size(); i++) {
        cout << vector[i] << ((i != vector.size() - 1) ? ", " : "");
    }
    cout << "}" << endl;
    return;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    os << "{";
    for (size_t i = 0; i < vec.size(); ++i) {
        os << vec[i];
        if (i != vec.size() - 1) {
            os << ", ";
        }
    }
    os << "}";
    return os;
}

bool generateNextCombination(vector<int>& currentVector, const vector<int>& finalVector) {
    size_t i = currentVector.size();

    while (i > 0) {
        --i;
        if (currentVector[i] < finalVector[i]) {
            ++currentVector[i];
            fill(currentVector.begin() + i + 1, currentVector.end(), 0);
            return true;
        }
    }
    return false;
}

std::atomic<double> globalRecord(1000.0);

double xVelFromExplosion(double projX, double powX) {
	return (1.0 - (projX-powX) / 8.0);
}

double xArrowFF(double startX, double startU, int ticks) {
	return startX + 100 * (1 - pow(0.99, ticks)) * startU;
}

int sumVectorElements(vector<int>& vector) {
	int sum = 0;
	for (int element : vector) {
		sum += element;
	}
	return sum;
}
template<typename T>
int numOfPopulatedElements(vector<vector<T>>& vector2D) {
	int populatedCount = 0;
	for (const auto& vec : vector2D) {
		if (!vec.empty()) {
			populatedCount++;
		}
	}
	// std::cout << populatedCount << std::endl;
	return populatedCount;
}

double decodeRatioU(const vector<int>& ratio, double* apBoostersArr, const int& size) {
	Tnt aPowerStart("-86228.49000000954 278.06125000119226 87930.5 -1.66683193393 19.99749994821824 0.0");
	Tnt aPowerFinish = aPowerStart;
	aPowerFinish.amount = 1296;

	Arrow centerArrowFinish("-86214.5 297.9500000087917 87930.5 0.0 0.8806000169232489 0.0");
	// std::cout << "testing ratio: " << ratio << std::endl;
	/*
	for (int i = 0; i < size; i++) {
		cout << " apBoostersArr[" << i << "]: " << apBoostersArr[i] << " * " << ratio[i] << std::endl;
	}
	*/
	// if the power affects the arrow
	aPowerFinish.setX(
		aPowerStart.getU() + aPowerStart.getX()
		+ apBoostersArr[0] * ratio[0]
		+ apBoostersArr[1] * ratio[1]
		+ apBoostersArr[2] * ratio[2]
		+ apBoostersArr[3] * ratio[3]
		+ apBoostersArr[4] * ratio[4]
		+ apBoostersArr[5] * ratio[5]
		+ apBoostersArr[6] * ratio[6]
		+ apBoostersArr[7] * ratio[7]
		+ apBoostersArr[8] * ratio[8]
		+ apBoostersArr[9] * ratio[9]
		+ apBoostersArr[10] * ratio[10]
		+ apBoostersArr[11] * ratio[11]);
	aPowerFinish.setY(
		aPowerStart.getV()
		+ aPowerStart.getY()
		+ Tnt::gravity);
	// restart the arrow finish
	centerArrowFinish.explosion(&aPowerFinish); // this can be obviously optimized
	centerArrowFinish.freefall(17);
	return centerArrowFinish.getU();
}

double decodeRatioX(const vector<int>& ratio, double* apBoostersArr, const int& size) {
	Tnt aPowerStart("-86228.49000000954 278.06125000119226 87930.5 -1.66683193393 19.99749994821824 0.0");
	Tnt aPowerFinish = aPowerStart;
	aPowerFinish.amount = 1296;

	Arrow centerArrowFinish("-86214.5 297.9500000087917 87930.5 0.0 0.8806000169232489 0.0");
	// std::cout << "testing ratio: " << ratio << std::endl;
	/*
	for (int i = 0; i < size; i++) {
		cout << " apBoostersArr[" << i << "]: " << apBoostersArr[i] << " * " << ratio[i] << std::endl;
	}
	*/
	// if the power affects the arrow
	aPowerFinish.setX(
		aPowerStart.getU() + aPowerStart.getX()
		+ apBoostersArr[0] * ratio[0]
		+ apBoostersArr[1] * ratio[1]
		+ apBoostersArr[2] * ratio[2]
		+ apBoostersArr[3] * ratio[3]
		+ apBoostersArr[4] * ratio[4]
		+ apBoostersArr[5] * ratio[5]
		+ apBoostersArr[6] * ratio[6]
		+ apBoostersArr[7] * ratio[7]
		+ apBoostersArr[8] * ratio[8]
		+ apBoostersArr[9] * ratio[9]
		+ apBoostersArr[10] * ratio[10]
		+ apBoostersArr[11] * ratio[11]);
	aPowerFinish.setY(
		aPowerStart.getV() 
		+ aPowerStart.getY() 
		+ Tnt::gravity);
	// restart the arrow finish
	centerArrowFinish.explosion(&aPowerFinish); // this can be obviously optimized
	centerArrowFinish.freefall(17);
	return centerArrowFinish.getX();
}

double fromRatioGetU(const vector<double>& ratio) {
	double velocity = 0.0;
	for (size_t i = 0; i < 12; i++) {
		velocity += ratio[i] * ratio[i + 12];
	}
	return velocity;
}

void runSimulation1(vector<int> parameters, mutex& mtx, unsigned int id) {
	std::cout << "starting id " << id << std::endl;

	Arrow centerArrowStart("-86214.5 297.9500000087917 87930.5 0.0 0.8806000169232489 0.0");
	Arrow centerArrowFinish = centerArrowStart;
	
	Tnt arrowPower("-86214.5 297.9500000087917 87930.5 0.0 0.8806000169232489 0.0");
	arrowPower.setX(centerArrowStart.getX() - 8.0);
	arrowPower.setY(centerArrowStart.getY() + Arrow::explosion_height - Tnt::explosion_height);
	double exitVel;
	arrowPower.amount = 1296;
	double step = 0.001;
	int counter = 0;
	while (arrowPower.getX() < centerArrowStart.getX()) {
		arrowPower.setX(arrowPower.getX() + step);
		centerArrowFinish = centerArrowStart;

		centerArrowFinish.explosion(&arrowPower);
		exitVel = centerArrowFinish.getU();
		centerArrowFinish.freefall(17);
		//std::cout << centerArrowStart.getX()- arrowPower.getX() << ", " << counter << std::endl;
		//std::cout << centerArrowStart.getX() - arrowPower.getX() << ", " << centerArrowFinish.getX() - centerX << std::endl;
		std::cout << centerArrowStart.getX() - arrowPower.getX() << ", " << exitVel << std::endl;
		counter++;
	}

}

void runSimulationArrows(vector<int> parameters, vector<vector<double>>& ratios,  mutex& mtxText, mutex& mtxOutput, unsigned int id) {
	int maxRangeIndex = ratios.size();
    {
        std::lock_guard<std::mutex> lock(mtxText);
        std::cout << "Starting [" << id << "] parameters: " << parameters << std::endl;
    }
	for (int param : parameters) {
		{
			std::lock_guard<std::mutex> lock(mtxText);
			std::cout << "Running Simulation [" << id << "] with current parameter: " << param << std::endl;
		}
	
		// Find the most optimal location for the arrows
		const int numOfArrowPowerAdj = 12;
		Tnt a1("-86234.66023655233 278.0 87930.5");
		Tnt a2("-86234.92907985931 278.0 87930.5");
		Tnt a3("-86235.21242460351 278.0 87930.5");
		Tnt a4("-86235.50970078005 278.0 87930.5");
		Tnt a5("-86235.82027033666 278.0 87930.5");
		Tnt a6("-86236.14342772469 278.0 87930.5");
		Tnt a7("-86234.53621773126 278.0 87930.5");
		Tnt a8("-86234.8832512168 278.0 87930.5");
		Tnt a9("-86235.2447882072 278.0 87930.5");
		Tnt a10("-86235.62025871217 278.0 87930.5");
		Tnt a11("-86236.00902470213 278.0 87930.5");
		Tnt a12("-86236.41038066245 278.0 87930.5");
		Tnt apBoosters[] = { a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12 };

		Tnt aPowerStart("-86228.49000000954 278.06125000119226 87930.5 -1.66683193393 19.99749994821824 0.0");
		Tnt aPowerFinish = aPowerStart;  // initialzies z
		{
			cout << "Arrow Power location without any booster: " << std::endl;
			Tnt aPowerCopy = aPowerStart;
			aPowerCopy.print();
			aPowerCopy.freefall(1);
			aPowerCopy.print();
			cout << endl;
		}

		// X-gains for the arrow power's 12 boosters
		double* apBoostersArr = aPowerStart.gainArray(apBoosters, numOfArrowPowerAdj, 1.0, direction::xDir);
		for (int i = 0; i < numOfArrowPowerAdj; i++) {
			std::cout << "aexposureX[" << i << "]: " << apBoostersArr[i] << std::endl;
		} std::cout << endl;
		// Y-gains for the arrow power's 12 boosters
		double* apBoostersArrY = aPowerStart.gainArray(apBoosters, numOfArrowPowerAdj, 1.0, direction::yDir);
		for (int i = 0; i < numOfArrowPowerAdj; i++) {
			std::cout << "aexposureY[" << i << "]: " << apBoostersArrY[i] << std::endl;
		} std::cout << endl;

		// Spawn Arrows
		Arrow centerArrowStart("-86214.5 297.9500000087917 87930.5 0.0 0.8806000169232489 0.0");
		Arrow centerArrowFinish = centerArrowStart;


		// divide 1273.424088 into 20k blocks
		// Initializing stuff
		long int counter = 0;
		int usableCounter = 0;
		int dispenserSectionSize = 36;
		int dispenserSections = 36;
		double aPowerLocX; // powers horizontal coordinate (this is probably just going to be x)
		double aPowerLocMotStart = aPowerStart.getX() + aPowerStart.getU(); // saves on an addition
		double aPowerLocY = aPowerStart.getY() + aPowerStart.getV() + Tnt::gravity; // y coordinate of the powers
		std::cout << "aPower yLoc: " << aPowerLocY << std::endl;
		aPowerFinish.setY(aPowerLocY);
		aPowerFinish.amount = 1296;

		// stuff to save/output
		unsigned int index;
		unsigned long int swaps = 0;
		double currentScore, existingScore;

		// Initializing ratios
		vector<int> startVector = { param, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		int max = 6;
		vector<int> finalVector = { param, max, max, max, max, max, max, max, max, max, max, max + 1 };
		vector<int> currentVector = startVector;
		do {
			if (currentVector[11] >= max) { continue; }
			counter++;
			if (counter % 1000000000 == 0) {
				{
					std::lock_guard<std::mutex> lock(mtxOutput);
					std::cout << "id " << id << " at " << currentVector << " swaps: " << swaps << endl;
				}
			}
			aPowerLocX = aPowerLocMotStart
				+ apBoostersArr[0] * currentVector[0]
				+ apBoostersArr[1] * currentVector[1]
				+ apBoostersArr[2] * currentVector[2]
				+ apBoostersArr[3] * currentVector[3]
				+ apBoostersArr[4] * currentVector[4]
				+ apBoostersArr[5] * currentVector[5]
				+ apBoostersArr[6] * currentVector[6]
				+ apBoostersArr[7] * currentVector[7]
				+ apBoostersArr[8] * currentVector[8]
				+ apBoostersArr[9] * currentVector[9]
				+ apBoostersArr[10] * currentVector[10]
				+ apBoostersArr[11] * currentVector[11];
			if (aPowerLocX > centerX - 0.05018546 || aPowerLocX - centerX < -7.9) {// 20000.01vel // 0.88 was good 0.05 for max range
				continue;
			}
			aPowerFinish.setX(aPowerLocX);
			// restart the arrow finish
			centerArrowFinish = centerArrowStart;
			centerArrowFinish.explosion(&aPowerFinish); 
			// compute the index
			index = (int)round(centerArrowFinish.getU() * velToIndex); // applied velocity to corresponding index
			if (index >= maxRangeIndex) { continue; } // if too much range
			
			
			// Create new entry
			if (ratios[index].empty()) { // populate elements that are empty
				{
					std::lock_guard<std::mutex> lock(mtxOutput);
					ratios[index].resize(14);
					for (size_t i = 0; i < 12; i++) {
						ratios[index][i] = static_cast<double>(currentVector[i]);
					}
					ratios[index][12] = centerArrowFinish.getU();
					ratios[index][13] = aPowerFinish.getX();
					std::cout << "New element added at " << index << ": " << currentVector << " full: " << numOfPopulatedElements(ratios) << "/" << maxRangeIndex << " swaps: " << swaps << std::endl;
				}
				
				continue;
			}
			// Check for improvements
			existingScore = abs((ratios[index][12]*velToIndex) - index);
			currentScore = abs((centerArrowFinish.getU()*velToIndex) - index);
			if (currentScore < existingScore) {
				{
					std::lock_guard<std::mutex> lock(mtxOutput);
					// display the swap locations
					/*std::cout << "swapping[" << index << "]: "
						<< ratios[index][12] * velToIndex - centerX << " -> "
						<< centerArrowFinish.getU() * velToIndex - centerX << std::endl;
					*/
					// No need to resize
					for (size_t i = 0; i < 12; i++) {
						ratios[index][i] = static_cast<double>(currentVector[i]);
					}
					ratios[index][12] = centerArrowFinish.getU();
					ratios[index][13] = aPowerFinish.getX();
				}
				swaps++;
			}
			else {
				continue;
			}

		} while (generateNextCombination(currentVector, finalVector) && currentVector != finalVector);

		//delete[] hpBoostersArr;
		//delete[] hpBoostersArrY;
		delete[] apBoostersArr;
		delete[] apBoostersArrY;
		return;
    }
	{
		std::lock_guard<std::mutex> lock(mtxText);
		std::cout << "Finished [" << id << "] parameters: " << parameters << std::endl;
	}
}

void runSimulationHammer(vector<int> parameters, vector<vector<double>>& ratios, mutex& mtxText, mutex& mtxOutput, unsigned int id) {
	int maxRangeIndex = ratios.size();
	{
		std::lock_guard<std::mutex> lock(mtxText);
		std::cout << "Starting [" << id << "] parameters: " << parameters << std::endl;
	}
	for (int param : parameters) {
		{
			std::lock_guard<std::mutex> lock(mtxText);
			std::cout << "Running Simulation [" << id << "] with current parameter: " << param << std::endl;
		}

		// Find the most optimal location for the arrows
		const int numOfArrowPowerAdj = 12;
		Tnt a1("-86234.66023655233 278.0 87930.5");
		Tnt a2("-86234.92907985931 278.0 87930.5");
		Tnt a3("-86235.21242460351 278.0 87930.5");
		Tnt a4("-86235.50970078005 278.0 87930.5");
		Tnt a5("-86235.82027033666 278.0 87930.5");
		Tnt a6("-86236.14342772469 278.0 87930.5");
		Tnt a7("-86234.53621773126 278.0 87930.5");
		Tnt a8("-86234.8832512168 278.0 87930.5");
		Tnt a9("-86235.2447882072 278.0 87930.5");
		Tnt a10("-86235.62025871217 278.0 87930.5");
		Tnt a11("-86236.00902470213 278.0 87930.5");
		Tnt a12("-86236.41038066245 278.0 87930.5");
		Tnt apBoosters[] = { a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12 };

		Tnt aPowerStart("-86228.49000000954 278.06125000119226 87930.5 -1.66683193393 19.99749994821824 0.0");
		Tnt aPowerFinish = aPowerStart;  // initialzies z
		{
			cout << "Arrow Power location without any booster: " << std::endl;
			Tnt aPowerCopy = aPowerStart;
			aPowerCopy.print();
			aPowerCopy.freefall(1);
			aPowerCopy.print();
			cout << endl;
		}

		// X-gains for the arrow power's 12 boosters
		double* apBoostersArr = aPowerStart.gainArray(apBoosters, numOfArrowPowerAdj, 1.0, direction::xDir);
		for (int i = 0; i < numOfArrowPowerAdj; i++) {
			std::cout << "aexposureX[" << i << "]: " << apBoostersArr[i] << std::endl;
		} std::cout << endl;
		// Y-gains for the arrow power's 12 boosters
		double* apBoostersArrY = aPowerStart.gainArray(apBoosters, numOfArrowPowerAdj, 1.0, direction::yDir);
		for (int i = 0; i < numOfArrowPowerAdj; i++) {
			std::cout << "aexposureY[" << i << "]: " << apBoostersArrY[i] << std::endl;
		} std::cout << endl;

		// Spawn Arrows
		Arrow centerArrowStart("-86214.5 297.9500000087917 87930.5 0.0 0.8806000169232489 0.0");
		Arrow centerArrowFinish = centerArrowStart;


		// divide 1273.424088 into 20k blocks
		// Initializing stuff
		long int counter = 0;
		int usableCounter = 0;
		int dispenserSectionSize = 36;
		int dispenserSections = 36;
		double aPowerLocX; // powers horizontal coordinate (this is probably just going to be x)
		double aPowerLocMotStart = aPowerStart.getX() + aPowerStart.getU(); // saves on an addition
		double aPowerLocY = aPowerStart.getY() + aPowerStart.getV() + Tnt::gravity; // y coordinate of the powers
		std::cout << "aPower yLoc: " << aPowerLocY << std::endl;
		aPowerFinish.setY(aPowerLocY);
		aPowerFinish.amount = 1296;

		// stuff to save/output
		unsigned int index;
		unsigned long int swaps = 0;
		double currentScore, existingScore;

		// Initializing ratios
		vector<int> startVector = { param, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		int max = 6;
		vector<int> finalVector = { param, max, max, max, max, max, max, max, max, max, max, max + 1 };
		vector<int> currentVector = startVector;
		do {
			if (currentVector[11] >= max) { continue; }
			counter++;
			if (counter % 1000000000 == 0) {
				{
					std::lock_guard<std::mutex> lock(mtxOutput);
					std::cout << "id " << id << " at " << currentVector << " swaps: " << swaps << endl;
				}
			}
			aPowerLocX = aPowerLocMotStart
				+ apBoostersArr[0] * currentVector[0]
				+ apBoostersArr[1] * currentVector[1]
				+ apBoostersArr[2] * currentVector[2]
				+ apBoostersArr[3] * currentVector[3]
				+ apBoostersArr[4] * currentVector[4]
				+ apBoostersArr[5] * currentVector[5]
				+ apBoostersArr[6] * currentVector[6]
				+ apBoostersArr[7] * currentVector[7]
				+ apBoostersArr[8] * currentVector[8]
				+ apBoostersArr[9] * currentVector[9]
				+ apBoostersArr[10] * currentVector[10]
				+ apBoostersArr[11] * currentVector[11];
			if (aPowerLocX > centerX - 0.05018546 || aPowerLocX - centerX < -7.9) {// 20000.01vel // 0.88 was good 0.05 for max range
				continue;
			}
			aPowerFinish.setX(aPowerLocX);
			// restart the arrow finish
			centerArrowFinish = centerArrowStart;
			centerArrowFinish.explosion(&aPowerFinish);
			// compute the index
			index = (int)round(centerArrowFinish.getU() * velToIndex); // applied velocity to corresponding index
			if (index >= maxRangeIndex) { continue; } // if too much range


			// Create new entry
			if (ratios[index].empty()) { // populate elements that are empty
				{
					std::lock_guard<std::mutex> lock(mtxOutput);
					ratios[index].resize(14);
					for (size_t i = 0; i < 12; i++) {
						ratios[index][i] = static_cast<double>(currentVector[i]);
					}
					ratios[index][12] = centerArrowFinish.getU();
					ratios[index][13] = aPowerFinish.getX();
					std::cout << "New element added at " << index << ": " << currentVector << " full: " << numOfPopulatedElements(ratios) << "/" << maxRangeIndex << " swaps: " << swaps << std::endl;
				}

				continue;
			}
			// Check for improvements
			existingScore = abs((ratios[index][12] * velToIndex) - index);
			currentScore = abs((centerArrowFinish.getU() * velToIndex) - index);
			if (currentScore < existingScore) {
				{
					std::lock_guard<std::mutex> lock(mtxOutput);
					// display the swap locations
					/*std::cout << "swapping[" << index << "]: "
						<< ratios[index][12] * velToIndex - centerX << " -> "
						<< centerArrowFinish.getU() * velToIndex - centerX << std::endl;
					*/
					// No need to resize
					for (size_t i = 0; i < 12; i++) {
						ratios[index][i] = static_cast<double>(currentVector[i]);
					}
					ratios[index][12] = centerArrowFinish.getU();
					ratios[index][13] = aPowerFinish.getX();
				}
				swaps++;
			}
			else {
				continue;
			}

		} while (generateNextCombination(currentVector, finalVector) && currentVector != finalVector);

		//delete[] hpBoostersArr;
		//delete[] hpBoostersArrY;
		delete[] apBoostersArr;
		delete[] apBoostersArrY;
		return;
	}
	{
		std::lock_guard<std::mutex> lock(mtxText);
		std::cout << "Finished [" << id << "] parameters: " << parameters << std::endl;
	}
}

int main() {
    cout << setprecision(17);
    cout << setfill(' ');
    //cout << fixed;

	// Initialize mutex and outputRatios for arrows and hammer
    mutex mtxText;
	mutex mtxOutput;
    vector<thread> threads;
	const int maxRangeIndex = 2001;// 9001
	vector<vector<double>> arrowRatios(maxRangeIndex); // index actually corresponds to the range, index 20,000 is 19,999 blocks

	// Obtain number of threads
	unsigned int numThreads = thread::hardware_concurrency()-3;
    if (numThreads == 0) {
        numThreads = 2; // Default to 2 if hardware_concurrency returns 0
    }
    std::cout << "Number of threads: " << numThreads << std::endl;

	// Begin parameters to find arrow ratios
    vector<vector<int>> parameters(numThreads);
    int begin = 0;
    int end = 12;
    int counter = 0;
    for (int i = begin; i <= end; i++) {
        parameters[counter].push_back(i);
        counter++;
        if (counter == numThreads) { counter = 0; }
    }
    for (unsigned int i = 0; i < numThreads; ++i) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        threads.emplace_back(runSimulationArrows, parameters[i], ref(arrowRatios), ref(mtxText), ref(mtxOutput), i);
    }
	// joins all threads
    for (auto& th : threads) { th.join(); }

	// Output the ratios found for arrows
	for (int i = 0; i < arrowRatios.size(); i++) {
		cout << "range " << i << ": " << arrowRatios[i];
		if (!arrowRatios[i].empty())
			cout << "\tpositions: " << arrowRatios[i][12] * velToIndex + centerX << std::endl;
		else
			std::cout << std::endl;
	}


	// Clear or reset variables for hammer ratios




    return 0;
}
