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
		if (i == 11 || i==14) os << "\t";
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

void runSimulationArrows(vector<int> parameters, vector<vector<double>>& ratios,  mutex& mtx, unsigned int id) {
	int maxRangeIndex = ratios.size();
    {
        std::lock_guard<std::mutex> lock(mtx);
        std::cout << "Starting [" << id << "] parameters: " << parameters << std::endl;
    }
	for (int param : parameters) {
		{
			std::lock_guard<std::mutex> lock(mtx);
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


		double* apBoostersArr = nullptr;
		double* apBoostersArrY = nullptr;
		{
			std::lock_guard<std::mutex> lock(mtx);
		
			// X-gains for the arrow power's 12 boosters
			apBoostersArr = aPowerStart.gainArray(apBoosters, numOfArrowPowerAdj, 1.0, direction::xDir);
			for (int i = 0; i < numOfArrowPowerAdj; i++) {
				std::cout << "aexposureX[" << i << "]: " << apBoostersArr[i] << std::endl;
			} std::cout << endl;
			// Y-gains for the arrow power's 12 boosters
			apBoostersArrY = aPowerStart.gainArray(apBoosters, numOfArrowPowerAdj, 1.0, direction::yDir);
			for (int i = 0; i < numOfArrowPowerAdj; i++) {
				std::cout << "aexposureY[" << i << "]: " << apBoostersArrY[i] << std::endl;
			} std::cout << endl;
		}
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
		{
			std::lock_guard<std::mutex> lock(mtx);
			std::cout << "id " << id << " starts at " << currentVector << " full: " << numOfPopulatedElements(ratios) << "/" << maxRangeIndex << " swaps: " << swaps << endl;
		}
		do {
			if (currentVector[11] > max) { continue; }
			counter++;
			if (counter % 1000000000 == 0) {
				{
					std::lock_guard<std::mutex> lock(mtx);
					std::cout << "id " << id << " at " << currentVector << " full: " << numOfPopulatedElements(ratios) << "/" << maxRangeIndex << " swaps: " << swaps << endl;
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
					std::lock_guard<std::mutex> lock(mtx);
					ratios[index].resize(14);
					for (size_t i = 0; i < 12; i++) {
						ratios[index][i] = static_cast<double>(currentVector[i]);
					}
					ratios[index][12] = centerArrowFinish.getU();
					ratios[index][13] = aPowerFinish.getX();
					//std::cout << "New element added at " << index << ": " << currentVector << " full: " << numOfPopulatedElements(ratios) << "/" << maxRangeIndex << " swaps: " << swaps << std::endl;
				}
				
				continue;
			}
			// Check for improvements
			existingScore = abs((ratios[index][12]*velToIndex) - index);
			currentScore = abs((centerArrowFinish.getU()*velToIndex) - index);
			if (currentScore < existingScore) {
				{
					std::lock_guard<std::mutex> lock(mtx);
					// display the swap locations
					/*std::cout << "swapping[" << index << "]: "
						<< ratios[index][12] * velToIndex - centerX << " -> "
						<< centerArrowFinish.getU() * velToIndex - centerX << std::endl;
					*/
					// No need to resize
					for (size_t i = 0; i < 12; i++) {
						ratios[index][i] = static_cast<double>(currentVector[i]);
					}
					ratios[index][12] = centerArrowFinish.getU(); // initial velocity
					ratios[index][13] = aPowerFinish.getX(); //final position of the arrow power
				}
				swaps++;
			}
			else {
				continue;
			}

		} while (generateNextCombination(currentVector, finalVector) && currentVector != finalVector);

		delete[] apBoostersArr;
		delete[] apBoostersArrY;
		
    }
	{
		std::lock_guard<std::mutex> lock(mtx);
		std::cout << "Finished [" << id << "] parameters: " << parameters << std::endl;
	}
	return;
}

void runSimulationHammer(vector<int> parameters, vector<vector<double>>& arrowRatios, vector<vector<double>>& hammerRatios, mutex& mtx, unsigned int id) {
	int maxRangeIndex = arrowRatios.size();
	{
		std::lock_guard<std::mutex> lock(mtx);
		std::cout << "Starting [" << id << "] parameters: " << parameters << std::endl;
	}
	for (int param : parameters) {
		{
			std::lock_guard<std::mutex> lock(mtx);
			std::cout << "Running Simulation [" << id << "] with current parameter: " << param << std::endl;
		}

		// Find the most optimal location for the hammer
		const int numOfHammerPowerAdj = 12;
		Tnt h1("-86234.78586161736 278.0 87930.5");
		Tnt h2("-86235.03129429705 278.0 87930.5");
		Tnt h3("-86235.29123018472 278.0 87930.5");
		Tnt h4("-86235.56509943453 278.0 87930.5");
		Tnt h5("-86235.85226417152 278.0 87930.5");
		Tnt h6("-86236.15201904518 278.0 87930.5");
		Tnt h7("-86233.91205295271 278.0 87930.5");
		Tnt h8("-86234.25911295973 278.0 87930.5");
		Tnt h9("-86234.62067852802 278.0 87930.5");
		Tnt h10("-86234.99617983929 278.0 87930.5");
		Tnt h11("-86235.38497905637 278.0 87930.5");
		Tnt h12("-86235.78637088025 278.0 87930.5");
		Tnt hpBoosters[] = { h1, h2, h3, h4, h5, h6, h7, h8, h9, h10, h11, h12 };

		Tnt hPowerStart("-86228.49000000954 278.06125000119226 87930.5 -3.59085052738 38.31272374472534 0.0");
		Tnt hPowerFinish = hPowerStart;  // initialzies z
		{
			cout << "Hammer Power location without any booster: " << std::endl;
			Tnt hPowerCopy = hPowerStart;
			hPowerCopy.print();
			hPowerCopy.freefall(1);
			hPowerCopy.print();
			cout << endl;
		}


		double* hpBoostersArr = nullptr;
		double* hpBoostersArrY = nullptr;
		{
			std::lock_guard<std::mutex> lock(mtx);
			// X-gains for the arrow power's 12 boosters
			hpBoostersArr = hPowerStart.gainArray(hpBoosters, numOfHammerPowerAdj, 1.0, direction::xDir);
			for (int i = 0; i < numOfHammerPowerAdj; i++) {
				std::cout << "hexposureX[" << i << "]: " << hpBoostersArr[i] << std::endl;
			} std::cout << endl;
			// Y-gains for the hammer power's 12 boosters, temporary, just testing
			hpBoostersArrY = hPowerStart.gainArray(hpBoosters, numOfHammerPowerAdj, 1.0, direction::yDir);
			for (int i = 0; i < numOfHammerPowerAdj; i++) {
				std::cout << "hexposureY[" << i << "]: " << hpBoostersArrY[i] << std::endl;
			} std::cout << endl;
		}
		

		// Spawn Hammer
		Tnt hammerStart("-86214.5 316.3952237679115 87930.5 0.0 -0.37519614825445 0.0"); // vels are: -2.7200464103316335E-17 -0.40689222528933305 0.0
		hammerStart.amount = 1440;
		Tnt hammerFinish = hammerStart;

		// Initializing stuff
		long int counter = 0;
		int usableCounter = 0;
		int dispenserSectionSize = 36;
		int dispenserSections = 36;
		double hPowerLocX; // powers horizontal coordinate (this is probably just going to be x)
		double hPowerLocMotStart = hPowerStart.getX() + hPowerStart.getU(); // saves on an addition
		double hPowerLocY = hPowerStart.getY() + hPowerStart.getV() + Tnt::gravity; // y coordinate of the powers
		
		hPowerFinish.setY(hPowerLocY);
		hPowerFinish.amount = 1296;
		// stuff to save/output
		unsigned int index;
		unsigned long int swaps = 0;
		double currentScore, existingScore;

		// Initializing ratios
		vector<int> startVector = { param, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		int max = 6;
		vector<int> finalVector = { param, max, max, max, max, max, max, max, max, max, max, max + 1 };
		vector<int> currentVector = startVector;
		{
			std::lock_guard<std::mutex> lock(mtx);
			std::cout << "id " << id << " starts at " << currentVector << " full: " << numOfPopulatedElements(hammerRatios) << "/" << maxRangeIndex << " swaps: " << swaps << endl;
		}
		do {

			if (currentVector[11] > max) { continue; }
			counter++;
			if (counter % 1000000000 == 0) {
				{
					std::lock_guard<std::mutex> lock(mtx);
					std::cout << "id " << id << " at " << currentVector << " full: " << numOfPopulatedElements(hammerRatios) << "/" << maxRangeIndex << " swaps: " << swaps << endl;
				}
			}
			hPowerLocX = hPowerLocMotStart
				+ hpBoostersArr[0] * currentVector[0]
				+ hpBoostersArr[1] * currentVector[1]
				+ hpBoostersArr[2] * currentVector[2]
				+ hpBoostersArr[3] * currentVector[3]
				+ hpBoostersArr[4] * currentVector[4]
				+ hpBoostersArr[5] * currentVector[5]
				+ hpBoostersArr[6] * currentVector[6]
				+ hpBoostersArr[7] * currentVector[7]
				+ hpBoostersArr[8] * currentVector[8]
				+ hpBoostersArr[9] * currentVector[9]
				+ hpBoostersArr[10] * currentVector[10]
				+ hpBoostersArr[11] * currentVector[11];
			if (hPowerLocX > centerX - 0.05018546 || hPowerLocX - centerX < -7.9) {// 20000.01vel // 0.88 was good 0.05 for max range // bounds of this may change due to tnt higher drag
				continue;
			}
			hPowerFinish.setX(hPowerLocX);
			hammerFinish = hammerStart;
			
			hammerFinish.explosion(&hPowerFinish);
			hammerFinish.freefall(17);
			// compute the index
			index = (int)round(hammerFinish.getX() - centerX); // applied velocity to corresponding index
			if (index >= maxRangeIndex||index < 252) { continue; } // if too much range
			if (arrowRatios[index].empty()) { continue; }

			
			Arrow affectedArrow(arrowRatios[index][13], 305.30882660809811568, hammerFinish.getZ(), arrowRatios[index][12], -0.04298819155656364, 0.0);
			affectedArrow.explosion(&hammerFinish);
			// Create new entry
			
			if (hammerRatios[index].empty()) { // populate elements that are empty

				{
					std::lock_guard<std::mutex> lock(mtx);
					hammerRatios[index].resize(15);
					for (size_t i = 0; i < 12; i++) {
						hammerRatios[index][i] = static_cast<double>(currentVector[i]);
					}
					hammerRatios[index][12] = affectedArrow.getU();
					hammerRatios[index][13] = affectedArrow.getV();
					hammerRatios[index][14] = hammerFinish.getX();
					// std::cout << "New element added at " << index << ": " << currentVector << " full: " << numOfPopulatedElements(hammerRatios) << "/" << maxRangeIndex << " swaps: " << swaps << std::endl;
				}
				continue;
			}
			// Check for improvements
			existingScore = abs(hammerRatios[index][12]);
			currentScore = abs(affectedArrow.getU());
			if (currentScore < existingScore) {
				{
					std::lock_guard<std::mutex> lock(mtx);
					// display the swap locations
					/*std::cout << "swapping[" << index << "]: "
						<< ratios[index][12] * velToIndex - centerX << " -> "
						<< centerArrowFinish.getU() * velToIndex - centerX << std::endl;
					*/
					// No need to resize
					for (size_t i = 0; i < 12; i++) {
						hammerRatios[index][i] = static_cast<double>(currentVector[i]);
					}
					hammerRatios[index][12] = affectedArrow.getU();
					hammerRatios[index][13] = affectedArrow.getV();
					hammerRatios[index][14] = hammerFinish.getX();
				}
				swaps++;
			}
			else {
				continue;
			}

		} while (generateNextCombination(currentVector, finalVector) && currentVector != finalVector);

		delete[] hpBoostersArr;
		delete[] hpBoostersArrY;
		
	}
	{
		std::lock_guard<std::mutex> lock(mtx);
		std::cout << "Finished [" << id << "] parameters: " << parameters << std::endl;
	}
	return;
}

int main() {
    cout << setprecision(17);
    cout << setfill(' ');
    //cout << fixed;

	// Initialize mutex and outputRatios for arrows and hammer
	mutex mtx;
    vector<thread> threads;
	const int maxRangeIndex = 1001;// 9001
	vector<vector<double>> arrowRatios(maxRangeIndex); // index actually corresponds to the range, index 20,000 is 19,999 blocks
	vector<vector<double>> hammerRatios(maxRangeIndex);
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
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        threads.emplace_back(runSimulationArrows, parameters[i], ref(arrowRatios), ref(mtx), i);
    }
	// Join all threads
	for (auto& thread : threads) {
		if (thread.joinable()) {
			thread.join();
		}
	}

	// Output the ratios found for ARROWS
	for (int i = 0; i < arrowRatios.size(); i++) {
		cout <<  i << " " << arrowRatios[i];
		if (!arrowRatios[i].empty()) {
			cout << "\tpositions: " << arrowRatios[i][12] * velToIndex + centerX << std::endl;
			// swap booster location info into final position

			arrowRatios[i][13] = arrowRatios[i][12] * velToIndex + centerX; // replace arrows power location with final position
			arrowRatios[i][12] = pow(Arrow::drag, 17) * arrowRatios[i][12];	// replace intial velocity with final velocityV
		}
		else
			std::cout << std::endl;
	}

	// Small pause
	std::this_thread::sleep_for(std::chrono::milliseconds(4000));
	// Clear or reset variables for hammer ratios
	threads.clear();
	parameters.clear();
	parameters.resize(numThreads);

	// Setup new or same parameters
	begin = 0;
	end = 12;
	counter = 0;
	std::cout << "do we retunr here" << std::endl;
	for (int i = begin; i <= end; i++) {
		parameters[counter].push_back(i);
		counter++;
		if (counter == numThreads) { counter = 0; }
	}
	for (unsigned int i = 0; i < numThreads; ++i) {
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
		threads.emplace_back(runSimulationHammer, parameters[i], ref(arrowRatios), ref(hammerRatios), ref(mtx), i);
	}
	// Join all threads
	for (auto& thread : threads) {
		if (thread.joinable()) {
			thread.join();
		}
	}
	// Output the ratios found for HAMMER AND ARROWS
	for (int i = 0; i < hammerRatios.size(); i++) {
		cout << "Hammer " << i << ": " << hammerRatios[i];
		if (!hammerRatios[i].empty()&& !arrowRatios[i].empty()) {
			cout << "\tHposX: " << hammerRatios[i][14] << std::endl;
			cout << "\tArrow " << arrowRatios[i] << "\tAposX: " << arrowRatios[i][13] << "\tarrowU: " << hammerRatios[i][12] << "\tarrowV: " << hammerRatios[i][13] << std::endl;
		}
		else
			std::cout << std::endl;
	}
	std::cout << "DONE?" << std::endl;


    return 0;
}
