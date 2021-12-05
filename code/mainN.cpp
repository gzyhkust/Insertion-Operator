#include "insertionN.h"

int main(int argc, char **args) {

	readInput();
	initGrid();

	cout << "Start sharing..." << endl;
	clock_t tbegin = clock();
  	//timeDependentInsertion();
	test();
	clock_t tend = clock();
	cout << "Finish sharing \t time cost: " + to_string((tend - tbegin) / CLOCKS_PER_SEC) + "s" << endl;
	
	// ofstream of("./data/sharingn.txt");
	// of << (tend - tbegin) / CLOCKS_PER_SEC << "\n";
	// recordTrajectory();
	// freeMemory();
}