#pragma once
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include <algorithm>

#include <vector>

using namespace std;

vector<vector<vector<long>*>> from_file() {
	vector<vector<vector<long>*>> vals(101, vector<vector<long>*>(101, nullptr));
	long p = 0;
	ifstream infile("z4_table.txt");
	string line;
	while (getline(infile, line)) {
		istringstream iss(line);
		long c;
		iss >> c;
		if (c == -1) {
			continue;
		}
		long set_head;
		iss >> set_head;
		if (set_head == -1) {
			p = c;
			continue;
		}
		vector<long>* tmp = new vector<long>();
		do {
			tmp->push_back(set_head);
		} while (iss >> set_head);
		sort(tmp->begin(), tmp->end());
		vals[p - 2][c] = tmp;
	}
	cout << "done processing file!" << endl;
	return vals;
}

struct Lookup_Table {
	long lt[50][100][8];
} lookup_table;

__device__ long fmt_p(long p) {
	if (p == 2) {
		return 0;
	} else {
		return (p - 1)/2;
	}
}

Lookup_Table raw_lt() {
	auto a = from_file();
	Lookup_Table ret;
	// zero initialize table
	for (long x = 0; x < 50; x++) {
		for (long y = 0; y < 100; y++) {
			for (long z = 0; z < 8; z++) {
				ret.lt[x][y][z] = 0;
			}
		}
	}
	// account for p=2
	if (a[0][0] != nullptr) {
		for (long s = 0; s < 2; s++) {
			for (long i = 0; i < 8; i++) {
				ret.lt[0][s][i] = (*a[0][s])[i];
			}
		}
	}
	for (long p = 3; p < 100; p += 2) {
		if (a[p-2][0] != nullptr) {
			for (long s = 0; s < p; s++) {
				long i = 0;
				for (long val : *a[p-2][s]) {
					ret.lt[(p-1)/2][s][i] = val;
					i++;
				}
			}
		}
	}
	return ret;
}


vector<vector<vector<long>*>> from_file_z3() {
	vector<vector<vector<long>*>> vals(101, vector<vector<long>*>(101, nullptr));
	long p = 0;
	ifstream infile("z3_table.txt");
	string line;
	while (getline(infile, line)) {
		istringstream iss(line);
		long c;
		iss >> c;
		if (c == -1) {
			continue;
		}
		long set_head;
		iss >> set_head;
		if (set_head == -1) {
			p = c;
			continue;
		}
		vector<long>* tmp = new vector<long>();
		do {
			tmp->push_back(set_head);
		} while (iss >> set_head);
		sort(tmp->begin(), tmp->end());
		vals[p - 2][c] = tmp;
	}
	cout << "done processing file!" << endl;
	return vals;
}


Lookup_Table raw_lt_z3() {
	auto a = from_file_z3();
	Lookup_Table ret;
	// zero initialize table
	for (long x = 0; x < 50; x++) {
		for (long y = 0; y < 100; y++) {
			for (long z = 0; z < 8; z++) {
				ret.lt[x][y][z] = 0;
			}
		}
	}
	// account for p=2
	if (a[0][0] != nullptr) {
		for (long s = 0; s < 2; s++) {
			for (long i = 0; i < 8; i++) {
				ret.lt[0][s][i] = (*a[0][s])[i];
			}
		}
	}
	for (long p = 3; p < 100; p += 2) {
		if (a[p-2][0] != nullptr) {
			for (long s = 0; s < p; s++) {
				long i = 0;
				for (long val : *a[p-2][s]) {
					ret.lt[(p-1)/2][s][i] = val;
					i++;
				}
			}
		}
	}
	return ret;
}

__device__
long* first_zero(long* inp) {
	long* tmp = inp;
	while (*tmp != 0) {
		tmp++;
	}
	return tmp;
}
