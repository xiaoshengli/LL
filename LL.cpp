//Copyright (c) 2019 Xiaosheng Li (xli22@gmu.edu)
//Reference: Linear Time Motif Discovery in Time Series, SDM 2019
/*
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <math.h>
#include <functional>
#include <random>
#include <algorithm>

using namespace std;

#define precision 1e-5
#define INF 1e20

long inline fast_floor(const double x) { return x > 0 ? (long)x : (long)x - 1; }

//boost hash_combine
template <class T>
inline void hash_combine(std::size_t & seed, const T & v) {
	std::hash<T> hasher;
	seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

typedef unordered_set<int> myset;
typedef unordered_map<size_t, vector<int> > mymap;

vector<int> transformIndex(vector<int> &motifs, const int* index) {
	vector<int>::iterator it;
	vector<int> newMotifs;
	for (it = motifs.begin(); it != motifs.end(); ++it) {
		newMotifs.push_back(*(index + *it));
	}
	return newMotifs;
}

void printMotifs(vector<int> &motifs, const int l, const double* TS) {
	for (size_t i = 0; i < motifs.size(); i += 2) {
		if (motifs[i] < motifs[i + 1])
			cout << motifs[i] << ", " << motifs[i + 1] << endl;
		else
			cout << motifs[i+1] << ", " << motifs[i] << endl;
	}
}

double dist(const double* TS, const int j, const int k, const int l, const double* Tmu, const double* Tsigma){
	if (abs(j - k)<l) return INF;
	int i;
	double d = 0;
	const double mu1 = *(Tmu + j);
	const double mu2 = *(Tmu + k);
	const double sigma1 = *(Tsigma + j);
	const double sigma2 = *(Tsigma + k);
	const double* t1 = TS + j;
	const double* t2 = TS + k;

	double maxd = 0;
	for (i = 0; i<l; ++i) {
		const double d1 = (*(t1 + i) - mu1) / sigma1;
		const double d2 = (*(t2 + i) - mu2) / sigma2;
		d = abs(d1 - d2);
		if (maxd < d) maxd = d;
	}

	return maxd;

}

double dist2(const double* TS, const int j, const int k, const int l, const double* Tmu, const double* Tsigma, const double* currentT, const double min){
	if (abs(j - k)<l) return INF;
	int i;
	double d = 0;
	const double mu1 = *(Tmu + j);
	const double sigma1 = *(Tsigma + j);
	const double* t1 = TS + j;

	double maxd = 0;
	for (i = 0; i<l; ++i) {
		const double d2 = *(currentT + i);
		const double d1 = (*(t1 + i) - mu1) / sigma1;
		d = abs(d1 - d2);
		if (maxd < d) maxd = d;
		if (maxd > min) return min + 1;
	}
	return maxd;
}

void retrieveGrids(vector<mymap> &grids, myset &neighbors, const size_t* currentV){
	int g;
	const int p = grids.size();

	for (g = 0; g<p; ++g) {
		if (grids[g].find(*(currentV + g)) != grids[g].end()){
			neighbors.insert(grids[g][*(currentV + g)].begin(), grids[g][*(currentV + g)].end());
		}
	}
}

void insertGrids(vector<mymap> &grids, const int j, const size_t* currentV){
	int g;
	const int p = grids.size();

	for (g = 0; g<p; ++g) {
		if (grids[g].find(*(currentV + g)) == grids[g].end()){
			vector<int> subs;
			subs.push_back(j);
			grids[g].insert(make_pair(*(currentV + g), subs));
		}
		else {
			grids[g][*(currentV + g)].push_back(j);
		}
	}
}

double minDist(myset &neighbors, const double* TS, const int j, vector<int> &motifs, const double delta, const int l, const double* Tmu, const double* Tsigma, const double* currentT, const int* index) {
	myset::const_iterator it;
	double min = delta;
	double d;
	const int i = *(index + j);
	for (it = neighbors.cbegin(); it != neighbors.cend(); ++it) {
		d = dist2(TS, *(index + *it), i, l, Tmu, Tsigma, currentT, min);
		if (d < precision) d = precision;
		if (d <= min) {
			if (d < min) {
				motifs.clear();
				min = d;
			}
			motifs.push_back(j);
			motifs.push_back(*it);
		}
	}
	return min;
}

void statistics(const double* TS, double* Tmu, double* Tsigma, const int n, const int l){
	int i;
	double sum = 0;
	double sum2 = 0;
	double mu;
	const double* pt = TS;
	double* pcx = Tmu;
	double* pcx2 = Tsigma;
	double d, d2;
	for (i = 0; i < l; ++i) {
		d = *pt;
		sum += d;
		sum2 += d*d;
		pt++;
	}
	mu = sum / l;
	*(pcx++) = mu;
	*(pcx2++) = sqrt(sum2 / l - mu*mu);
	for (i = 1; i < n - l + 1; ++i) {
		d = *pt;
		d2 = *(pt - l);
		sum = sum + d - d2;
		sum2 = sum2 + d*d - d2*d2;
		mu = sum / l;
		*(pcx++) = mu;
		*(pcx2++) = sqrt(sum2 / l - mu*mu);
		pt++;
	}
}

void update(const double* TS, double* currentT, const int j, const int l, const double* Tmu, const double* Tsigma){
	int i;
	const double mu = *(Tmu + j);
	const double sigma = *(Tsigma + j);
	const double* t = TS + j;
	for (i = 0; i<l; ++i) {
		*(currentT + i) = (*(t + i) - mu) / sigma;
	}
}

void updateV(size_t* currentV, const int l, const int p, const double delta, const double* currentT){
	int i, g;
	long d;
	const double pdelta = p*delta;
	size_t temp;

	for (g = 0; g<p; ++g) {
		temp = 0;
		for (i = 0; i<l; ++i) {
			d = fast_floor((*(currentT + i) - g*delta) / pdelta);
			hash_combine(temp, d);
		}
		*(currentV + g) = temp;
	}
}

void buildGrids(vector<mymap> &grids, const double* TS, const int i, const double delta, const int l, const double* Tmu, const double* Tsigma, double* currentT, const int* index){
	int g;
	const int p = grids.size();
	for (g = 0; g<p; ++g) {
		grids[g].clear();
	}
	const double pdelta = p*delta;

	int j, k;
	long d;
	size_t temp;
	int q;

	for (j = 0; j < i; ++j) {
		q = *(index + j);
		update(TS, currentT, q, l, Tmu, Tsigma);
		for (g = 0; g<p; ++g) {
			temp = 0;
			for (k = 0; k<l; ++k) {
				d = fast_floor((*(currentT + k) - g*delta) / pdelta);
				hash_combine(temp, d);
			}

			if (grids[g].find(temp) == grids[g].end()){
				vector<int> subs;
				subs.push_back(j);
				grids[g].insert(make_pair(temp, subs));
			}
			else {
				grids[g][temp].push_back(j);
			}
		}
	}
}

inline void shuffleIndex(int* index, const int n, const int l){
	int i;
	for (i = 0; i < n - l + 1; ++i) {
		*(index + i) = i;
	}
	random_shuffle(index, index + n - l + 1);
}

inline int firstPair(const int* index, const int l) {
	int i = 1;
	const int j = *index;
	while (abs(*(index + i) - j) < l){
		i++;
	}
	return i;
}

void loadData(char* fileName, double* TS, const int n){
	FILE * fp = fopen(fileName, "r");
	double d;
	int i = 0;

	if (fp == NULL) {
		cout << "Error!! Input file " << fileName << " not found." << endl;
		exit(0);
	}

	while (fscanf(fp, "%lf", &d) != EOF && i < n){
		*(TS + i) = d;
		i++;
	}
	fclose(fp);
}

int main(int argc, char *argv[]) {

	myset neighbors;
	int i;
	int j;

	double d, delta;

	const int n = atoi(argv[2]);
	const int l = atoi(argv[3]);
	const int g1 = atoi(argv[4]);

	int g = g1;

	vector<mymap> grids(g);

	double tStart1, tEnd1;

	double* TS = (double *)malloc(sizeof(double) * n);
	int* index = (int *)malloc(sizeof(int) * (n - l + 1));

	srand((int)time(NULL));
	shuffleIndex(index, n, l);

	loadData(argv[1], TS, n);

	tStart1 = clock();

	double* Tmu = (double *)malloc(sizeof(double) * (n - l + 1));
	double* Tsigma = (double *)malloc(sizeof(double) * (n - l + 1));
	statistics(TS, Tmu, Tsigma, n, l);

	double* currentT = (double *)malloc(sizeof(double) * l);

	int fp = firstPair(index, l);

	delta = dist(TS, *index, *(index + fp), l, Tmu, Tsigma);
	buildGrids(grids, TS, fp + 1, delta, l, Tmu, Tsigma, currentT, index);

	vector<int> motifs;
	motifs.push_back(0);
	motifs.push_back(fp);
	if (delta < precision) delta = precision;
	size_t* currentV = (size_t *)malloc(sizeof(size_t) * g);

	for (i = fp + 1; i < n - l + 1; ++i) {

		j = *(index + i);
		neighbors.clear();
		update(TS, currentT, j, l, Tmu, Tsigma);
		updateV(currentV, l, g, delta, currentT);

		retrieveGrids(grids, neighbors, currentV);
		d = minDist(neighbors, TS, i, motifs, delta, l, Tmu, Tsigma, currentT, index);

		if (d >= delta) {
			insertGrids(grids, i, currentV);
		}
		else {
			delta = d;
			buildGrids(grids, TS, i + 1, delta, l, Tmu, Tsigma, currentT, index);
		}
	}
	
	if (argc > 5) {
		const int m1 = motifs[0];
		const int m2 = motifs[1];

		int index0 = *(index);
		int index1 = *(index + 1);
		*(index) = *(index + m1);
		*(index + 1) = *(index + m2);
		*(index + m1) = index0;
		*(index + m2) = index1;
		g = atoi(argv[5]);
		fp = 1;
		vector<mymap> grids2(g);

		delta = dist(TS, *index, *(index + fp), l, Tmu, Tsigma);
		buildGrids(grids2, TS, fp + 1, delta, l, Tmu, Tsigma, currentT, index);

		motifs.clear();
		motifs.push_back(0);
		motifs.push_back(fp);
		if (delta < precision) delta = precision;

		size_t* currentV2 = (size_t *)malloc(sizeof(size_t) * g);

		for (i = fp + 1; i < n - l + 1; ++i) {
			j = *(index + i);
			neighbors.clear();
			update(TS, currentT, j, l, Tmu, Tsigma);
			updateV(currentV2, l, g, delta, currentT);

			retrieveGrids(grids2, neighbors, currentV2);
			d = minDist(neighbors, TS, i, motifs, delta, l, Tmu, Tsigma, currentT, index);

			if (d >= delta) {
				insertGrids(grids2, i, currentV2);
			}
			else {
				delta = d;
				buildGrids(grids2, TS, i + 1, delta, l, Tmu, Tsigma, currentT, index);
			}
		}
		free(currentV2);
	}
	
	
	if (argc > 6) {
		const int mm1 = motifs[0];
		const int mm2 = motifs[1];

		int index0 = *(index);
		int index1 = *(index + 1);
		*(index) = *(index + mm1);
		*(index + 1) = *(index + mm2);
		*(index + mm1) = index0;
		*(index + mm2) = index1;

		g = atoi(argv[6]);
		fp = 1;
		vector<mymap> grids3(g);
		delta = dist(TS, *index, *(index + fp), l, Tmu, Tsigma);
		buildGrids(grids3, TS, fp + 1, delta, l, Tmu, Tsigma, currentT, index);

		motifs.clear();
		motifs.push_back(0);
		motifs.push_back(fp);
		if (delta < precision) delta = precision;

		size_t* currentV3 = (size_t *)malloc(sizeof(size_t) * g);

		for (i = fp + 1; i < n - l + 1; ++i) {
			j = *(index + i);
			neighbors.clear();
			update(TS, currentT, j, l, Tmu, Tsigma);
			updateV(currentV3, l, g, delta, currentT);

			retrieveGrids(grids3, neighbors, currentV3);
			d = minDist(neighbors, TS, i, motifs, delta, l, Tmu, Tsigma, currentT, index);

			if (d >= delta) {
				insertGrids(grids3, i, currentV3);
			}
			else {
				delta = d;
				buildGrids(grids3, TS, i + 1, delta, l, Tmu, Tsigma, currentT, index);
			}
		}
		free(currentV3);
	}

	if (delta == precision) delta = 0;

	tEnd1 = clock();

	cout << "time series length: " << n << " motif length:" << l << endl;
	cout << "motif distance(Chebyshev): " << delta << endl;
	cout << "motif: " << endl;
	vector<int> newMotifs = transformIndex(motifs, index);
	printMotifs(newMotifs, l, TS);
	cout << "The running time is: " << (tEnd1 - tStart1) / CLOCKS_PER_SEC << " seconds" << endl;

	free(TS);
	free(index);
	free(Tmu);
	free(Tsigma);
	free(currentT);
	free(currentV);	

	return 0;
}