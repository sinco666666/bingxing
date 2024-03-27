//////#include<iostream>
//////#include <windows.h>
//////#include <stdlib.h>
//////using namespace std;
//////const int N = 10240;
//////void initialize(float*a, float**b,int n)
//////{
//////	for (int i = 0; i < n; i++)
//////		a[i] = i + 10;
//////	for (int i = 0; i < n; i++)
//////		for (int j = 0; j < n; j++)
//////			b[i][j] = i + j;
//////}
//////void Trivial_algorithm(float *a,float **b,int n)
//////{
//////	int i = 0;
//////	int j = 0;
//////	float* sum;
//////	sum = new float[n];
//////	long long head, tail, freq;
//////	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//////	QueryPerformanceCounter((LARGE_INTEGER*)&head);
//////	for (i = 0; i < n; i++)
//////	{
//////		sum[i] = 0.0;
//////		for (j = 0; j < n; j++)
//////			sum[i] += b[j][i] * a[j];
//////	}
//////	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//////	cout << "Trivial_algorithm" << (tail - head) * 1000.0 / freq << "ms" << endl;
//////}
//////
//////void cache_algorithm(float* a, float** b, int n)
//////{
//////	int i = 0;
//////	int j = 0;
//////	float* sum;
//////	sum = new float[n];
//////	long long head, tail, freq;
//////	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//////	QueryPerformanceCounter((LARGE_INTEGER*)&head);
//////	for (i = 0; i < n; i++)
//////		sum[i] = 0.0;
//////	for (j = 0; j < n; j++)
//////		for (i = 0; i < n; i++)
//////			sum[i] += b[j][i] * a[j];
//////	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//////	cout << "cache_algorithm" << (tail - head) * 1000.0 / freq << "ms" << endl;
//////}
//////int main()
//////{
//////	int n = 10000;
//////	float* a;
//////	a = new float[n];
//////	float** b = new float*[n];
//////	for (int i = 0; i < n; i++)
//////		b[i] = new float[n];
//////	initialize(a, b, n);
//////	Trivial_algorithm(a, b, n);
//////	cache_algorithm(a, b, n);
//////	
//////}
////
//#include<iostream>
//#include<windows.h>
//#include<stdlib.h>
//using namespace std;
//
//const int N = 10240;
//const int TEST_TIMES = 10; // 定义测试次数
//
//void initialize(float* a, float** b, int n) {
//	for (int i = 0; i < n; i++)
//		a[i] = i + 10;
//	for (int i = 0; i < n; i++)
//		for (int j = 0; j < n; j++)
//			b[i][j] = i + j;
//}
//
//// 函数指针类型定义
//typedef void(*TestFunction)(float*, float**, int);
//
//void Trivial_algorithm(float* a, float** b, int n) {
//	float* sum = new float[n];
//	for (int i = 0; i < n; i++) {
//		sum[i] = 0.0;
//		for (int j = 0; j < n; j++)
//			sum[i] += b[j][i] * a[j];
//	}
//	delete[] sum; // 释放内存
//}
//
//void cache_algorithm(float* a, float** b, int n) {
//	float* sum = new float[n];
//	for (int i = 0; i < n; i++)
//		sum[i] = 0.0;
//	for (int j = 0; j < n; j++)
//		for (int i = 0; i < n; i++)
//			sum[i] += b[j][i] * a[j];
//	delete[] sum; // 释放内存
//}
//
//void unroll_algorithm(float* a, float** b, int n) {
//	float* sum = new float[n];
//	for (int i = 0; i < n; i++) {
//		sum[i] = 0.0;
//		// 主体循环，每次迭代处理10个元素
//		for (int j = 0; j < (n / 10) * 10; j += 10) {
//			sum[i] += b[j][i] * a[j];
//			sum[i] += b[j + 1][i] * a[j + 1];
//			sum[i] += b[j + 2][i] * a[j + 2];
//			sum[i] += b[j + 3][i] * a[j + 3];
//			sum[i] += b[j + 4][i] * a[j + 4];
//			sum[i] += b[j + 5][i] * a[j + 5];
//			sum[i] += b[j + 6][i] * a[j + 6];
//			sum[i] += b[j + 7][i] * a[j + 7];
//			sum[i] += b[j + 8][i] * a[j + 8];
//			sum[i] += b[j + 9][i] * a[j + 9];
//		}
//		// 处理因循环展开而剩余的元素
//		for (int j = (n / 10) * 10; j < n; j++) {
//			sum[i] += b[j][i] * a[j];
//		}
//	}
//	delete[] sum; // 释放内存
//}
//
//
//
//void run_and_measure(TestFunction func, float* a, float** b, int n, const string& funcName) {
//	long long totalTime = 0;
//	long long head, tail, freq;
//	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//
//	QueryPerformanceCounter((LARGE_INTEGER*)&head);
//	for (int i = 0; i < TEST_TIMES; i++)
//		func(a, b, n);
//	QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//	totalTime += tail - head;
//
//	double averageTime = (double)totalTime * 1000.0 / (freq * TEST_TIMES);
//	//cout << funcName << " average time: " << averageTime << "ms" << endl;
//	//cout << funcName << " with n = " << n << ", average time: " << averageTime << "ms" << endl;
//	cout <<averageTime<< endl;
//}
//
//int main() {
//	////int n = 10; // 使用N作为矩阵的大小
//	//for (int n = 1000; n <= 30000; n += 1000) {
//	//	float* a = new float[n];
//	//	float** b = new float* [n];
//	//	for (int i = 0; i < n; i++)
//	//		b[i] = new float[n];
//
//	//	initialize(a, b, n);
//
//	//	//run_and_measure(Trivial_algorithm, a, b, n, "Trivial_algorithm");
//	//	run_and_measure(cache_algorithm, a, b, n, "Cache_algorithm");
//	//	//run_and_measure(unroll_algorithm, a, b, n, "Unroll_algorithm");
//	//	//cout << endl;
//
//	//	// 清理资源s
//	//	delete[] a;
//	//	for (int i = 0; i < n; i++)
//	//		delete[] b[i];
//	//	delete[] b;
//	//}
//	//int n = 10000;
//	//float* a = new float[n];
//	//float** b = new float* [n];
//	//for (int i = 0; i < n; i++)
//	//	b[i] = new float[n];
//
//	//initialize(a, b, n);
//
//	////run_and_measure(Trivial_algorithm, a, b, n, "Trivial_algorithm");
//	////run_and_measure(cache_algorithm, a, b, n, "Cache_algorithm");
//	//run_and_measure(unroll_algorithm, a, b, n, "Unroll_algorithm");
//	////cout << endl;
//
//	//// 清理资源
//	//delete[] a;
//	//for (int i = 0; i < n; i++)
//	//	delete[] b[i];
//	//delete[] b;
//
//	return 0;
//}
//
