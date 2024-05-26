
//////////////////////////////////////////////////////低效的程序框架
//# include <iostream>
//# include <pthread.h>
//#include <windows.h>
//
//#include <pmmintrin.h>
//#include <xmmintrin.h>
//#include <emmintrin.h>
//#include <smmintrin.h>
//#include <tmmintrin.h>
//#include <nmmintrin.h>
//#include <immintrin.h> 、
//using namespace std;
//
//const int n = 100;
//float A[n][n];
//int worker_count = 7; 
//void init()
//{
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			A[i][j] = 0;
//		}
//		A[i][i] = 1.0;
//		for (int j = i + 1; j < n; j++)
//			A[i][j] = rand() % 100;
//	}
//
//	for (int i = 0; i < n; i++)
//	{
//		int k1 = rand() % n;
//		int k2 = rand() % n;
//		for (int j = 0; j < n; j++)
//		{
//			A[i][j] += A[0][j];
//			A[k1][j] += A[k2][j];
//		}
//	}
//}
//
//void f_ordinary()
//{
//	for (int k = 0; k < n; k++)
//	{
//		for (int j = k + 1; j < n; j++)
//		{
//			A[k][j] = A[k][j] * 1.0 / A[k][k];
//		}
//		A[k][k] = 1.0;
//
//		for (int i = k + 1; i < n; i++)
//		{
//			for (int j = k + 1; j < n; j++)
//			{
//				A[i][j] = A[i][j] - A[i][k] * A[k][j];
//			}
//			A[i][k] = 0;
//		}
//	}
//}
//
//
//struct threadParam_t
//{
//	int k; //消去的轮次
//	int t_id; // 线程 id
//};
//
//void* threadFunc(void* param) {
//	threadParam_t* p = (threadParam_t*)param;
//	int k = p->k; //消去的轮次
//	int t_id = p->t_id; //线程编号
//	int i = k + t_id + 1; //获取自己的计算任务
//
//	for (int j = k + 1; j < n; ++j) {
//		A[i][j] = A[i][j] - A[i][k] * A[k][j];
//	}
//	A[i][k] = 0;
//	pthread_exit(NULL);
//	return 0;
//}
//
//
//int main() {
//	init();
//
//	long long counter; // 记录次数
//	double seconds;
//	long long head, tail, freq, noww;
//	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//	QueryPerformanceCounter((LARGE_INTEGER*)&head); //开始计时
//
//	for (int k = 0; k < n; ++k) {
//		//主线程做除法操作
//		for (int j = k + 1; j < n; j++) {
//			A[k][j] = A[k][j] / A[k][k];
//		}
//		A[k][k] = 1.0;
//
//		//创建工作线程，进行消去操作
//		int worker_count = n - 1 - k; //工作线程数量
//		pthread_t* handles = new pthread_t[worker_count]; // 创建对应的 Handle
//		threadParam_t* params = new threadParam_t[worker_count]; // 创建对应的线程数据结构
//
//		//分配任务
//		for (int t_id = 0; t_id < worker_count; t_id++) {
//			params[t_id].k = k;
//			params[t_id].t_id = t_id;
//			pthread_create(&handles[t_id], NULL, threadFunc, (void*)&params[t_id]);
//		}
//
//		//主线程挂起等待所有的工作线程完成此轮消去工作
//		for (int t_id = 0; t_id < worker_count; t_id++) {
//			pthread_join(handles[t_id], NULL);
//		}
//
//		delete[] handles;
//		delete[] params;
//	}
//
//	QueryPerformanceCounter((LARGE_INTEGER*)&tail); //结束计时
//	seconds = (tail - head) * 1000.0 / freq; //单位 ms
//
//	cout << "pthread_dong: " << seconds << " ms" << endl;
//}

//////////////////////////////////////////////////////动态线程版本
# include <iostream>
# include <pthread.h>
#include <windows.h>

#include <pmmintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <tmmintrin.h>
#include <nmmintrin.h>
#include <immintrin.h> //AVX、AVX2
using namespace std;

const int n = 2000;
float A[n][n];
int worker_count = 7; //工作线程数量
void init()
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			A[i][j] = 0;
		}
		A[i][i] = 1.0;
		for (int j = i + 1; j < n; j++)
			A[i][j] = rand() % 100;
	}

	for (int i = 0; i < n; i++)
	{
		int k1 = rand() % n;
		int k2 = rand() % n;
		for (int j = 0; j < n; j++)
		{
			A[i][j] += A[0][j];
			A[k1][j] += A[k2][j];
		}
	}
}

void f_ordinary()
{
	for (int k = 0; k < n; k++)
	{
		for (int j = k + 1; j < n; j++)
		{
			A[k][j] = A[k][j] * 1.0 / A[k][k];
		}
		A[k][k] = 1.0;

		for (int i = k + 1; i < n; i++)
		{
			for (int j = k + 1; j < n; j++)
			{
				A[i][j] = A[i][j] - A[i][k] * A[k][j];
			}
			A[i][k] = 0;
		}
	}
}


struct threadParam_t
{
	int k; //消去的轮次
	int t_id; // 线程 id
};

void *threadFunc(void* param)
{

	__m256 va, vt, vx, vaij, vaik, vakj;

	threadParam_t* p = (threadParam_t*)param;
	int k = p->k; //消去的轮次
	int t_id = p->t_id; //线程编号
	int i = k + t_id + 1; //获取自己的计算任务
	for (int m = k + 1 + t_id; m < n; m += worker_count)
	{
		vaik = _mm256_set1_ps(A[k][k]);
		int j;
		for (j = k + 1; j + 8 <= n; j += 8)
		{
			vakj = _mm256_loadu_ps(&(A[k][j]));
			vaij = _mm256_loadu_ps(&(A[m][j]));
			vx = _mm256_mul_ps(vakj, vaik);
			vaij = _mm256_sub_ps(vaij, vx);

			_mm256_store_ps(&A[i][j], vaij);
		}
		for (; j < n; j++)
			A[m][j] = A[m][j] - A[m][k] * A[k][j];

		A[m][k] = 0;
	}


	pthread_exit(NULL);
	return 0;

}
void printMatrix() {
	cout << "Matrix after Gaussian Elimination:" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
}

int main()
{
	init();
	__m256 va2, vt2, vx2, vaij2, vaik2, vakj2;

	long long counter;// 记录次数
	double seconds;
	long long head, tail, freq, noww;
	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
	QueryPerformanceCounter((LARGE_INTEGER*)&head);//开始计时


	for (int k = 0; k < n; k++)
	{
		vt2 = _mm256_set1_ps(A[k][k]);
		int j;
		for (j = k + 1; j + 8 <= n; j += 8)
		{
			va2 = _mm256_loadu_ps(&(A[k][j]));
			va2 = _mm256_div_ps(va2, vt2);
			_mm256_store_ps(&(A[k][j]), va2);
		}

		for (; j < n; j++)
		{
			A[k][j] = A[k][j] * 1.0 / A[k][k];

		}
		A[k][k] = 1.0;

		//创建工作线程，进行消去操作

		pthread_t* handles = new pthread_t[worker_count];// 创建对应的 Handle
		threadParam_t* param = new threadParam_t[worker_count];// 创建对应的线程数据结构

		//分配任务
		for (int t_id = 0; t_id < worker_count; t_id++)
		{
			param[t_id].k = k;
			param[t_id].t_id = t_id;
		}
		//创建线程
		for (int t_id = 0; t_id < worker_count; t_id++)
			pthread_create(&handles[t_id], NULL, threadFunc, (void*)&param[t_id]);

		//主线程挂起等待所有的工作线程完成此轮消去工作
		for (int t_id = 0; t_id < worker_count; t_id++)
			pthread_join(handles[t_id], NULL);

	}


	QueryPerformanceCounter((LARGE_INTEGER*)&tail);//结束计时
	seconds = (tail - head) * 1000.0 / freq;//单位 ms

	cout << "pthread_dong: " << seconds << " ms" << endl;

}


/////////////////////////////////////////////动态测试代码
//#include <iostream>
//#include <pthread.h>
//#include <windows.h>
//
//#include <pmmintrin.h>
//#include <xmmintrin.h>
//#include <emmintrin.h>
//#include <smmintrin.h>
//#include <tmmintrin.h>
//#include <nmmintrin.h>
//#include <immintrin.h> //AVX、AVX2
//
//using namespace std;
//
//const int max_n = 6000; // 最大矩阵规模
//const int step = 256;   // 矩阵规模递增步长
//const int worker_count = 3; // 工作线程数量
//int n;
//
//float A[max_n][max_n]; // 全局数组 A，矩阵规模最大为 max_n
//
//void init(int n)
//{
//    for (int i = 0; i < n; i++)
//    {
//        for (int j = 0; j < n; j++)
//        {
//            A[i][j] = 0;
//        }
//        A[i][i] = 1.0;
//        for (int j = i + 1; j < n; j++)
//            A[i][j] = rand() % 100;
//    }
//
//    for (int i = 0; i < n; i++)
//    {
//        int k1 = rand() % n;
//        int k2 = rand() % n;
//        for (int j = 0; j < n; j++)
//        {
//            A[i][j] += A[0][j];
//            A[k1][j] += A[k2][j];
//        }
//    }
//}
//
//void f_ordinary(int n)
//{
//    for (int k = 0; k < n; k++)
//    {
//        for (int j = k + 1; j < n; j++)
//        {
//            A[k][j] = A[k][j] * 1.0 / A[k][k];
//        }
//        A[k][k] = 1.0;
//
//        for (int i = k + 1; i < n; i++)
//        {
//            for (int j = k + 1; j < n; j++)
//            {
//                A[i][j] = A[i][j] - A[i][k] * A[k][j];
//            }
//            A[i][k] = 0;
//        }
//    }
//}
//
//struct threadParam_t
//{
//    int k;     // 消去的轮次
//    int t_id;  // 线程 id
//};
//
//void* threadFunc(void* param)
//{
//    __m256 va, vt, vx, vaij, vaik, vakj;
//
//    threadParam_t* p = (threadParam_t*)param;
//    int k = p->k;       // 消去的轮次
//    int t_id = p->t_id; // 线程编号
//    int i = k + t_id + 1; // 获取自己的计算任务
//    for (int m = k + 1 + t_id; m < n; m += worker_count)
//    {
//        vaik = _mm256_set1_ps(A[k][k]);
//        int j;
//        for (j = k + 1; j + 8 <= n; j += 8)
//        {
//            vakj = _mm256_loadu_ps(&(A[k][j]));
//            vaij = _mm256_loadu_ps(&(A[m][j]));
//            vx = _mm256_mul_ps(vakj, vaik);
//            vaij = _mm256_sub_ps(vaij, vx);
//
//            _mm256_store_ps(&A[i][j], vaij);
//        }
//        for (; j < n; j++)
//            A[m][j] = A[m][j] - A[m][k] * A[k][j];
//
//        A[m][k] = 0;
//    }
//
//    pthread_exit(NULL);
//    return 0;
//}
//
//int main()
//{
//    for (int n = step; n <= max_n; n += step)
//    {
//        double parallel_time;
//
//        init(n);
//        long long head, tail, freq;
//        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//        QueryPerformanceCounter((LARGE_INTEGER*)&head); // 开始计时
//
//        // 进行并行高斯消元
//        for (int k = 0; k < n; k++)
//        {
//            __m256 vt2 = _mm256_set1_ps(A[k][k]);
//            int j;
//            for (j = k + 1; j + 8 <= n; j += 8)
//            {
//                __m256 va2 = _mm256_loadu_ps(&(A[k][j]));
//                va2 = _mm256_div_ps(va2, vt2);
//                _mm256_store_ps(&(A[k][j]), va2);
//            }
//
//            for (; j < n; j++)
//            {
//                A[k][j] = A[k][j] * 1.0 / A[k][k];
//            }
//            A[k][k] = 1.0;
//
//            // 创建工作线程，进行消去操作
//            pthread_t* handles = new pthread_t[worker_count]; // 创建对应的 Handle
//            threadParam_t* param = new threadParam_t[worker_count]; // 创建对应的线程数据结构
//
//            // 分配任务
//            for (int t_id = 0; t_id < worker_count; t_id++)
//            {
//                param[t_id].k = k;
//                param[t_id].t_id = t_id;
//            }
//            // 创建线程
//            for (int t_id = 0; t_id < worker_count; t_id++)
//                pthread_create(&handles[t_id], NULL, threadFunc, (void*)&param[t_id]);
//
//            // 主线程挂起等待所有的工作线程完成此轮消去工作
//            for (int t_id = 0; t_id < worker_count; t_id++)
//                pthread_join(handles[t_id], NULL);
//        }
//
//        QueryPerformanceCounter((LARGE_INTEGER*)&tail); // 结束计时
//        parallel_time = (tail - head) * 1000.0 / freq; // 单位 ms
//
//        cout << "Matrix size: " << n << "x" << n << endl;
//        cout << "pthread_dong: " << parallel_time << " ms" << endl;
//        cout << endl;
//
//        // 测试串行计算时间
//        init(n);
//        long long head1, tail1, freq1;
//        double seconds1;
//        QueryPerformanceFrequency((LARGE_INTEGER*)&freq1);
//        QueryPerformanceCounter((LARGE_INTEGER*)&head1); // 开始计时
//
//        f_ordinary(n);
//
//        QueryPerformanceCounter((LARGE_INTEGER*)&tail1); // 结束计时
//        seconds1 = (tail1 - head1) * 1000.0 / freq1; // 单位 ms
//
//        cout << "f_ordinary: " << seconds1 << " ms" << endl;
//
//        // 计算加速比
//        double speedup = seconds1 / parallel_time;
//        cout << "Speedup: " << speedup << endl;
//        cout << endl;
//    }
//
//    return 0;
//}





//////////////////////////////////////////动态
//#include <iostream>
//#include <pthread.h>
//#include <windows.h>
//
//#include <pmmintrin.h>
//#include <xmmintrin.h>
//#include <emmintrin.h>
//#include <smmintrin.h>
//#include <tmmintrin.h>
//#include <nmmintrin.h>
//#include <immintrin.h> //AVX、AVX2
//using namespace std;
//int n;
//const int max_n = 6000; // 最大矩阵规模
//const int step = 256;   // 矩阵规模递增步长
//int worker_count = 7; // 工作线程数量
//float A[max_n][max_n]; // 全局数组 A，矩阵规模最大为 max_n
//
//void init(int n)
//{
//    for (int i = 0; i < n; i++)
//    {
//        for (int j = 0; j < n; j++)
//        {
//            A[i][j] = 0;
//        }
//        A[i][i] = 1.0;
//        for (int j = i + 1; j < n; j++)
//            A[i][j] = rand() % 100;
//    }
//
//    for (int i = 0; i < n; i++)
//    {
//        int k1 = rand() % n;
//        int k2 = rand() % n;
//        for (int j = 0; j < n; j++)
//        {
//            A[i][j] += A[0][j];
//            A[k1][j] += A[k2][j];
//        }
//    }
//}
//
//// 其余函数（f_ordinary, threadParam_t, threadFunc）保持不变
//void f_ordinary()
//{
//    for (int k = 0; k < n; k++)
//    {
//        for (int j = k + 1; j < n; j++)
//        {
//            A[k][j] = A[k][j] * 1.0 / A[k][k];
//        }
//        A[k][k] = 1.0;
//
//        for (int i = k + 1; i < n; i++)
//        {
//            for (int j = k + 1; j < n; j++)
//            {
//                A[i][j] = A[i][j] - A[i][k] * A[k][j];
//            }
//            A[i][k] = 0;
//        }
//    }
//}
//
//
//struct threadParam_t
//{
//    int k; //消去的轮次
//    int t_id; // 线程 id
//};
//
//void* threadFunc(void* param)
//{
//
//    __m256 va, vt, vx, vaij, vaik, vakj;
//
//    threadParam_t* p = (threadParam_t*)param;
//    int k = p->k; //消去的轮次
//    int t_id = p->t_id; //线程编号
//    int i = k + t_id + 1; //获取自己的计算任务
//    for (int m = k + 1 + t_id; m < n; m += worker_count)
//    {
//        vaik = _mm256_set_ps(A[m][k], A[m][k], A[m][k], A[m][k], A[m][k], A[m][k], A[m][k], A[m][k]);
//        int j;
//        for (j = k + 1; j + 8 <= n; j += 8)
//        {
//            vakj = _mm256_loadu_ps(&(A[k][j]));
//            vaij = _mm256_loadu_ps(&(A[m][j]));
//            vx = _mm256_mul_ps(vakj, vaik);
//            vaij = _mm256_sub_ps(vaij, vx);
//
//            _mm256_store_ps(&A[i][j], vaij);
//        }
//        for (; j < n; j++)
//            A[m][j] = A[m][j] - A[m][k] * A[k][j];
//
//        A[m][k] = 0;
//    }
//
//
//    pthread_exit(NULL);
//    return 0;
//}
//
//int main()
//{
//    for (int n = step; n <= max_n; n += step)
//    {
//        init(n); // 使用当前矩阵大小进行初始化
//        __m256 va2, vt2, vx2, vaij2, vaik2, vakj2;
//
//        long long counter; // 记录次数
//        double seconds;
//        long long head, tail, freq, noww;
//        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//        QueryPerformanceCounter((LARGE_INTEGER*)&head); // 开始计时
//
//        for (int k = 0; k < n; k++)
//        {
//            vt2 = _mm256_set_ps(A[k][k], A[k][k], A[k][k], A[k][k], A[k][k], A[k][k], A[k][k], A[k][k]);
//            int j;
//            for (j = k + 1; j + 8 <= n; j += 8)
//            {
//                va2 = _mm256_loadu_ps(&(A[k][j]));
//                va2 = _mm256_div_ps(va2, vt2);
//                _mm256_store_ps(&(A[k][j]), va2);
//            }
//
//            for (; j < n; j++)
//            {
//                A[k][j] = A[k][j] * 1.0 / A[k][k];
//
//            }
//            A[k][k] = 1.0;
//
//            //创建工作线程，进行消去操作
//
//            pthread_t* handles = new pthread_t[worker_count];// 创建对应的 Handle
//            threadParam_t* param = new threadParam_t[worker_count];// 创建对应的线程数据结构
//
//            //分配任务
//            for (int t_id = 0; t_id < worker_count; t_id++)
//            {
//                param[t_id].k = k;
//                param[t_id].t_id = t_id;
//            }
//            //创建线程
//            for (int t_id = 0; t_id < worker_count; t_id++)
//                pthread_create(&handles[t_id], NULL, threadFunc, (void*)&param[t_id]);
//
//            //主线程挂起等待所有的工作线程完成此轮消去工作
//            for (int t_id = 0; t_id < worker_count; t_id++)
//                pthread_join(handles[t_id], NULL);
//
//        }
//
//        QueryPerformanceCounter((LARGE_INTEGER*)&tail); // 结束计时
//        seconds = (tail - head) * 1000.0 / freq; // 单位 ms
//
//        cout << "Matrix size: " << n << "x" << n << endl;
//        cout << "pthread_dong: " << seconds << " ms" << endl;
//        cout << endl;
//    }
//}



///////////////////////////////////////////静态线程 + 信号量同步版本
//#include <iostream>
//#include <pthread.h>
//#include <semaphore.h>
//#include <windows.h>
//
//#include <pmmintrin.h>
//#include <xmmintrin.h>
//#include <emmintrin.h>
//#include <smmintrin.h>
//#include <tmmintrin.h>
//#include <nmmintrin.h>
//#include <immintrin.h> //AVX、AVX2
//using namespace std;
//
//const int n = 2000;
//float A[n][n];
//const int NUM_THREADS = 5; // 工作线程数量
//
//void init()
//{
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			A[i][j] = 0;
//		}
//		A[i][i] = 1.0;
//		for (int j = i + 1; j < n; j++)
//			A[i][j] = rand() % 100;
//	}
//
//	for (int i = 0; i < n; i++)
//	{
//		int k1 = rand() % n;
//		int k2 = rand() % n;
//		for (int j = 0; j < n; j++)
//		{
//			A[i][j] += A[0][j];
//			A[k1][j] += A[k2][j];
//		}
//	}
//}
//
//void f_ordinary()
//{
//	for (int k = 0; k < n; k++)
//	{
//		for (int j = k + 1; j < n; j++)
//		{
//			A[k][j] = A[k][j] * 1.0 / A[k][k];
//		}
//		A[k][k] = 1.0;
//
//		for (int i = k + 1; i < n; i++)
//		{
//			for (int j = k + 1; j < n; j++)
//			{
//				A[i][j] = A[i][j] - A[i][k] * A[k][j];
//			}
//			A[i][k] = 0;
//		}
//	}
//}
//struct threadParam_t
//{
//    int k;     // 消去的轮次
//    int t_id;  // 线程 id
//};
//
//sem_t sem_main;
//sem_t sem_workerstart[NUM_THREADS];
//sem_t sem_workerend[NUM_THREADS];
//
//void* threadFunc(void* param)
//{
//    __m256 va, vt, vx, vaij, vaik, vakj;
//
//    threadParam_t* p = (threadParam_t*)param;
//    int k = p->k;
//    int t_id = p->t_id;
//    int i = k + t_id + 1;
//
//    for (int m = k + 1 + t_id; m < n; m += NUM_THREADS)
//    {
//        vaik = _mm256_set1_ps(A[k][k]);
//        int j;
//        for (j = k + 1; j + 8 <= n; j += 8)
//        {
//            vakj = _mm256_loadu_ps(&(A[k][j]));
//            vaij = _mm256_loadu_ps(&(A[m][j]));
//            vx = _mm256_mul_ps(vakj, vaik);
//            vaij = _mm256_sub_ps(vaij, vx);
//            _mm256_store_ps(&A[i][j], vaij);
//        }
//        for (; j < n; j++)
//            A[m][j] = A[m][j] - A[m][k] * A[k][j];
//
//        A[m][k] = 0;
//    }
//
//    sem_post(&sem_main); // 唤醒主线程
//    sem_wait(&sem_workerend[t_id]); // 等待主线程唤醒进入下一轮
//
//    pthread_exit(NULL);
//    return 0;
//}
//
//void printMatrix() {
//	cout << "Matrix after Gaussian Elimination:" << endl;
//	for (int i = 0; i < n; i++) {
//		for (int j = 0; j < n; j++) {
//			cout << A[i][j] << " ";
//		}
//		cout << endl;
//	}
//}
//
//int main()
//{
//    init();
//    long long counter;// 记录次数
//    double seconds;
//    long long head, tail, freq, noww;
//    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//    QueryPerformanceCounter((LARGE_INTEGER*)&head);//开始计时
//
//    sem_init(&sem_main, 0, 0);
//    for (int i = 0; i < NUM_THREADS; ++i)
//    {
//        sem_init(&sem_workerstart[i], 0, 0);
//        sem_init(&sem_workerend[i], 0, 0);
//    }
//
//    pthread_t handles[NUM_THREADS];
//    threadParam_t param[NUM_THREADS];
//
//    for (int k = 0; k < n; k++)
//    {
//        // 主线程进行除法操作
//        __m256 vt2 = _mm256_set1_ps(A[k][k]);
//        int j;
//        for (j = k + 1; j + 8 <= n; j += 8)
//        {
//            __m256 va2 = _mm256_loadu_ps(&(A[k][j]));
//            va2 = _mm256_div_ps(va2, vt2);
//            _mm256_store_ps(&(A[k][j]), va2);
//        }
//        for (; j < n; j++)
//            A[k][j] = A[k][j] * 1.0 / A[k][k];
//        A[k][k] = 1.0;
//
//        // 分配任务
//        for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//        {
//            param[t_id].k = k;
//            param[t_id].t_id = t_id;
//        }
//
//        // 创建线程
//        for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//            pthread_create(&handles[t_id], NULL, threadFunc, (void*)&param[t_id]);
//
//        // 主线程挂起等待所有的工作线程完成此轮消去工作
//        for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//            sem_wait(&sem_main);
//
//        // 唤醒所有工作线程进入下一轮
//        for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//            sem_post(&sem_workerend[t_id]);
//
//        // 等待所有工作线程结束
//        for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//            pthread_join(handles[t_id], NULL);
//    }
//
//    // 销毁信号量
//    sem_destroy(&sem_main);
//    for (int i = 0; i < NUM_THREADS; ++i)
//    {
//        sem_destroy(&sem_workerstart[i]);
//        sem_destroy(&sem_workerend[i]);
//    }
//
//	QueryPerformanceCounter((LARGE_INTEGER*)&tail);//结束计时
//	seconds = (tail - head) * 1000.0 / freq;//单位 ms
//
//cout << "pthread_jing: " << seconds << " ms" << endl;
//    //printMatrix();
//  return 0;
//}

///////////////////////////////////////////静态线程 + 信号量同步版本 + 三重循环全部纳入线程函数
//#include <iostream>
//#include <pthread.h>
//#include <semaphore.h>
//#include <windows.h>
//
//#include <pmmintrin.h>
//#include <xmmintrin.h>
//#include <emmintrin.h>
//#include <smmintrin.h>
//#include <tmmintrin.h>
//#include <nmmintrin.h>
//#include <immintrin.h> //AVX、AVX2
//
//
//using namespace std;
//
//const int n = 5000;
//float A[n][n];
//const int NUM_THREADS = 7; // 工作线程数量
//
//void init()
//{
//    for (int i = 0; i < n; i++)
//    {
//        for (int j = 0; j < n; j++)
//        {
//            A[i][j] = 0;
//        }
//        A[i][i] = 1.0;
//        for (int j = i + 1; j < n; j++)
//            A[i][j] = rand() % 100;
//    }
//
//    for (int i = 0; i < n; i++)
//    {
//        int k1 = rand() % n;
//        int k2 = rand() % n;
//        for (int j = 0; j < n; j++)
//        {
//            A[i][j] += A[0][j];
//            A[k1][j] += A[k2][j];
//        }
//    }
//}
//
//void f_ordinary()
//{
//    for (int k = 0; k < n; k++)
//    {
//        for (int j = k + 1; j < n; j++)
//        {
//            A[k][j] = A[k][j] * 1.0 / A[k][k];
//        }
//        A[k][k] = 1.0;
//
//        for (int i = k + 1; i < n; i++)
//        {
//            for (int j = k + 1; j < n; j++)
//            {
//                A[i][j] = A[i][j] - A[i][k] * A[k][j];
//            }
//            A[i][k] = 0;
//        }
//    }
//}
//
//struct threadParam_t {
//    int t_id; // 线程 id
//};
//
//// 信号量定义
//sem_t sem_leader;
//sem_t sem_Division[NUM_THREADS - 1];
//sem_t sem_Elimination[NUM_THREADS - 1];
//
//void* threadFunc(void* param) {
//    threadParam_t* p = (threadParam_t*)param;
//    int t_id = p->t_id;
//
//    for (int k = 0; k < n; ++k) {
//        if (t_id == 0) {
//            // 主线程进行除法操作
//            for (int j = k + 1; j < n; j++) {
//                A[k][j] = A[k][j] / A[k][k];
//            }
//            A[k][k] = 1.0;
//            // 唤醒其他线程进行消去操作
//            for (int i = 0; i < NUM_THREADS - 1; ++i) {
//                sem_post(&sem_Division[i]);
//            }
//        }
//        else {
//            // 工作线程等待除法完成
//            sem_wait(&sem_Division[t_id - 1]);
//        }
//
//        // 消去操作
//        for (int i = k + 1 + t_id; i < n; i += NUM_THREADS) {
//            for (int j = k + 1; j < n; ++j) {
//                A[i][j] -= A[i][k] * A[k][j];
//            }
//            A[i][k] = 0.0;
//        }
//
//        // 通知主线程已完成消去操作
//        if (t_id == 0) {
//            for (int i = 0; i < NUM_THREADS - 1; ++i) {
//                sem_wait(&sem_leader);
//            }
//            // 通知其他线程进入下一轮
//            for (int i = 0; i < NUM_THREADS - 1; ++i) {
//                sem_post(&sem_Elimination[i]);
//            }
//        }
//        else {
//            sem_post(&sem_leader); // 通知 leader, 已完成消去任务
//            sem_wait(&sem_Elimination[t_id - 1]); // 等待通知，进入下一轮
//        }
//    }
//
//    pthread_exit(NULL);
//    return 0;
//}
//void printMatrix() {
//    cout << "Matrix after Gaussian Elimination:" << endl;
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < n; j++) {
//            cout << A[i][j] << " ";
//        }
//        cout << endl;
//    }
//}
//int main() {
//    init();
//    long long counter;// 记录次数
//    double seconds;
//    long long head, tail, freq, noww;
//    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//    QueryPerformanceCounter((LARGE_INTEGER*)&head);//开始计时
//
//    // 初始化信号量
//    sem_init(&sem_leader, 0, 0);
//    for (int i = 0; i < NUM_THREADS - 1; ++i) {
//        sem_init(&sem_Division[i], 0, 0);
//        sem_init(&sem_Elimination[i], 0, 0);
//    }
//
//    pthread_t handles[NUM_THREADS];
//    threadParam_t params[NUM_THREADS];
//
//    // 创建线程
//    for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
//        params[t_id].t_id = t_id;
//        pthread_create(&handles[t_id], NULL, threadFunc, (void*)&params[t_id]);
//    }
//
//    // 等待线程结束
//    for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
//        pthread_join(handles[t_id], NULL);
//    }
//
//    // 销毁所有信号量
//    sem_destroy(&sem_leader);
//    for (int i = 0; i < NUM_THREADS - 1; ++i) {
//        sem_destroy(&sem_Division[i]);
//        sem_destroy(&sem_Elimination[i]);
//    }
//
//    QueryPerformanceCounter((LARGE_INTEGER*)&tail);//结束计时
//    seconds = (tail - head) * 1000.0 / freq;//单位 ms
//
//    cout << "pthread_jing: " << seconds << " ms" << endl;
//    //printMatrix();
//    return 0;
//}

//#include <iostream>
//#include <pthread.h>
//#include <semaphore.h>
//#include <windows.h>
//
//#include <pmmintrin.h>
//#include <xmmintrin.h>
//#include <emmintrin.h>
//#include <smmintrin.h>
//#include <tmmintrin.h>
//#include <nmmintrin.h>
//#include <immintrin.h> //AVX、AVX2
//
//using namespace std;
//
//const int max_n = 5500; // 最大矩阵规模
//const int step = 256;   // 矩阵规模递增步长
//int n; // 当前矩阵大小
//float A[max_n][max_n];
//const int NUM_THREADS = 7; // 工作线程数量
//
//void init(int n)
//{
//    for (int i = 0; i < n; i++)
//    {
//        for (int j = 0; j < n; j++)
//        {
//            A[i][j] = 0;
//        }
//        A[i][i] = 1.0;
//        for (int j = i + 1; j < n; j++)
//            A[i][j] = rand() % 100;
//    }
//
//    for (int i = 0; i < n; i++)
//    {
//        int k1 = rand() % n;
//        int k2 = rand() % n;
//        for (int j = 0; j < n; j++)
//        {
//            A[i][j] += A[0][j];
//            A[k1][j] += A[k2][j];
//        }
//    }
//}
//
//struct threadParam_t {
//    int t_id; // 线程 id
//};
//
//sem_t sem_leader;
//sem_t sem_Division[NUM_THREADS - 1];
//sem_t sem_Elimination[NUM_THREADS - 1];
//
//void* threadFunc(void* param) {
//    threadParam_t* p = (threadParam_t*)param;
//    int t_id = p->t_id;
//
//    for (int k = 0; k < n; ++k) {
//        if (t_id == 0) {
//            for (int j = k + 1; j < n; j++) {
//                A[k][j] = A[k][j] / A[k][k];
//            }
//            A[k][k] = 1.0;
//            for (int i = 0; i < NUM_THREADS - 1; ++i) {
//                sem_post(&sem_Division[i]);
//            }
//        }
//        else {
//            sem_wait(&sem_Division[t_id - 1]);
//        }
//
//        for (int i = k + 1 + t_id; i < n; i += NUM_THREADS) {
//            for (int j = k + 1; j < n; ++j) {
//                A[i][j] -= A[i][k] * A[k][j];
//            }
//            A[i][k] = 0.0;
//        }
//
//        if (t_id == 0) {
//            for (int i = 0; i < NUM_THREADS - 1; ++i) {
//                sem_wait(&sem_leader);
//            }
//            for (int i = 0; i < NUM_THREADS - 1; ++i) {
//                sem_post(&sem_Elimination[i]);
//            }
//        }
//        else {
//            sem_post(&sem_leader);
//            sem_wait(&sem_Elimination[t_id - 1]);
//        }
//    }
//
//    pthread_exit(NULL);
//    return 0;
//}
//
//int main() {
//    for (n = step; n <= max_n; n += step) {
//        init(n);
//        long long head, tail, freq;
//        double seconds;
//        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//        QueryPerformanceCounter((LARGE_INTEGER*)&head);//开始计时
//
//        // 初始化信号量
//        sem_init(&sem_leader, 0, 0);
//        for (int i = 0; i < NUM_THREADS - 1; ++i) {
//            sem_init(&sem_Division[i], 0, 0);
//            sem_init(&sem_Elimination[i], 0, 0);
//        }
//
//        pthread_t handles[NUM_THREADS];
//        threadParam_t params[NUM_THREADS];
//
//        // 创建线程
//        for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
//            params[t_id].t_id = t_id;
//            pthread_create(&handles[t_id], NULL, threadFunc, (void*)&params[t_id]);
//        }
//
//        // 等待线程结束
//        for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
//            pthread_join(handles[t_id], NULL);
//        }
//
//        // 销毁所有信号量
//        sem_destroy(&sem_leader);
//        for (int i = 0; i < NUM_THREADS - 1; ++i) {
//            sem_destroy(&sem_Division[i]);
//            sem_destroy(&sem_Elimination[i]);
//        }
//
//        QueryPerformanceCounter((LARGE_INTEGER*)&tail);//结束计时
//        seconds = (tail - head) * 1000.0 / freq;//单位 ms
//
//        cout << "Matrix size: " << n << "x" << n << endl;
//        cout << "pthread_jing: " << seconds << " ms" << endl << endl;
//    }
//
//    return 0;
//}

///////////////////////////////////////////////////////////////静态线程 +barrier 同步
//# include <iostream>
//# include <pthread.h>
//#include <semaphore.h>
//#include <windows.h>
//#include <pmmintrin.h>
//#include <xmmintrin.h>
//#include <emmintrin.h>
//#include <smmintrin.h>
//#include <tmmintrin.h>
//#include <nmmintrin.h>
//#include <immintrin.h> //AVX、AVX2
//using namespace std;
//
//const int n = 2000;
//float A[n][n];
//int NUM_THREADS = 7;
//void init()
//{
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			A[i][j] = 0;
//		}
//		A[i][i] = 1.0;
//		for (int j = i + 1; j < n; j++)
//			A[i][j] = rand() % 100;
//	}
//
//	for (int i = 0; i < n; i++)
//	{
//		int k1 = rand() % n;
//		int k2 = rand() % n;
//		for (int j = 0; j < n; j++)
//		{
//			A[i][j] += A[0][j];
//			A[k1][j] += A[k2][j];
//		}
//	}
//}
//
//void f_ordinary()
//{
//	for (int k = 0; k < n; k++)
//	{
//		for (int j = k + 1; j < n; j++)
//		{
//			A[k][j] = A[k][j] * 1.0 / A[k][k];
//		}
//		A[k][k] = 1.0;
//
//		for (int i = k + 1; i < n; i++)
//		{
//			for (int j = k + 1; j < n; j++)
//			{
//				A[i][j] = A[i][j] - A[i][k] * A[k][j];
//			}
//			A[i][k] = 0;
//		}
//	}
//}
//
//
//struct threadParam_t
//{
//	int t_id; //线程 id
//};
//
////barrier 定义
//pthread_barrier_t barrier_Divsion;
//pthread_barrier_t barrier_Elimination;
//
//
////线程函数定义
//void* threadFunc(void* param)
//{
//	__m256 va2, vt2, vx2, vaij2, vaik2, vakj2;
//
//
//
//	threadParam_t* p = (threadParam_t*)param;
//	int t_id = p->t_id;
//
//	for (int k = 0; k < n; ++k)
//	{
//		//vt2 = _mm256_set_ps(A[k][k], A[k][k], A[k][k], A[k][k], A[k][k], A[k][k], A[k][k], A[k][k]);
//		vt2 = _mm256_set1_ps(A[k][k]);
//		if (t_id == 0)
//		{
//			int j;
//			for (j = k + 1; j + 8 <= n; j += 8)
//			{
//				va2 = _mm256_loadu_ps(&(A[k][j]));
//				va2 = _mm256_div_ps(va2, vt2);
//				_mm256_store_ps(&(A[k][j]), va2);
//			}
//
//			for (; j < n; j++)
//			{
//				A[k][j] = A[k][j] * 1.0 / A[k][k];
//			}
//			A[k][k] = 1.0;
//		}
//
//		//第一个同步点
//		pthread_barrier_wait(&barrier_Divsion);
//
//		//循环划分任务（同学们可以尝试多种任务划分方式）
//		for (int i = k + 1 + t_id; i < n; i += NUM_THREADS)
//		{
//			//消去
//			//vaik2 = _mm256_set_ps(A[i][k], A[i][k], A[i][k], A[i][k], A[i][k], A[i][k], A[i][k], A[i][k]);
//			vaik2 = _mm256_set1_ps(A[i][k]);
//			int j;
//			for (j = k + 1; j + 8 <= n; j += 8)
//			{
//				vakj2 = _mm256_loadu_ps(&(A[k][j]));
//				vaij2 = _mm256_loadu_ps(&(A[i][j]));
//				vx2 = _mm256_mul_ps(vakj2, vaik2);
//				vaij2 = _mm256_sub_ps(vaij2, vx2);
//
//				_mm256_store_ps(&A[i][j], vaij2);
//			}
//			for (; j < n; j++)
//				A[i][j] = A[i][j] - A[i][k] * A[k][j];
//
//			A[i][k] = 0;
//		}
//		// 第二个同步点
//		pthread_barrier_wait(&barrier_Elimination);
//
//	}
//	pthread_exit(NULL);
//	return 0;
//}
//
//
//
//int main()
//{
//	init();
//	long long counter;// 记录次数
//	double seconds;
//	long long head, tail, freq, noww;
//	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//	QueryPerformanceCounter((LARGE_INTEGER*)&head);//开始计时
//
//
//	//初始化barrier
//	pthread_barrier_init(&barrier_Divsion, NULL, NUM_THREADS);
//	pthread_barrier_init(&barrier_Elimination, NULL, NUM_THREADS);
//
//
//	//创建线程
//	pthread_t* handles = new pthread_t[NUM_THREADS];// 创建对应的 Handle
//	threadParam_t* param = new threadParam_t[NUM_THREADS];// 创建对应的线程数据结构
//	for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//	{
//		param[t_id].t_id = t_id;
//		pthread_create(&handles[t_id], NULL, threadFunc, (void*)&param[t_id]);
//	}
//
//
//	for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//		pthread_join(handles[t_id], NULL);
//
//	//销毁所有的 barrier
//	pthread_barrier_destroy(&barrier_Divsion);
//	pthread_barrier_destroy(&barrier_Elimination);
//
//
//
//	QueryPerformanceCounter((LARGE_INTEGER*)&tail);//结束计时
//	seconds = (tail - head) * 1000.0 / freq;//单位 ms
//
//	cout << "pthread_neon_4: " << seconds << " ms" << endl;
//}

//




//#include <iostream>
//#include <pthread.h>
//#include <semaphore.h>
//#include <windows.h>
//#include <pmmintrin.h>
//#include <xmmintrin.h>
//#include <emmintrin.h>
//#include <smmintrin.h>
//#include <tmmintrin.h>
//#include <nmmintrin.h>
//#include <immintrin.h> //AVX、AVX2
//using namespace std;
//
//const int max_n = 5500; // 最大矩阵规模
//const int step = 256;   // 矩阵规模递增步长
//int n; // 当前矩阵规模，动态指定
//float A[max_n][max_n]; // 根据最大规模预分配空间
//int NUM_THREADS = 14;
//
//void init(int n)
//{
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			A[i][j] = 0;
//		}
//		A[i][i] = 1.0;
//		for (int j = i + 1; j < n; j++)
//			A[i][j] = rand() % 100;
//	}
//
//	for (int i = 0; i < n; i++)
//	{
//		int k1 = rand() % n;
//		int k2 = rand() % n;
//		for (int j = 0; j < n; j++)
//		{
//			A[i][j] += A[0][j];
//			A[k1][j] += A[k2][j];
//		}
//	}
//}
//
//
//struct threadParam_t
//{
//	int t_id; //线程 id
//};
//
////barrier 定义
//pthread_barrier_t barrier_Divsion;
//pthread_barrier_t barrier_Elimination;
//
//
////线程函数定义
//void* threadFunc(void* param)
//{
//	__m256 va2, vt2, vx2, vaij2, vaik2, vakj2;
//
//
//
//	threadParam_t* p = (threadParam_t*)param;
//	int t_id = p->t_id;
//
//	for (int k = 0; k < n; ++k)
//	{
//		vt2 = _mm256_set1_ps(A[k][k]);
//		if (t_id == 0)
//		{
//			int j;
//			for (j = k + 1; j + 8 <= n; j += 8)
//			{
//				va2 = _mm256_loadu_ps(&(A[k][j]));
//				va2 = _mm256_div_ps(va2, vt2);
//				_mm256_store_ps(&(A[k][j]), va2);
//			}
//
//			for (; j < n; j++)
//			{
//				A[k][j] = A[k][j] * 1.0 / A[k][k];
//			}
//			A[k][k] = 1.0;
//		}
//
//		//第一个同步点
//		pthread_barrier_wait(&barrier_Divsion);
//
//		//循环划分任务
//		for (int i = k + 1 + t_id; i < n; i += NUM_THREADS)
//		{
//			//消去
//			vaik2 = _mm256_set1_ps(A[i][k]);
//			int j;
//			for (j = k + 1; j + 8 <= n; j += 8)
//			{
//				vakj2 = _mm256_loadu_ps(&(A[k][j]));
//				vaij2 = _mm256_loadu_ps(&(A[i][j]));
//				vx2 = _mm256_mul_ps(vakj2, vaik2);
//				vaij2 = _mm256_sub_ps(vaij2, vx2);
//
//				_mm256_store_ps(&A[i][j], vaij2);
//			}
//			for (; j < n; j++)
//				A[i][j] = A[i][j] - A[i][k] * A[k][j];
//
//			A[i][k] = 0;
//		}
//		// 第二个同步点
//		pthread_barrier_wait(&barrier_Elimination);
//
//	}
//	pthread_exit(NULL);
//	return 0;
//}
//
//
//int main()
//{
//	for (n = step; n <= max_n; n += step) // 循环测试不同矩阵大小
//	{
//		init(n); // 根据当前矩阵大小初始化
//		long long counter;// 记录次数
//		double seconds;
//		long long head, tail, freq, noww;
//		QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//		QueryPerformanceCounter((LARGE_INTEGER*)&head);//开始计时
//
//		// 初始化barrier
//		pthread_barrier_init(&barrier_Divsion, NULL, NUM_THREADS);
//		pthread_barrier_init(&barrier_Elimination, NULL, NUM_THREADS);
//
//		// 创建线程
//		pthread_t* handles = new pthread_t[NUM_THREADS]; // 创建对应的 Handle
//		threadParam_t* param = new threadParam_t[NUM_THREADS]; // 创建对应的线程数据结构
//		for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//		{
//			param[t_id].t_id = t_id;
//			pthread_create(&handles[t_id], NULL, threadFunc, (void*)&param[t_id]);
//		}
//
//		for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//			pthread_join(handles[t_id], NULL);
//
//		// 销毁所有的 barrier
//		pthread_barrier_destroy(&barrier_Divsion);
//		pthread_barrier_destroy(&barrier_Elimination);
//
//		QueryPerformanceCounter((LARGE_INTEGER*)&tail); // 结束计时
//		seconds = (tail - head) * 1000.0 / freq; // 单位 ms
//
//		cout << "Matrix size: " << n << "x" << n << endl;
//		cout << "pthread_neon_4: " << seconds << " ms" << endl << endl;
//	}
//}



























//#include <iostream>
//#include <pthread.h>
//#include <semaphore.h>
//#include <windows.h>
//#include <pmmintrin.h>
//#include <xmmintrin.h>
//#include <emmintrin.h>
//#include <smmintrin.h>
//#include <tmmintrin.h>
//#include <nmmintrin.h>
//#include <immintrin.h> //AVX、AVX2
//using namespace std;
//
//const int max_n = 5500; // 最大矩阵规模
//const int step = 256;   // 矩阵规模递增步长
//int n; // 当前矩阵规模，动态指定
//float A[max_n][max_n]; // 根据最大规模预分配空间
//const int NUM_THREADS = 7;
//
//void init(int n)
//{
//	for (int i = 0; i < n; i++)
//	{
//		for (int j = 0; j < n; j++)
//		{
//			A[i][j] = 0;
//		}
//		A[i][i] = 1.0;
//		for (int j = i + 1; j < n; j++)
//			A[i][j] = rand() % 100;
//	}
//
//	for (int i = 0; i < n; i++)
//	{
//		int k1 = rand() % n;
//		int k2 = rand() % n;
//		for (int j = 0; j < n; j++)
//		{
//			A[i][j] += A[0][j];
//			A[k1][j] += A[k2][j];
//		}
//	}
//}
//
//struct threadParam_t
//{
//    int k;     // 消去的轮次
//    int t_id;  // 线程 id
//};
//
//sem_t sem_main;
//sem_t sem_workerstart[NUM_THREADS];
//sem_t sem_workerend[NUM_THREADS];
//
//void* threadFunc(void* param)
//{
//    __m256 va, vt, vx, vaij, vaik, vakj;
//
//    threadParam_t* p = (threadParam_t*)param;
//    int k = p->k;
//    int t_id = p->t_id;
//    int i = k + t_id + 1;
//
//    for (int m = k + 1 + t_id; m < n; m += NUM_THREADS)
//    {
//        vaik = _mm256_set1_ps(A[k][k]);
//        int j;
//        for (j = k + 1; j + 8 <= n; j += 8)
//        {
//            vakj = _mm256_loadu_ps(&(A[k][j]));
//            vaij = _mm256_loadu_ps(&(A[m][j]));
//            vx = _mm256_mul_ps(vakj, vaik);
//            vaij = _mm256_sub_ps(vaij, vx);
//            _mm256_store_ps(&A[i][j], vaij);
//        }
//        for (; j < n; j++)
//            A[m][j] = A[m][j] - A[m][k] * A[k][j];
//
//        A[m][k] = 0;
//    }
//
//    sem_post(&sem_main); // 唤醒主线程
//    sem_wait(&sem_workerend[t_id]); // 等待主线程唤醒进入下一轮
//
//    pthread_exit(NULL);
//    return 0;
//}
//
//void printMatrix() {
//    cout << "Matrix after Gaussian Elimination:" << endl;
//    for (int i = 0; i < n; i++) {
//        for (int j = 0; j < n; j++) {
//            cout << A[i][j] << " ";
//        }
//        cout << endl;
//    }
//}
//
//
//int main()
//{
//	for (n = step; n <= max_n; n += step) // 循环测试不同矩阵大小
//	{
//		init(n); // 根据当前矩阵大小初始化
//		long long counter;// 记录次数
//		double seconds;
//		long long head, tail, freq, noww;
//		QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//		QueryPerformanceCounter((LARGE_INTEGER*)&head);//开始计时
//
//        sem_init(&sem_main, 0, 0);
//        for (int i = 0; i < NUM_THREADS; ++i)
//        {
//            sem_init(&sem_workerstart[i], 0, 0);
//            sem_init(&sem_workerend[i], 0, 0);
//        }
//
//        pthread_t handles[NUM_THREADS];
//        threadParam_t param[NUM_THREADS];
//
//        for (int k = 0; k < n; k++)
//        {
//            // 主线程进行除法操作
//            __m256 vt2 = _mm256_set1_ps(A[k][k]);
//            int j;
//            for (j = k + 1; j + 8 <= n; j += 8)
//            {
//                __m256 va2 = _mm256_loadu_ps(&(A[k][j]));
//                va2 = _mm256_div_ps(va2, vt2);
//                _mm256_store_ps(&(A[k][j]), va2);
//            }
//            for (; j < n; j++)
//                A[k][j] = A[k][j] * 1.0 / A[k][k];
//            A[k][k] = 1.0;
//
//            // 分配任务
//            for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//            {
//                param[t_id].k = k;
//                param[t_id].t_id = t_id;
//            }
//
//            // 创建线程
//            for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//                pthread_create(&handles[t_id], NULL, threadFunc, (void*)&param[t_id]);
//
//            // 主线程挂起等待所有的工作线程完成此轮消去工作
//            for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//                sem_wait(&sem_main);
//
//            // 唤醒所有工作线程进入下一轮
//            for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//                sem_post(&sem_workerend[t_id]);
//
//            // 等待所有工作线程结束
//            for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//                pthread_join(handles[t_id], NULL);
//        }
//
//        // 销毁信号量
//        sem_destroy(&sem_main);
//        for (int i = 0; i < NUM_THREADS; ++i)
//        {
//            sem_destroy(&sem_workerstart[i]);
//            sem_destroy(&sem_workerend[i]);
//        }
//
//		QueryPerformanceCounter((LARGE_INTEGER*)&tail); // 结束计时
//		seconds = (tail - head) * 1000.0 / freq; // 单位 ms
//
//		cout << "Matrix size: " << n << "x" << n << endl;
//		cout << "pthread_neon_4: " << seconds << " ms" << endl << endl;
//	}
//}