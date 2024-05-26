
//////////////////////////////////////////////////////��Ч�ĳ�����
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
//#include <immintrin.h> ��
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
//	int k; //��ȥ���ִ�
//	int t_id; // �߳� id
//};
//
//void* threadFunc(void* param) {
//	threadParam_t* p = (threadParam_t*)param;
//	int k = p->k; //��ȥ���ִ�
//	int t_id = p->t_id; //�̱߳��
//	int i = k + t_id + 1; //��ȡ�Լ��ļ�������
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
//	long long counter; // ��¼����
//	double seconds;
//	long long head, tail, freq, noww;
//	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//	QueryPerformanceCounter((LARGE_INTEGER*)&head); //��ʼ��ʱ
//
//	for (int k = 0; k < n; ++k) {
//		//���߳�����������
//		for (int j = k + 1; j < n; j++) {
//			A[k][j] = A[k][j] / A[k][k];
//		}
//		A[k][k] = 1.0;
//
//		//���������̣߳�������ȥ����
//		int worker_count = n - 1 - k; //�����߳�����
//		pthread_t* handles = new pthread_t[worker_count]; // ������Ӧ�� Handle
//		threadParam_t* params = new threadParam_t[worker_count]; // ������Ӧ���߳����ݽṹ
//
//		//��������
//		for (int t_id = 0; t_id < worker_count; t_id++) {
//			params[t_id].k = k;
//			params[t_id].t_id = t_id;
//			pthread_create(&handles[t_id], NULL, threadFunc, (void*)&params[t_id]);
//		}
//
//		//���̹߳���ȴ����еĹ����߳���ɴ�����ȥ����
//		for (int t_id = 0; t_id < worker_count; t_id++) {
//			pthread_join(handles[t_id], NULL);
//		}
//
//		delete[] handles;
//		delete[] params;
//	}
//
//	QueryPerformanceCounter((LARGE_INTEGER*)&tail); //������ʱ
//	seconds = (tail - head) * 1000.0 / freq; //��λ ms
//
//	cout << "pthread_dong: " << seconds << " ms" << endl;
//}

//////////////////////////////////////////////////////��̬�̰߳汾
# include <iostream>
# include <pthread.h>
#include <windows.h>

#include <pmmintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <tmmintrin.h>
#include <nmmintrin.h>
#include <immintrin.h> //AVX��AVX2
using namespace std;

const int n = 2000;
float A[n][n];
int worker_count = 7; //�����߳�����
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
	int k; //��ȥ���ִ�
	int t_id; // �߳� id
};

void *threadFunc(void* param)
{

	__m256 va, vt, vx, vaij, vaik, vakj;

	threadParam_t* p = (threadParam_t*)param;
	int k = p->k; //��ȥ���ִ�
	int t_id = p->t_id; //�̱߳��
	int i = k + t_id + 1; //��ȡ�Լ��ļ�������
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

	long long counter;// ��¼����
	double seconds;
	long long head, tail, freq, noww;
	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
	QueryPerformanceCounter((LARGE_INTEGER*)&head);//��ʼ��ʱ


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

		//���������̣߳�������ȥ����

		pthread_t* handles = new pthread_t[worker_count];// ������Ӧ�� Handle
		threadParam_t* param = new threadParam_t[worker_count];// ������Ӧ���߳����ݽṹ

		//��������
		for (int t_id = 0; t_id < worker_count; t_id++)
		{
			param[t_id].k = k;
			param[t_id].t_id = t_id;
		}
		//�����߳�
		for (int t_id = 0; t_id < worker_count; t_id++)
			pthread_create(&handles[t_id], NULL, threadFunc, (void*)&param[t_id]);

		//���̹߳���ȴ����еĹ����߳���ɴ�����ȥ����
		for (int t_id = 0; t_id < worker_count; t_id++)
			pthread_join(handles[t_id], NULL);

	}


	QueryPerformanceCounter((LARGE_INTEGER*)&tail);//������ʱ
	seconds = (tail - head) * 1000.0 / freq;//��λ ms

	cout << "pthread_dong: " << seconds << " ms" << endl;

}


/////////////////////////////////////////////��̬���Դ���
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
//#include <immintrin.h> //AVX��AVX2
//
//using namespace std;
//
//const int max_n = 6000; // �������ģ
//const int step = 256;   // �����ģ��������
//const int worker_count = 3; // �����߳�����
//int n;
//
//float A[max_n][max_n]; // ȫ������ A�������ģ���Ϊ max_n
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
//    int k;     // ��ȥ���ִ�
//    int t_id;  // �߳� id
//};
//
//void* threadFunc(void* param)
//{
//    __m256 va, vt, vx, vaij, vaik, vakj;
//
//    threadParam_t* p = (threadParam_t*)param;
//    int k = p->k;       // ��ȥ���ִ�
//    int t_id = p->t_id; // �̱߳��
//    int i = k + t_id + 1; // ��ȡ�Լ��ļ�������
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
//        QueryPerformanceCounter((LARGE_INTEGER*)&head); // ��ʼ��ʱ
//
//        // ���в��и�˹��Ԫ
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
//            // ���������̣߳�������ȥ����
//            pthread_t* handles = new pthread_t[worker_count]; // ������Ӧ�� Handle
//            threadParam_t* param = new threadParam_t[worker_count]; // ������Ӧ���߳����ݽṹ
//
//            // ��������
//            for (int t_id = 0; t_id < worker_count; t_id++)
//            {
//                param[t_id].k = k;
//                param[t_id].t_id = t_id;
//            }
//            // �����߳�
//            for (int t_id = 0; t_id < worker_count; t_id++)
//                pthread_create(&handles[t_id], NULL, threadFunc, (void*)&param[t_id]);
//
//            // ���̹߳���ȴ����еĹ����߳���ɴ�����ȥ����
//            for (int t_id = 0; t_id < worker_count; t_id++)
//                pthread_join(handles[t_id], NULL);
//        }
//
//        QueryPerformanceCounter((LARGE_INTEGER*)&tail); // ������ʱ
//        parallel_time = (tail - head) * 1000.0 / freq; // ��λ ms
//
//        cout << "Matrix size: " << n << "x" << n << endl;
//        cout << "pthread_dong: " << parallel_time << " ms" << endl;
//        cout << endl;
//
//        // ���Դ��м���ʱ��
//        init(n);
//        long long head1, tail1, freq1;
//        double seconds1;
//        QueryPerformanceFrequency((LARGE_INTEGER*)&freq1);
//        QueryPerformanceCounter((LARGE_INTEGER*)&head1); // ��ʼ��ʱ
//
//        f_ordinary(n);
//
//        QueryPerformanceCounter((LARGE_INTEGER*)&tail1); // ������ʱ
//        seconds1 = (tail1 - head1) * 1000.0 / freq1; // ��λ ms
//
//        cout << "f_ordinary: " << seconds1 << " ms" << endl;
//
//        // ������ٱ�
//        double speedup = seconds1 / parallel_time;
//        cout << "Speedup: " << speedup << endl;
//        cout << endl;
//    }
//
//    return 0;
//}





//////////////////////////////////////////��̬
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
//#include <immintrin.h> //AVX��AVX2
//using namespace std;
//int n;
//const int max_n = 6000; // �������ģ
//const int step = 256;   // �����ģ��������
//int worker_count = 7; // �����߳�����
//float A[max_n][max_n]; // ȫ������ A�������ģ���Ϊ max_n
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
//// ���ຯ����f_ordinary, threadParam_t, threadFunc�����ֲ���
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
//    int k; //��ȥ���ִ�
//    int t_id; // �߳� id
//};
//
//void* threadFunc(void* param)
//{
//
//    __m256 va, vt, vx, vaij, vaik, vakj;
//
//    threadParam_t* p = (threadParam_t*)param;
//    int k = p->k; //��ȥ���ִ�
//    int t_id = p->t_id; //�̱߳��
//    int i = k + t_id + 1; //��ȡ�Լ��ļ�������
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
//        init(n); // ʹ�õ�ǰ�����С���г�ʼ��
//        __m256 va2, vt2, vx2, vaij2, vaik2, vakj2;
//
//        long long counter; // ��¼����
//        double seconds;
//        long long head, tail, freq, noww;
//        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//        QueryPerformanceCounter((LARGE_INTEGER*)&head); // ��ʼ��ʱ
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
//            //���������̣߳�������ȥ����
//
//            pthread_t* handles = new pthread_t[worker_count];// ������Ӧ�� Handle
//            threadParam_t* param = new threadParam_t[worker_count];// ������Ӧ���߳����ݽṹ
//
//            //��������
//            for (int t_id = 0; t_id < worker_count; t_id++)
//            {
//                param[t_id].k = k;
//                param[t_id].t_id = t_id;
//            }
//            //�����߳�
//            for (int t_id = 0; t_id < worker_count; t_id++)
//                pthread_create(&handles[t_id], NULL, threadFunc, (void*)&param[t_id]);
//
//            //���̹߳���ȴ����еĹ����߳���ɴ�����ȥ����
//            for (int t_id = 0; t_id < worker_count; t_id++)
//                pthread_join(handles[t_id], NULL);
//
//        }
//
//        QueryPerformanceCounter((LARGE_INTEGER*)&tail); // ������ʱ
//        seconds = (tail - head) * 1000.0 / freq; // ��λ ms
//
//        cout << "Matrix size: " << n << "x" << n << endl;
//        cout << "pthread_dong: " << seconds << " ms" << endl;
//        cout << endl;
//    }
//}



///////////////////////////////////////////��̬�߳� + �ź���ͬ���汾
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
//#include <immintrin.h> //AVX��AVX2
//using namespace std;
//
//const int n = 2000;
//float A[n][n];
//const int NUM_THREADS = 5; // �����߳�����
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
//    int k;     // ��ȥ���ִ�
//    int t_id;  // �߳� id
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
//    sem_post(&sem_main); // �������߳�
//    sem_wait(&sem_workerend[t_id]); // �ȴ����̻߳��ѽ�����һ��
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
//    long long counter;// ��¼����
//    double seconds;
//    long long head, tail, freq, noww;
//    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//    QueryPerformanceCounter((LARGE_INTEGER*)&head);//��ʼ��ʱ
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
//        // ���߳̽��г�������
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
//        // ��������
//        for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//        {
//            param[t_id].k = k;
//            param[t_id].t_id = t_id;
//        }
//
//        // �����߳�
//        for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//            pthread_create(&handles[t_id], NULL, threadFunc, (void*)&param[t_id]);
//
//        // ���̹߳���ȴ����еĹ����߳���ɴ�����ȥ����
//        for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//            sem_wait(&sem_main);
//
//        // �������й����߳̽�����һ��
//        for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//            sem_post(&sem_workerend[t_id]);
//
//        // �ȴ����й����߳̽���
//        for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//            pthread_join(handles[t_id], NULL);
//    }
//
//    // �����ź���
//    sem_destroy(&sem_main);
//    for (int i = 0; i < NUM_THREADS; ++i)
//    {
//        sem_destroy(&sem_workerstart[i]);
//        sem_destroy(&sem_workerend[i]);
//    }
//
//	QueryPerformanceCounter((LARGE_INTEGER*)&tail);//������ʱ
//	seconds = (tail - head) * 1000.0 / freq;//��λ ms
//
//cout << "pthread_jing: " << seconds << " ms" << endl;
//    //printMatrix();
//  return 0;
//}

///////////////////////////////////////////��̬�߳� + �ź���ͬ���汾 + ����ѭ��ȫ�������̺߳���
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
//#include <immintrin.h> //AVX��AVX2
//
//
//using namespace std;
//
//const int n = 5000;
//float A[n][n];
//const int NUM_THREADS = 7; // �����߳�����
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
//    int t_id; // �߳� id
//};
//
//// �ź�������
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
//            // ���߳̽��г�������
//            for (int j = k + 1; j < n; j++) {
//                A[k][j] = A[k][j] / A[k][k];
//            }
//            A[k][k] = 1.0;
//            // ���������߳̽�����ȥ����
//            for (int i = 0; i < NUM_THREADS - 1; ++i) {
//                sem_post(&sem_Division[i]);
//            }
//        }
//        else {
//            // �����̵߳ȴ��������
//            sem_wait(&sem_Division[t_id - 1]);
//        }
//
//        // ��ȥ����
//        for (int i = k + 1 + t_id; i < n; i += NUM_THREADS) {
//            for (int j = k + 1; j < n; ++j) {
//                A[i][j] -= A[i][k] * A[k][j];
//            }
//            A[i][k] = 0.0;
//        }
//
//        // ֪ͨ���߳��������ȥ����
//        if (t_id == 0) {
//            for (int i = 0; i < NUM_THREADS - 1; ++i) {
//                sem_wait(&sem_leader);
//            }
//            // ֪ͨ�����߳̽�����һ��
//            for (int i = 0; i < NUM_THREADS - 1; ++i) {
//                sem_post(&sem_Elimination[i]);
//            }
//        }
//        else {
//            sem_post(&sem_leader); // ֪ͨ leader, �������ȥ����
//            sem_wait(&sem_Elimination[t_id - 1]); // �ȴ�֪ͨ��������һ��
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
//    long long counter;// ��¼����
//    double seconds;
//    long long head, tail, freq, noww;
//    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//    QueryPerformanceCounter((LARGE_INTEGER*)&head);//��ʼ��ʱ
//
//    // ��ʼ���ź���
//    sem_init(&sem_leader, 0, 0);
//    for (int i = 0; i < NUM_THREADS - 1; ++i) {
//        sem_init(&sem_Division[i], 0, 0);
//        sem_init(&sem_Elimination[i], 0, 0);
//    }
//
//    pthread_t handles[NUM_THREADS];
//    threadParam_t params[NUM_THREADS];
//
//    // �����߳�
//    for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
//        params[t_id].t_id = t_id;
//        pthread_create(&handles[t_id], NULL, threadFunc, (void*)&params[t_id]);
//    }
//
//    // �ȴ��߳̽���
//    for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
//        pthread_join(handles[t_id], NULL);
//    }
//
//    // ���������ź���
//    sem_destroy(&sem_leader);
//    for (int i = 0; i < NUM_THREADS - 1; ++i) {
//        sem_destroy(&sem_Division[i]);
//        sem_destroy(&sem_Elimination[i]);
//    }
//
//    QueryPerformanceCounter((LARGE_INTEGER*)&tail);//������ʱ
//    seconds = (tail - head) * 1000.0 / freq;//��λ ms
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
//#include <immintrin.h> //AVX��AVX2
//
//using namespace std;
//
//const int max_n = 5500; // �������ģ
//const int step = 256;   // �����ģ��������
//int n; // ��ǰ�����С
//float A[max_n][max_n];
//const int NUM_THREADS = 7; // �����߳�����
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
//    int t_id; // �߳� id
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
//        QueryPerformanceCounter((LARGE_INTEGER*)&head);//��ʼ��ʱ
//
//        // ��ʼ���ź���
//        sem_init(&sem_leader, 0, 0);
//        for (int i = 0; i < NUM_THREADS - 1; ++i) {
//            sem_init(&sem_Division[i], 0, 0);
//            sem_init(&sem_Elimination[i], 0, 0);
//        }
//
//        pthread_t handles[NUM_THREADS];
//        threadParam_t params[NUM_THREADS];
//
//        // �����߳�
//        for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
//            params[t_id].t_id = t_id;
//            pthread_create(&handles[t_id], NULL, threadFunc, (void*)&params[t_id]);
//        }
//
//        // �ȴ��߳̽���
//        for (int t_id = 0; t_id < NUM_THREADS; t_id++) {
//            pthread_join(handles[t_id], NULL);
//        }
//
//        // ���������ź���
//        sem_destroy(&sem_leader);
//        for (int i = 0; i < NUM_THREADS - 1; ++i) {
//            sem_destroy(&sem_Division[i]);
//            sem_destroy(&sem_Elimination[i]);
//        }
//
//        QueryPerformanceCounter((LARGE_INTEGER*)&tail);//������ʱ
//        seconds = (tail - head) * 1000.0 / freq;//��λ ms
//
//        cout << "Matrix size: " << n << "x" << n << endl;
//        cout << "pthread_jing: " << seconds << " ms" << endl << endl;
//    }
//
//    return 0;
//}

///////////////////////////////////////////////////////////////��̬�߳� +barrier ͬ��
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
//#include <immintrin.h> //AVX��AVX2
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
//	int t_id; //�߳� id
//};
//
////barrier ����
//pthread_barrier_t barrier_Divsion;
//pthread_barrier_t barrier_Elimination;
//
//
////�̺߳�������
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
//		//��һ��ͬ����
//		pthread_barrier_wait(&barrier_Divsion);
//
//		//ѭ����������ͬѧ�ǿ��Գ��Զ������񻮷ַ�ʽ��
//		for (int i = k + 1 + t_id; i < n; i += NUM_THREADS)
//		{
//			//��ȥ
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
//		// �ڶ���ͬ����
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
//	long long counter;// ��¼����
//	double seconds;
//	long long head, tail, freq, noww;
//	QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//	QueryPerformanceCounter((LARGE_INTEGER*)&head);//��ʼ��ʱ
//
//
//	//��ʼ��barrier
//	pthread_barrier_init(&barrier_Divsion, NULL, NUM_THREADS);
//	pthread_barrier_init(&barrier_Elimination, NULL, NUM_THREADS);
//
//
//	//�����߳�
//	pthread_t* handles = new pthread_t[NUM_THREADS];// ������Ӧ�� Handle
//	threadParam_t* param = new threadParam_t[NUM_THREADS];// ������Ӧ���߳����ݽṹ
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
//	//�������е� barrier
//	pthread_barrier_destroy(&barrier_Divsion);
//	pthread_barrier_destroy(&barrier_Elimination);
//
//
//
//	QueryPerformanceCounter((LARGE_INTEGER*)&tail);//������ʱ
//	seconds = (tail - head) * 1000.0 / freq;//��λ ms
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
//#include <immintrin.h> //AVX��AVX2
//using namespace std;
//
//const int max_n = 5500; // �������ģ
//const int step = 256;   // �����ģ��������
//int n; // ��ǰ�����ģ����ָ̬��
//float A[max_n][max_n]; // ��������ģԤ����ռ�
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
//	int t_id; //�߳� id
//};
//
////barrier ����
//pthread_barrier_t barrier_Divsion;
//pthread_barrier_t barrier_Elimination;
//
//
////�̺߳�������
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
//		//��һ��ͬ����
//		pthread_barrier_wait(&barrier_Divsion);
//
//		//ѭ����������
//		for (int i = k + 1 + t_id; i < n; i += NUM_THREADS)
//		{
//			//��ȥ
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
//		// �ڶ���ͬ����
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
//	for (n = step; n <= max_n; n += step) // ѭ�����Բ�ͬ�����С
//	{
//		init(n); // ���ݵ�ǰ�����С��ʼ��
//		long long counter;// ��¼����
//		double seconds;
//		long long head, tail, freq, noww;
//		QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//		QueryPerformanceCounter((LARGE_INTEGER*)&head);//��ʼ��ʱ
//
//		// ��ʼ��barrier
//		pthread_barrier_init(&barrier_Divsion, NULL, NUM_THREADS);
//		pthread_barrier_init(&barrier_Elimination, NULL, NUM_THREADS);
//
//		// �����߳�
//		pthread_t* handles = new pthread_t[NUM_THREADS]; // ������Ӧ�� Handle
//		threadParam_t* param = new threadParam_t[NUM_THREADS]; // ������Ӧ���߳����ݽṹ
//		for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//		{
//			param[t_id].t_id = t_id;
//			pthread_create(&handles[t_id], NULL, threadFunc, (void*)&param[t_id]);
//		}
//
//		for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//			pthread_join(handles[t_id], NULL);
//
//		// �������е� barrier
//		pthread_barrier_destroy(&barrier_Divsion);
//		pthread_barrier_destroy(&barrier_Elimination);
//
//		QueryPerformanceCounter((LARGE_INTEGER*)&tail); // ������ʱ
//		seconds = (tail - head) * 1000.0 / freq; // ��λ ms
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
//#include <immintrin.h> //AVX��AVX2
//using namespace std;
//
//const int max_n = 5500; // �������ģ
//const int step = 256;   // �����ģ��������
//int n; // ��ǰ�����ģ����ָ̬��
//float A[max_n][max_n]; // ��������ģԤ����ռ�
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
//    int k;     // ��ȥ���ִ�
//    int t_id;  // �߳� id
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
//    sem_post(&sem_main); // �������߳�
//    sem_wait(&sem_workerend[t_id]); // �ȴ����̻߳��ѽ�����һ��
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
//	for (n = step; n <= max_n; n += step) // ѭ�����Բ�ͬ�����С
//	{
//		init(n); // ���ݵ�ǰ�����С��ʼ��
//		long long counter;// ��¼����
//		double seconds;
//		long long head, tail, freq, noww;
//		QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//		QueryPerformanceCounter((LARGE_INTEGER*)&head);//��ʼ��ʱ
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
//            // ���߳̽��г�������
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
//            // ��������
//            for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//            {
//                param[t_id].k = k;
//                param[t_id].t_id = t_id;
//            }
//
//            // �����߳�
//            for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//                pthread_create(&handles[t_id], NULL, threadFunc, (void*)&param[t_id]);
//
//            // ���̹߳���ȴ����еĹ����߳���ɴ�����ȥ����
//            for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//                sem_wait(&sem_main);
//
//            // �������й����߳̽�����һ��
//            for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//                sem_post(&sem_workerend[t_id]);
//
//            // �ȴ����й����߳̽���
//            for (int t_id = 0; t_id < NUM_THREADS; t_id++)
//                pthread_join(handles[t_id], NULL);
//        }
//
//        // �����ź���
//        sem_destroy(&sem_main);
//        for (int i = 0; i < NUM_THREADS; ++i)
//        {
//            sem_destroy(&sem_workerstart[i]);
//            sem_destroy(&sem_workerend[i]);
//        }
//
//		QueryPerformanceCounter((LARGE_INTEGER*)&tail); // ������ʱ
//		seconds = (tail - head) * 1000.0 / freq; // ��λ ms
//
//		cout << "Matrix size: " << n << "x" << n << endl;
//		cout << "pthread_neon_4: " << seconds << " ms" << endl << endl;
//	}
//}