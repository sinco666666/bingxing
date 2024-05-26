# include <iostream>
# include <arm_neon.h> // use Neon
# include <pthread.h>
# include <semaphore.h>
# include <sys/time.h>
using namespace std;

int n; // ����Ĵ�С
float** A; // �޸�Ϊָ��ָ���ָ��
int NUM_THREADS = 14;

void init(int size)
{
    n = size; // �������size��ֵ��ȫ�ֱ���n
    A = new float* [n];
    for (int i = 0; i < n; i++)
    {
        A[i] = new float[n];
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



int main()
{
    const int max_n = 5500; // �������С
    const int step = 256; // �����С���Ӳ���

    for (int size = step; size <= max_n; size += step)
    {
        init(size); // ��ʼ����ǰ��С�ľ���

        struct timeval head, tail;
        double seconds;
        gettimeofday(&head, NULL); // ��ʼ��ʱ

        f_ordinary();

        gettimeofday(&tail, NULL); // ������ʱ
        seconds = ((tail.tv_sec - head.tv_sec) * 1000000 + (tail.tv_usec - head.tv_usec)) / 1000.0; // ��λ ms
        cout << "Matrix size: " << n << "x" << n << endl;
        cout << "pthread_neon_4: " << seconds << " ms" << endl << endl;

        // ����̬������ڴ�
        for (int i = 0; i < n; i++)
        {
            delete[] A[i];
        }
        delete[] A;
    }
}