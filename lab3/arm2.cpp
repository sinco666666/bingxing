#include <iostream>
#include <arm_neon.h> // use Neon
#include <pthread.h>
#include <semaphore.h>
#include <sys/time.h>
using namespace std;

const int max_n = 5500; // 最大矩阵规模
const int step = 256; // 矩阵规模递增步长
int NUM_THREADS = 7;
float** A;

void init(int n)
{
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

struct threadParam_t
{
    int t_id; //线程 id
    int n; // 矩阵大小
};

//信号量定义
sem_t sem_leader;
sem_t* sem_Divsion; // 每个线程有自己专属的信号量
sem_t* sem_Elimination;

//垂直穿插划分
void* threadFunc_vertica2(void* param)
{

    threadParam_t* p = (threadParam_t*)param;
    int t_id = p->t_id;
    int n = p->n;

    for (int k = 0; k < n; ++k)
    {

        if (t_id == 0)
        {
            for (int j = k + 1; j < n; j++)
            {
                A[k][j] = A[k][j] * 1.0 / A[k][k];
            }
            A[k][k] = 1.0;
        }
        else
        {
            sem_wait(&sem_Divsion[t_id - 1]); // 阻塞，等待完成除法操作
        }

        // t_id 为 0 的线程唤醒其它工作线程，进行消去操作
        if (t_id == 0)
        {
            for (int i = 0; i < NUM_THREADS - 1; ++i)
                sem_post(&sem_Divsion[i]);
        }


        //循环划分任务
        for (int i = k + 1; i < n; i++)
        {
            for (int j = k + 1 + t_id; j < n; j += NUM_THREADS)
            {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
            A[i][k] = 0;
        }


        if (t_id == 0)
        {
            for (int i = 0; i < NUM_THREADS - 1; ++i)
                sem_wait(&sem_leader); // 等待其它 worker 完成消去

            for (int i = 0; i < NUM_THREADS - 1; ++i)
                sem_post(&sem_Elimination[i]); // 通知其它 worker 进入下一轮
        }
        else
        {
            sem_post(&sem_leader);// 通知 leader, 已完成消去任务
            sem_wait(&sem_Elimination[t_id - 1]); // 等待通知，进入下一轮
        }
    }
    pthread_exit(NULL);
    return NULL;
}
int main()
{
    // 初始化信号量和其他准备工作，比如：
    sem_Divsion = new sem_t[NUM_THREADS - 1];
    sem_Elimination = new sem_t[NUM_THREADS - 1];

    for (int n = step; n <= max_n; n += step)
    {
        init(n);
        struct timeval head, tail;
        double seconds;
        gettimeofday(&head, NULL); //开始计时

        sem_init(&sem_leader, 0, 0);
        for (int i = 0; i < NUM_THREADS - 1; ++i)
        {
            sem_init(&sem_Divsion[i], 0, 0);
            sem_init(&sem_Elimination[i], 0, 0);
        }

        //创建线程
        pthread_t* handles = new pthread_t[NUM_THREADS];// 创建对应的 Handle
        threadParam_t* param = new threadParam_t[NUM_THREADS];// 创建对应的线程数据结构
        for (int t_id = 0; t_id < NUM_THREADS; t_id++)
        {
            param[t_id].t_id = t_id;
            param[t_id].n = n;
            pthread_create(&handles[t_id], NULL, threadFunc_vertica2, (void*)&param[t_id]);
        }

        for (int t_id = 0; t_id < NUM_THREADS; t_id++)
            pthread_join(handles[t_id], NULL);

        gettimeofday(&tail, NULL);//结束计时
        seconds = ((tail.tv_sec - head.tv_sec) * 1000000 + (tail.tv_usec - head.tv_usec)) / 1000.0;//单位 ms
        cout << "Matrix size: " << n << "x" << n << endl;
        cout << "pthread_neon_3: " << seconds << " ms" << endl;

        // 清理动态分配的内存
        for (int i = 0; i < n; i++) {
            delete[] A[i];
        }
        delete[] A;

        //销毁所有信号量
        sem_destroy(&sem_leader);
        for (int i = 0; i < NUM_THREADS - 1; ++i)
        {
            sem_destroy(&sem_Divsion[i]);
            sem_destroy(&sem_Elimination[i]);
        }
    }

    delete[] sem_Divsion;
    delete[] sem_Elimination;
}