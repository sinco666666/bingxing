//#include<iostream>
//#include<fstream>
//#include<sstream>
//#include <ctime>
//#include <ratio>
//#include <chrono>
//#include<pthread.h>
//#include<semaphore.h>
//using namespace std;
//const int col = 8399, elinenum = 4535, num_thread = 8; //列数、被消元行数
//int tmp = 0;
//int bytenum = (col - 1) / 32 + 1;   //每个实例中的byte型数组数
//typedef struct {
//    int t_id;   //线程id
//}threadparam_t;
//pthread_barrier_t barrier;
//pthread_barrier_t barrier2;
//pthread_barrier_t barrier3;
//class bitmatrix {
//public:
//    int mycol;    //首项
//    int* mybyte;
//    bitmatrix() {    //初始化
//        mycol = -1;
//        mybyte = new int[bytenum];
//        for (int i = 0; i < bytenum; i++)
//            mybyte[i] = 0;
//    }
//    bool isnull() {  //判断当前行是否为空行
//        if (mycol == -1)return 1;
//        return 0;
//    }
//    void insert(int x) { //数据读入
//        if (mycol == -1)mycol = x;
//        int a = x / 32, b = x % 32;
//        mybyte[a] |= (1 << b);
//    }
//};
//bitmatrix* eliminer = new bitmatrix[col], * eline = new bitmatrix[elinenum];
//void readdata() {
//    ifstream ifs;
//    ifs.open("x8399.txt");  //消元子
//    string temp;
//    while (getline(ifs, temp)) {
//        istringstream ss(temp);
//        int x;
//        int trow = 0;
//        while (ss >> x) {
//            if (!trow)trow = x;    //第一个读入元素代表行号
//            eliminer[trow].insert(x);
//        }
//    }
//    ifs.close();
//    ifstream ifs2;
//    ifs2.open("b8399.txt");     //被消元行,读入方式与消元子不同
//    int trow = 0;
//    while (getline(ifs2, temp)) {
//        istringstream ss(temp);
//        int x;
//        while (ss >> x) {
//            eline[trow].insert(x);
//        }
//        trow++;
//    }
//    ifs2.close();
//}
//void* threadfunc(void* param) {
//    threadparam_t* p = (threadparam_t*)param;
//    int t_id = p->t_id;
//    for (int i = col - 1; i >= 0; i--)
//        if (!eliminer[i].isnull()) {
//            for (int j = t_id; j < elinenum; j += num_thread) {
//                if (eline[j].mycol == i) {
//                    for (int k = 0; k < bytenum; k++)
//                        eline[j].mybyte[k] ^= eliminer[i].mybyte[k];
//                    bool f = 1;
//                    for (int p = bytenum - 1; p >= 0 && f; p--)
//                        for (int k = 31; k >= 0 && f; k--)
//                            if ((eline[j].mybyte[p] & (1 << k)) != 0) {
//                                eline[j].mycol = p * 32 + k;
//                                f = 0;
//                            }
//                    if (f)eline[j].mycol = -1;
//                }
//            }
//        }
//        else {
//            pthread_barrier_wait(&barrier);
//            if (t_id == 0)
//                for (int j = 0; j < elinenum; j++) {
//                    if (eline[j].mycol == i) {
//                        eliminer[i] = eline[j];
//                        tmp = j + 1;
//                        break;
//                    }
//                    tmp = j + 2;
//                }
//            pthread_barrier_wait(&barrier2);
//            int temp = t_id;
//            while (temp < tmp)temp += num_thread;
//            for (int j = temp; j < elinenum; j += num_thread) {
//                if (eline[j].mycol == i) {
//                    for (int k = 0; k < bytenum; k++)
//                        eline[j].mybyte[k] ^= eliminer[i].mybyte[k];
//                    bool f = 1;
//                    for (int p = bytenum - 1; p >= 0 && f; p--)
//                        for (int k = 31; k >= 0 && f; k--)
//                            if ((eline[j].mybyte[p] & (1 << k)) != 0) {
//                                eline[j].mycol = p * 32 + k;
//                                f = 0;
//                            }
//                    if (f)eline[j].mycol = -1;
//                }
//            }
//        }
//    pthread_exit(NULL);
//    return 0;
//}
//void dowork() {  //消元
//    pthread_barrier_init(&barrier, NULL, num_thread);
//    pthread_barrier_init(&barrier2, NULL, num_thread);
//    pthread_barrier_init(&barrier3, NULL, num_thread);
//    //创建线程
//    pthread_t handles[num_thread];
//    threadparam_t param[num_thread];
//    for (int t_id = 0; t_id < num_thread; t_id++) {
//        param[t_id].t_id = t_id;
//        pthread_create(&handles[t_id], NULL, threadfunc, (void*)&param[t_id]);
//    }
//    for (int i = 0; i < num_thread; i++)
//        pthread_join(handles[i], NULL);
//    pthread_barrier_destroy(&barrier);
//    pthread_barrier_destroy(&barrier2);
//    pthread_barrier_destroy(&barrier3);
//}
//void printres() { //打印结果
//    for (int i = 0; i < elinenum; i++) {
//        if (eline[i].isnull()) { puts(""); continue; }   //空行的特殊情况
//        for (int j = bytenum - 1; j >= 0; j--) {
//            for (int k = 31; k >= 0; k--)
//                if ((eline[i].mybyte[j] & (1 << k)) != 0) {     //一个错误调了半小时，谨记当首位为1时>>不等于除法！
//                    printf("%d ", j * 32 + k);
//                }
//        }
//        puts("");
//    }
//}
//int main() {
//    readdata();
//    using namespace std::chrono;
//    high_resolution_clock::time_point t1 = high_resolution_clock::now();
//    dowork();
//    high_resolution_clock::time_point t2 = high_resolution_clock::now();
//    std::cout << "serial: " << duration_cast<duration<double>>(t2 - t1).count() << std::endl;
//    //printres();
//    system("pause");
//    return 0;
//}
//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <ctime>
//#include <ratio>
//#include <chrono>
//#include <pthread.h>
//#include <semaphore.h>
//
//using namespace std;
//
//const int col = 37960, elinenum = 14921, num_thread = 8; //列数、被消元行数
//int tmp = 0;
//int bytenum = (col - 1) / 32 + 1;   //每个实例中的byte型数组数
//
//typedef struct {
//    int t_id;   //线程id
//} threadparam_t;
//
//pthread_barrier_t barrier;
//pthread_barrier_t barrier2;
//pthread_barrier_t barrier3;
//sem_t semaphore;
//
//class bitmatrix {
//public:
//    int mycol;    //首项
//    int* mybyte;
//
//    bitmatrix() {    //初始化
//        mycol = -1;
//        mybyte = new int[bytenum];
//        for (int i = 0; i < bytenum; i++)
//            mybyte[i] = 0;
//    }
//
//    bool isnull() {  //判断当前行是否为空行
//        if (mycol == -1) return 1;
//        return 0;
//    }
//
//    void insert(int x) { //数据读入
//        if (mycol == -1) mycol = x;
//        int a = x / 32, b = x % 32;
//        mybyte[a] |= (1 << b);
//    }
//};
//
//bitmatrix* eliminer = new bitmatrix[col], * eline = new bitmatrix[elinenum];
//
//void readdata() {
//    ifstream ifs;
//    ifs.open("x37960.txt");  //消元子
//    string temp;
//    while (getline(ifs, temp)) {
//        istringstream ss(temp);
//        int x;
//        int trow = 0;
//        while (ss >> x) {
//            if (!trow) trow = x;    //第一个读入元素代表行号
//            eliminer[trow].insert(x);
//        }
//    }
//    ifs.close();
//
//    ifstream ifs2;
//    ifs2.open("b37960.txt");     //被消元行,读入方式与消元子不同
//    int trow = 0;
//    while (getline(ifs2, temp)) {
//        istringstream ss(temp);
//        int x;
//        while (ss >> x) {
//            eline[trow].insert(x);
//        }
//        trow++;
//    }
//    ifs2.close();
//}
//
//void* threadfunc_barrier(void* param) {
//    threadparam_t* p = (threadparam_t*)param;
//    int t_id = p->t_id;
//
//    for (int i = col - 1; i >= 0; i--) {
//        if (!eliminer[i].isnull()) {
//            for (int j = t_id; j < elinenum; j += num_thread) {
//                if (eline[j].mycol == i) {
//                    for (int k = 0; k < bytenum; k++)
//                        eline[j].mybyte[k] ^= eliminer[i].mybyte[k];
//                    bool f = 1;
//                    for (int p = bytenum - 1; p >= 0 && f; p--)
//                        for (int k = 31; k >= 0 && f; k--)
//                            if ((eline[j].mybyte[p] & (1 << k)) != 0) {
//                                eline[j].mycol = p * 32 + k;
//                                f = 0;
//                            }
//                    if (f) eline[j].mycol = -1;
//                }
//            }
//        }
//        else {
//            pthread_barrier_wait(&barrier);
//            if (t_id == 0) {
//                for (int j = 0; j < elinenum; j++) {
//                    if (eline[j].mycol == i) {
//                        eliminer[i] = eline[j];
//                        tmp = j + 1;
//                        break;
//                    }
//                    tmp = j + 2;
//                }
//            }
//            pthread_barrier_wait(&barrier2);
//            int temp = t_id;
//            while (temp < tmp) temp += num_thread;
//            for (int j = temp; j < elinenum; j += num_thread) {
//                if (eline[j].mycol == i) {
//                    for (int k = 0; k < bytenum; k++)
//                        eline[j].mybyte[k] ^= eliminer[i].mybyte[k];
//                    bool f = 1;
//                    for (int p = bytenum - 1; p >= 0 && f; p--)
//                        for (int k = 31; k >= 0 && f; k--)
//                            if ((eline[j].mybyte[p] & (1 << k)) != 0) {
//                                eline[j].mycol = p * 32 + k;
//                                f = 0;
//                            }
//                    if (f) eline[j].mycol = -1;
//                }
//            }
//        }
//    }
//    pthread_exit(NULL);
//    return 0;
//}
//
//void* threadfunc_semaphore(void* param) {
//    threadparam_t* p = (threadparam_t*)param;
//    int t_id = p->t_id;
//
//    for (int i = col - 1; i >= 0; i--) {
//        if (!eliminer[i].isnull()) {
//            for (int j = t_id; j < elinenum; j += num_thread) {
//                if (eline[j].mycol == i) {
//                    for (int k = 0; k < bytenum; k++)
//                        eline[j].mybyte[k] ^= eliminer[i].mybyte[k];
//                    bool f = 1;
//                    for (int p = bytenum - 1; p >= 0 && f; p--)
//                        for (int k = 31; k >= 0 && f; k--)
//                            if ((eline[j].mybyte[p] & (1 << k)) != 0) {
//                                eline[j].mycol = p * 32 + k;
//                                f = 0;
//                            }
//                    if (f) eline[j].mycol = -1;
//                }
//            }
//        }
//        else {
//            sem_wait(&semaphore);
//            if (t_id == 0) {
//                for (int j = 0; j < elinenum; j++) {
//                    if (eline[j].mycol == i) {
//                        eliminer[i] = eline[j];
//                        tmp = j + 1;
//                        break;
//                    }
//                    tmp = j + 2;
//                }
//            }
//            sem_post(&semaphore);
//            int temp = t_id;
//            while (temp < tmp) temp += num_thread;
//            for (int j = temp; j < elinenum; j += num_thread) {
//                if (eline[j].mycol == i) {
//                    for (int k = 0; k < bytenum; k++)
//                        eline[j].mybyte[k] ^= eliminer[i].mybyte[k];
//                    bool f = 1;
//                    for (int p = bytenum - 1; p >= 0 && f; p--)
//                        for (int k = 31; k >= 0 && f; k--)
//                            if ((eline[j].mybyte[p] & (1 << k)) != 0) {
//                                eline[j].mycol = p * 32 + k;
//                                f = 0;
//                            }
//                    if (f) eline[j].mycol = -1;
//                }
//            }
//        }
//    }
//    pthread_exit(NULL);
//    return 0;
//}
//
//void* threadfunc_dynamic(void* param) {
//    threadparam_t* p = (threadparam_t*)param;
//    int t_id = p->t_id;
//
//    for (int i = col - 1; i >= 0; i--) {
//        if (!eliminer[i].isnull()) {
//            for (int j = t_id; j < elinenum; j += num_thread) {
//                if (eline[j].mycol == i) {
//                    for (int k = 0; k < bytenum; k++)
//                        eline[j].mybyte[k] ^= eliminer[i].mybyte[k];
//                    bool f = 1;
//                    for (int p = bytenum - 1; p >= 0 && f; p--)
//                        for (int k = 31; k >= 0 && f; k--)
//                            if ((eline[j].mybyte[p] & (1 << k)) != 0) {
//                                eline[j].mycol = p * 32 + k;
//                                f = 0;
//                            }
//                    if (f) eline[j].mycol = -1;
//                }
//            }
//        }
//        else {
//            if (t_id == 0) {
//                for (int j = 0; j < elinenum; j++) {
//                    if (eline[j].mycol == i) {
//                        eliminer[i] = eline[j];
//                        tmp = j + 1;
//                        break;
//
//                    }
//                    tmp = j + 2;
//                }
//            }
//            int temp = t_id;
//            while (temp < tmp) temp += num_thread;
//            for (int j = temp; j < elinenum; j += num_thread) {
//                if (eline[j].mycol == i) {
//                    for (int k = 0; k < bytenum; k++)
//                        eline[j].mybyte[k] ^= eliminer[i].mybyte[k];
//                    bool f = 1;
//                    for (int p = bytenum - 1; p >= 0 && f; p--)
//                        for (int k = 31; k >= 0 && f; k--)
//                            if ((eline[j].mybyte[p] & (1 << k)) != 0) {
//                                eline[j].mycol = p * 32 + k;
//                                f = 0;
//                            }
//                    if (f) eline[j].mycol = -1;
//                }
//            }
//        }
//    }
//    pthread_exit(NULL);
//    return 0;
//}
//
//void dowork_barrier() {
//    pthread_barrier_init(&barrier, NULL, num_thread);
//    pthread_barrier_init(&barrier2, NULL, num_thread);
//    pthread_barrier_init(&barrier3, NULL, num_thread);
//
//    pthread_t handles[num_thread];
//    threadparam_t param[num_thread];
//
//    for (int t_id = 0; t_id < num_thread; t_id++) {
//        param[t_id].t_id = t_id;
//        pthread_create(&handles[t_id], NULL, threadfunc_barrier, (void*)&param[t_id]);
//    }
//    for (int i = 0; i < num_thread; i++)
//        pthread_join(handles[i], NULL);
//
//    pthread_barrier_destroy(&barrier);
//    pthread_barrier_destroy(&barrier2);
//    pthread_barrier_destroy(&barrier3);
//}
//
//void dowork_semaphore() {
//    sem_init(&semaphore, 0, 1);
//
//    pthread_t handles[num_thread];
//    threadparam_t param[num_thread];
//
//    for (int t_id = 0; t_id < num_thread; t_id++) {
//        param[t_id].t_id = t_id;
//        pthread_create(&handles[t_id], NULL, threadfunc_semaphore, (void*)&param[t_id]);
//    }
//    for (int i = 0; i < num_thread; i++)
//        pthread_join(handles[i], NULL);
//
//    sem_destroy(&semaphore);
//}
//
//void dowork_dynamic() {
//    pthread_t handles[num_thread];
//    threadparam_t param[num_thread];
//
//    for (int t_id = 0; t_id < num_thread; t_id++) {
//        param[t_id].t_id = t_id;
//        pthread_create(&handles[t_id], NULL, threadfunc_dynamic, (void*)&param[t_id]);
//    }
//    for (int i = 0; i < num_thread; i++)
//        pthread_join(handles[i], NULL);
//}
//
//void printres() { //打印结果
//    for (int i = 0; i < elinenum; i++) {
//        if (eline[i].isnull()) { puts(""); continue; }   //空行的特殊情况
//        for (int j = bytenum - 1; j >= 0; j--) {
//            for (int k = 31; k >= 0; k--)
//                if ((eline[i].mybyte[j] & (1 << k)) != 0) {
//                    printf("%d ", j * 32 + k);
//                }
//        }
//        puts("");
//    }
//}
//
//int main() {
//    readdata();
//    using namespace std::chrono;
//
//    
//    high_resolution_clock::time_point t1 = high_resolution_clock::now();
//    dowork_barrier();
//    high_resolution_clock::time_point t2 = high_resolution_clock::now();
//    std::cout << "Barrier synchronization: " << duration_cast<duration<double>>(t2 - t1).count() << " seconds" << std::endl;
//    
//
//    readdata(); // 重新加载数据
//    t1 = high_resolution_clock::now();
//    dowork_semaphore();
//    t2 = high_resolution_clock::now();
//    std::cout << "Semaphore synchronization: " << duration_cast<duration<double>>(t2 - t1).count() << " seconds" << std::endl;
//    
//
//    
//    readdata(); // 重新加载数据
//    t1 = high_resolution_clock::now();
//    dowork_dynamic();
//    t2 = high_resolution_clock::now();
//    std::cout << "Dynamic threading: " << duration_cast<duration<double>>(t2 - t1).count() << " seconds" << std::endl;
//    
//    //printres();
//    return 0;
//}
// 
























//openmp
//#include<iostream>
//#include<fstream>
//#include<sstream>
//#include <ctime>
//#include <ratio>
//#include <chrono>
//#include<omp.h>
//#include <emmintrin.h>
//#include <immintrin.h>
//using namespace std;
//const int col = 130, elinenum = 8, num_thread = 8;
//bool parallel = 1; //列数、被消元行数
//int tmp = 0;
//int bytenum = (col - 1) / 32 + 1;   //每个实例中的byte型数组数
//class bitmatrix {
//public:
//    int mycol;    //首项
//    int* mybyte;
//    bitmatrix() {    //初始化
//        mycol = -1;
//        mybyte = (int*)_aligned_malloc(bytenum * sizeof(int), 32);
//        for (int i = 0; i < bytenum; i++)
//            mybyte[i] = 0;
//    }
//    ~bitmatrix() {
//        _aligned_free(mybyte);
//    }
//    bool isnull() {  //判断当前行是否为空行
//        if (mycol == -1) return true;
//        return false;
//    }
//    void insert(int x) { //数据读入
//        if (mycol == -1) mycol = x;
//        int a = x / 32, b = x % 32;
//        mybyte[a] |= (1 << b);
//    }
//    void doxor(bitmatrix& b) {  //两行做异或操作，由于结果留在本实例中，只有被消元行能执行这一操作,且异或操作后要更新首项
//        for (int i = 0; i < bytenum; i++)
//            mybyte[i] ^= b.mybyte[i];
//        update_mycol();
//    }
//    void doxor2(bitmatrix& b) {  // SSE版本异或
//        int i = 0;
//        for (; i + 4 <= bytenum; i += 4) {
//            __m128i byte1 = _mm_load_si128((__m128i*) & mybyte[i]);
//            __m128i byte2 = _mm_load_si128((__m128i*) & b.mybyte[i]);
//            byte1 = _mm_xor_si128(byte1, byte2);
//            _mm_store_si128((__m128i*) & mybyte[i], byte1);
//        }
//        for (; i < bytenum; i++)
//            mybyte[i] ^= b.mybyte[i];
//        update_mycol();
//    }
//    void doxor3(bitmatrix& b) {  // AVX版本异或
//        int i = 0;
//        for (; i + 8 <= bytenum; i += 8) {
//            __m256i byte1 = _mm256_load_si256((__m256i*) & mybyte[i]);
//            __m256i byte2 = _mm256_load_si256((__m256i*) & b.mybyte[i]);
//            byte1 = _mm256_xor_si256(byte1, byte2);
//            _mm256_store_si256((__m256i*) & mybyte[i], byte1);
//        }
//        for (; i < bytenum; i++)
//            mybyte[i] ^= b.mybyte[i];
//        update_mycol();
//    }
//private:
//    void update_mycol() {
//        for (int i = bytenum - 1; i >= 0; i--)
//            for (int j = 31; j >= 0; j--)
//                if ((mybyte[i] & (1 << j)) != 0) {
//                    mycol = i * 32 + j;
//                    return;
//                }
//        mycol = -1;
//    }
//};
//
//
//bitmatrix* eliminer = new bitmatrix[col], * eline = new bitmatrix[elinenum];
//void readdata() {
//    ifstream ifs;
//    ifs.open("x130.txt");  //消元子
//    string temp;
//    while (getline(ifs, temp)) {
//        istringstream ss(temp);
//        int x;
//        int trow = 0;
//        while (ss >> x) {
//            if (!trow)trow = x;    //第一个读入元素代表行号
//            eliminer[trow].insert(x);
//        }
//    }
//    ifs.close();
//    ifstream ifs2;
//    ifs2.open("b130.txt");     //被消元行,读入方式与消元子不同
//    int trow = 0;
//    while (getline(ifs2, temp)) {
//        istringstream ss(temp);
//        int x;
//        while (ss >> x) {
//            eline[trow].insert(x);
//        }
//        trow++;
//    }
//    ifs2.close();
//}
//
//void dowork() {  //串行消元--被消元行->消元子
//    for (int i = 0; i < elinenum; i++) {
//        while (!eline[i].isnull()) {  //只要被消元行非空，循环处理
//            int tcol = eline[i].mycol;  //被消元行的首项
//            if (!eliminer[tcol].isnull())    //如果存在对应消元子
//                eline[i].doxor(eliminer[tcol]);
//            else {
//                eliminer[tcol] = eline[i];    //由于被消元行升格为消元子后不参与后续处理，可以直接用=来浅拷贝
//                break;
//            }
//        }
//    }
//}
//
//
//void dowork1() {     //openmp消元
//    int i, j;
//#pragma omp parallel if(parallel),num_threads(num_thread),private(i,j)
//    for (i = col - 1; i >= 0; i--)
//        if (!eliminer[i].isnull()) {
//#pragma omp for 
//            for (j = 0; j < elinenum; j++) {
//                if (eline[j].mycol == i)
//                    eline[j].doxor(eliminer[i]);
//            }
//        }
//        else {
//#pragma omp barrier
//#pragma omp single
//            for (j = 0; j < elinenum; j++) {
//                if (eline[j].mycol == i) {
//                    eliminer[i] = eline[j];
//                    tmp = j + 1;
//                    break;
//                }
//                tmp = j + 2;
//            }
//#pragma omp for
//            for (j = tmp; j < elinenum; j++) {
//                if (eline[j].mycol == i)
//                    eline[j].doxor(eliminer[i]);
//            }
//        }
//}
//
//void dowork2() {     // openmp + SSE消元
//    int i, j;
//#pragma omp parallel if(parallel), num_threads(num_thread), private(i, j)
//    for (i = col - 1; i >= 0; i--) {
//        if (!eliminer[i].isnull()) {
//#pragma omp for 
//            for (j = 0; j < elinenum; j++) {
//                if (eline[j].mycol == i)
//                    eline[j].doxor2(eliminer[i]);   // SSE异或
//            }
//        }
//        else {
//#pragma omp barrier
//#pragma omp single
//            for (j = 0; j < elinenum; j++) {
//                if (eline[j].mycol == i) {
//                    eliminer[i] = eline[j];
//                    tmp = j + 1;
//                    break;
//                }
//                tmp = j + 2;
//            }
//#pragma omp for
//            for (j = tmp; j < elinenum; j++) {
//                if (eline[j].mycol == i)
//                    eline[j].doxor2(eliminer[i]);
//            }
//        }
//    }
//}
//
//void dowork3() {     // openmp + AVX消元
//    int i, j;
//#pragma omp parallel if(parallel), num_threads(num_thread), private(i, j)
//    for (i = col - 1; i >= 0; i--) {
//        if (!eliminer[i].isnull()) {
//#pragma omp for 
//            for (j = 0; j < elinenum; j++) {
//                if (eline[j].mycol == i)
//                    eline[j].doxor3(eliminer[i]);   // AVX异或
//            }
//        }
//        else {
//#pragma omp barrier
//#pragma omp single
//            for (j = 0; j < elinenum; j++) {
//                if (eline[j].mycol == i) {
//                    eliminer[i] = eline[j];
//                    tmp = j + 1;
//                    break;
//                }
//                tmp = j + 2;
//            }
//#pragma omp for
//            for (j = tmp; j < elinenum; j++) {
//                if (eline[j].mycol == i)
//                    eline[j].doxor3(eliminer[i]);
//            }
//        }
//    }
//}
//
//
//
//void printres() { //打印结果
//    for (int i = 0; i < elinenum; i++) {
//        if (eline[i].isnull()) { puts(""); continue; }   //空行的特殊情况
//        for (int j = bytenum - 1; j >= 0; j--) {
//            for (int k = 31; k >= 0; k--)
//                if ((eline[i].mybyte[j] & (1 << k)) != 0) {     //一个错误调了半小时，谨记当首位为1时>>不等于除法！
//                    printf("%d ", j * 32 + k);
//                }
//        }
//        puts("");
//    }
//}
//int main() {
//    readdata();
//    using namespace std::chrono;
//    /*high_resolution_clock::time_point t1 = high_resolution_clock::now();
//    dowork1();
//    high_resolution_clock::time_point t2 = high_resolution_clock::now();
//    std::cout << "openmp: " << duration_cast<duration<double>>(t2 - t1).count() << std::endl;*/
//    
//    readdata();
//    high_resolution_clock::time_point t3 = high_resolution_clock::now();
//    dowork2();
//    high_resolution_clock::time_point t4 = high_resolution_clock::now();
//    std::cout << "openmp+sse: " << duration_cast<duration<double>>(t4 - t3).count() << std::endl;
//
//
//    readdata();
//    high_resolution_clock::time_point t5 = high_resolution_clock::now();
//    dowork3();
//    high_resolution_clock::time_point t6 = high_resolution_clock::now();
//    std::cout << "openmp+AVX: " << duration_cast<duration<double>>(t6 - t5).count() << std::endl;
//    //printres();
//    //system("pause");
//    return 0;
//}