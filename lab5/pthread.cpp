//#include<iostream>
//#include<fstream>
//#include<sstream>
//#include <ctime>
//#include <ratio>
//#include <chrono>
//#include<pthread.h>
//#include<semaphore.h>
//using namespace std;
//const int col = 254, elinenum = 53, num_thread = 8; //����������Ԫ����
//int tmp = 0;
//int bytenum = (col - 1) / 32 + 1;   //ÿ��ʵ���е�byte��������
//typedef struct {
//    int t_id;   //�߳�id
//}threadparam_t;
//pthread_barrier_t barrier;
//pthread_barrier_t barrier2;
//pthread_barrier_t barrier3;
//class bitmatrix {
//public:
//    int mycol;    //����
//    int* mybyte;
//    bitmatrix() {    //��ʼ��
//        mycol = -1;
//        mybyte = new int[bytenum];
//        for (int i = 0; i < bytenum; i++)
//            mybyte[i] = 0;
//    }
//    bool isnull() {  //�жϵ�ǰ���Ƿ�Ϊ����
//        if (mycol == -1)return 1;
//        return 0;
//    }
//    void insert(int x) { //���ݶ���
//        if (mycol == -1)mycol = x;
//        int a = x / 32, b = x % 32;
//        mybyte[a] |= (1 << b);
//    }
//};
//bitmatrix* eliminer = new bitmatrix[col], * eline = new bitmatrix[elinenum];
//void readdata() {
//    ifstream ifs;
//    ifs.open("x254.txt");  //��Ԫ��
//    string temp;
//    while (getline(ifs, temp)) {
//        istringstream ss(temp);
//        int x;
//        int trow = 0;
//        while (ss >> x) {
//            if (!trow)trow = x;    //��һ������Ԫ�ش����к�
//            eliminer[trow].insert(x);
//        }
//    }
//    ifs.close();
//    ifstream ifs2;
//    ifs2.open("b254.txt");     //����Ԫ��,���뷽ʽ����Ԫ�Ӳ�ͬ
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
//void dowork() {  //��Ԫ
//    pthread_barrier_init(&barrier, NULL, num_thread);
//    pthread_barrier_init(&barrier2, NULL, num_thread);
//    pthread_barrier_init(&barrier3, NULL, num_thread);
//    //�����߳�
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
//void printres() { //��ӡ���
//    for (int i = 0; i < elinenum; i++) {
//        if (eline[i].isnull()) { puts(""); continue; }   //���е��������
//        for (int j = bytenum - 1; j >= 0; j--) {
//            for (int k = 31; k >= 0; k--)
//                if ((eline[i].mybyte[j] & (1 << k)) != 0) {     //һ��������˰�Сʱ�����ǵ���λΪ1ʱ>>�����ڳ�����
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
//    printres();
//    system("pause");
//    return 0;
//}