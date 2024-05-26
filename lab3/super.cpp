//#include <pmmintrin.h>
//#include <xmmintrin.h>
//#include <emmintrin.h>
//#include <smmintrin.h>
//#include <tmmintrin.h>
//#include <nmmintrin.h>
//#include <iostream>
//#include <sstream>
//#include <fstream>
//#include <immintrin.h> 
//#include <windows.h>
//#include <cstring> 
//#include <vector>
//using namespace std;
//
//// 数据规模和数组将在 main 中动态定义
//unsigned int Act[37960][1188]; // 最大规模的数组预设大小
//unsigned int Pas[37960][1188];
//
//int Num;
//int pasNum;
//int lieNum;
//
//void init_A(const string& filename)
//{
//    ifstream infile(filename);
//    char fin[10000] = { 0 };
//    int index;
//
//    while (infile.getline(fin, sizeof(fin)))
//    {
//        stringstream line(fin);
//        int biaoji = 0;
//        unsigned int a;
//
//        while (line >> a)
//        {
//            if (biaoji == 0)
//            {
//                index = a;
//                biaoji = 1;
//            }
//            int k = a % 32;
//            int j = a / 32;
//
//            int temp = 1 << k;
//            Act[index][Num - 1 - j] |= temp; // Use bitwise OR to set the bit
//        }
//        Act[index][Num] = 1; // Mark this row as non-empty
//    }
//}
//
//void init_P(const string& filename)
//{
//    ifstream infile(filename);
//    char fin[10000] = { 0 };
//    int index = 0;
//
//    while (infile.getline(fin, sizeof(fin)))
//    {
//        stringstream line(fin);
//        int biaoji = 0;
//        unsigned int a;
//
//        while (line >> a)
//        {
//            if (biaoji == 0)
//            {
//                Pas[index][Num] = a;
//                biaoji = 1;
//            }
//            int k = a % 32;
//            int j = a / 32;
//
//            int temp = 1 << k;
//            Pas[index][Num - 1 - j] |= temp;
//        }
//        index++;
//    }
//}
//
//
//void f_serial() {
//    int i;
//    for (i = lieNum - 1; i >= 0; i -= 8) {
//        for (int j = 0; j < pasNum; j++) {
//            int idxPas = Pas[j][Num];
//            if (idxPas <= i && idxPas >= i - 7) {
//                int index = idxPas;
//                if (Act[index][Num] == 1) {
//                    for (int k = 0; k < Num; k++) {
//                        Pas[j][k] ^= Act[index][k];
//                    }
//                    int S_num = 0;
//                    for (int num = 0; num < Num; num++) {
//                        if (Pas[j][num] != 0) {
//                            for (int bit = 0; bit < 32; bit++) {
//                                if ((Pas[j][num] & (1 << bit)) != 0) {
//                                    S_num = bit + num * 32;
//                                    break;
//                                }
//                            }
//                            if (S_num != 0) break;
//                        }
//                    }
//                    Pas[j][Num] = S_num;
//                }
//                else {
//                    memcpy(Act[index], Pas[j], Num * sizeof(Pas[0][0]));
//                    Act[index][Num] = 1;
//                }
//            }
//        }
//    }
//}
//
//// 定义线程参数结构体
//struct ThreadData {
//    int startIdx; // 开始的索引
//    int endIdx;   // 结束的索引
//};
//
//DWORD WINAPI ThreadFunc(LPVOID lpParam) {
//    ThreadData* pData = (ThreadData*)lpParam;
//    int startIdx = pData->startIdx;  // 线程处理的起始索引
//    int endIdx = pData->endIdx;      // 线程处理的结束索引
//
//    // 以下仅为示例逻辑，实际逻辑应根据f_cache_optimized或其他函数的处理方式进行修改
//    for (int i = startIdx; i <= endIdx; i++) {
//        // 假设我们在这里处理Act和Pas数组
//        for (int j = 0; j < pasNum; j++) {
//            unsigned int idxPas = Pas[j][Num];
//            if (idxPas <= i && idxPas >= i - 7) {
//                unsigned int index = idxPas;
//                unsigned int* pasPtr = Pas[j];
//                unsigned int* actPtr = Act[index];
//                if (actPtr[Num] == 1) {
//                    for (int k = 0; k < Num; k++) {
//                        pasPtr[k] ^= actPtr[k];
//                    }
//
//                    unsigned int S_num = 0;
//                    for (int num = 0; num < Num; num++) {
//                        if (pasPtr[num] != 0) {
//                            for (int bit = 0; bit < 32; bit++) {
//                                if ((pasPtr[num] & (1 << bit)) != 0) {
//                                    S_num = bit + num * 32;
//                                    break;
//                                }
//                            }
//                            if (S_num != 0) break;
//                        }
//                    }
//                    pasPtr[Num] = S_num;
//                }
//                else {
//                    memcpy(actPtr, pasPtr, Num * sizeof(unsigned int)); // Use cached pointers
//                    actPtr[Num] = 1;
//                }
//            }
//        }
//    }
//
//    return 0;
//}
//
//void f_parallel() {
//    int threadCount = 7; // 或根据需要设置线程数量
//    std::vector<HANDLE> threads(threadCount);
//    std::vector<ThreadData> threadData(threadCount);
//
//    int partSize = lieNum / threadCount;
//    for (int i = 0; i < threadCount; ++i) {
//        threadData[i].startIdx = i * partSize;
//        threadData[i].endIdx = (i + 1) * partSize - 1;
//        if (i == threadCount - 1) threadData[i].endIdx = lieNum - 1; // 最后一个线程处理剩余部分
//
//        threads[i] = CreateThread(NULL, 0, ThreadFunc, &threadData[i], 0, NULL);
//    }
//
//    // 等待所有线程完成
//    WaitForMultipleObjects(threadCount, threads.data(), TRUE, INFINITE);
//}
//
//
//
//int main()
//{
//    // Data sizes, corresponding pasNum and Num values
//    const int dataSize[] = { 130, 254, 562, 1011, 2362, 3799, 8399, 23045, 37960 };
//    const int pasNums[] = { 8, 53, 53, 263, 453, 1953, 4535, 14325, 14921 };
//    const int Nums[] = { 5, 8, 18, 32, 74, 119, 263, 721, 1187 };
//    const int dataCount = sizeof(dataSize) / sizeof(dataSize[0]);
//
//    // Loop over each data size for testing
//    for (int i = 0; i < dataCount; ++i)
//    {
//        lieNum = dataSize[i];
//        pasNum = pasNums[i];
//        Num = Nums[i];
//
//        string xFilename = "x" + to_string(dataSize[i]) + ".txt";
//        string bFilename = "b" + to_string(dataSize[i]) + ".txt";
//
//        init_A(xFilename);
//        init_P(bFilename);
//
//        // Performance measurement variables
//        double seconds;
//        long long head, tail, freq;
//
//        // Get frequency
//        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//
//        // 测试 f_parallel
//        QueryPerformanceCounter((LARGE_INTEGER*)&head);
//        //f_parallel();
//        f_serial();
//        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//        seconds = (tail - head) * 1000.0 / freq;
//        cout << "Data size: " << dataSize[i] << ", f_parallel() took " << seconds << " ms\n";
//
//        // Reset arrays for the next test
//        memset(Act, 0, sizeof(unsigned int) * (Num + 1));
//        memset(Pas, 0, sizeof(unsigned int) * (Num + 1));
//    }
//
//    return 0;
//}










//
//////////////////////////////////////////静态线程 + 信号量同步版本
//#include <pmmintrin.h>
//#include <xmmintrin.h>
//#include <emmintrin.h>
//#include <smmintrin.h>
//#include <tmmintrin.h>
//#include <nmmintrin.h>
//#include <iostream>
//#include <sstream>
//#include <fstream>
//#include <immintrin.h> 
//#include <windows.h>
//#include <cstring> 
//#include <vector>
//#include <semaphore.h>
//using namespace std;
//
//// 数据规模和数组将在 main 中动态定义
//unsigned int Act[37960][1188]; // 最大规模的数组预设大小
//unsigned int Pas[37960][1188];
//
//int Num;
//int pasNum;
//int lieNum;
//
//void init_A(const string& filename)
//{
//    ifstream infile(filename);
//    char fin[10000] = { 0 };
//    int index;
//
//    while (infile.getline(fin, sizeof(fin)))
//    {
//        stringstream line(fin);
//        int biaoji = 0;
//        unsigned int a;
//
//        while (line >> a)
//        {
//            if (biaoji == 0)
//            {
//                index = a;
//                biaoji = 1;
//            }
//            int k = a % 32;
//            int j = a / 32;
//
//            int temp = 1 << k;
//            Act[index][Num - 1 - j] |= temp; // Use bitwise OR to set the bit
//        }
//        Act[index][Num] = 1; // Mark this row as non-empty
//    }
//}
//
//void init_P(const string& filename)
//{
//    ifstream infile(filename);
//    char fin[10000] = { 0 };
//    int index = 0;
//
//    while (infile.getline(fin, sizeof(fin)))
//    {
//        stringstream line(fin);
//        int biaoji = 0;
//        unsigned int a;
//
//        while (line >> a)
//        {
//            if (biaoji == 0)
//            {
//                Pas[index][Num] = a;
//                biaoji = 1;
//            }
//            int k = a % 32;
//            int j = a / 32;
//
//            int temp = 1 << k;
//            Pas[index][Num - 1 - j] |= temp;
//        }
//        index++;
//    }
//}
//
//
//void f_serial() {
//    int i;
//    for (i = lieNum - 1; i >= 0; i -= 8) {
//        for (int j = 0; j < pasNum; j++) {
//            int idxPas = Pas[j][Num];
//            if (idxPas <= i && idxPas >= i - 7) {
//                int index = idxPas;
//                if (Act[index][Num] == 1) {
//                    for (int k = 0; k < Num; k++) {
//                        Pas[j][k] ^= Act[index][k];
//                    }
//                    int S_num = 0;
//                    for (int num = 0; num < Num; num++) {
//                        if (Pas[j][num] != 0) {
//                            for (int bit = 0; bit < 32; bit++) {
//                                if ((Pas[j][num] & (1 << bit)) != 0) {
//                                    S_num = bit + num * 32;
//                                    break;
//                                }
//                            }
//                            if (S_num != 0) break;
//                        }
//                    }
//                    Pas[j][Num] = S_num;
//                }
//                else {
//                    memcpy(Act[index], Pas[j], Num * sizeof(Pas[0][0]));
//                    Act[index][Num] = 1;
//                }
//            }
//        }
//    }
//}
//
//
//
//// 定义线程参数结构体
//struct ThreadData {
//    int startIdx; // 开始的索引
//    int endIdx;   // 结束的索引
//};
//
//// 全局信号量
//HANDLE sem_main;
//vector<HANDLE> sem_workerstart;
//vector<HANDLE> sem_workerend;
//
//DWORD WINAPI ThreadFuncSemaphore(LPVOID lpParam) {
//    ThreadData* pData = (ThreadData*)lpParam;
//    int startIdx = pData->startIdx;  // 线程处理的起始索引
//    int endIdx = pData->endIdx;      // 线程处理的结束索引
//
//    // 示例处理逻辑
//    for (int i = startIdx; i <= endIdx; i++) {
//        // 处理逻辑...
//        // 在示例中，我们将使用 f_serial 的逻辑作为处理逻辑的基础
//        for (int j = 0; j < pasNum; j++) {
//            unsigned int idxPas = Pas[j][Num];
//            if (idxPas <= i && idxPas >= i - 7) {
//                unsigned int index = idxPas;
//                unsigned int* pasPtr = Pas[j];
//                unsigned int* actPtr = Act[index];
//                if (actPtr[Num] == 1) {
//                    for (int k = 0; k < Num; k++) {
//                        pasPtr[k] ^= actPtr[k];
//                    }
                //    // 更新S_num等逻辑...
                //    int S_num = 0;
                //    for (int num = 0; num < Num; num++) {
                //        if (pasPtr[num] != 0) {
                //            for (int bit = 0; bit < 32; bit++) {
                //                if ((pasPtr[num] & (1 << bit)) != 0) {
                //                    S_num = bit + num * 32;
                //                    break;
                //                }
                //            }
                //            if (S_num != 0) break;
                //        }
                //    }
                //    pasPtr[Num] = S_num;
                //}
//                else {
//                    memcpy(actPtr, pasPtr, Num * sizeof(unsigned int));
//                    actPtr[Num] = 1;
//                }
//            }
//        }
//    }
//
//    // 信号量同步结束，通知主线程
//    ReleaseSemaphore(sem_main, 1, NULL);
//
//    return 0;
//}
//
//void f_parallel_semaphore() {
//    int threadCount = 7; // 或根据需要设置线程数量
//    std::vector<HANDLE> threads(threadCount);
//    std::vector<ThreadData> threadData(threadCount);
//
//    // 创建信号量
//    sem_main = CreateSemaphore(NULL, 0, threadCount, NULL);
//    sem_workerstart.resize(threadCount);
//    sem_workerend.resize(threadCount);
//    for (int i = 0; i < threadCount; ++i) {
//        sem_workerstart[i] = CreateSemaphore(NULL, 0, 1, NULL);
//        sem_workerend[i] = CreateSemaphore(NULL, 0, 1, NULL);
//    }
//
//    int partSize = lieNum / threadCount;
//    for (int i = 0; i < threadCount; ++i) {
//        threadData[i].startIdx = i * partSize;
//        threadData[i].endIdx = (i + 1) * partSize - 1;
//        if (i == threadCount - 1) threadData[i].endIdx = lieNum - 1; // 最后一个线程处理剩余部分
//
//        threads[i] = CreateThread(NULL, 0, ThreadFuncSemaphore, &threadData[i], 0, NULL);
//    }
//
//    // 等待所有线程完成
//    WaitForMultipleObjects(threadCount, threads.data(), TRUE, INFINITE);
//
//    // 清理资源
//    CloseHandle(sem_main);
//    for (int i = 0; i < threadCount; ++i) {
//        CloseHandle(sem_workerstart[i]);
//        CloseHandle(sem_workerend[i]);
//    }
//}
//
//
//int main()
//{
//    // Data sizes, corresponding pasNum and Num values
//    const int dataSize[] = { 130, 254, 562, 1011, 2362, 3799, 8399, 23045, 37960 };
//    const int pasNums[] = { 8, 53, 53, 263, 453, 1953, 4535, 14325, 14921 };
//    const int Nums[] = { 5, 8, 18, 32, 74, 119, 263, 721, 1187 };
//    const int dataCount = sizeof(dataSize) / sizeof(dataSize[0]);
//
//    // Loop over each data size for testing
//    for (int i = 0; i < dataCount; ++i)
//    {
//        lieNum = dataSize[i];
//        pasNum = pasNums[i];
//        Num = Nums[i];
//
//        string xFilename = "x" + to_string(dataSize[i]) + ".txt";
//        string bFilename = "b" + to_string(dataSize[i]) + ".txt";
//
//        init_A(xFilename);
//        init_P(bFilename);
//
//        // Performance measurement variables
//        double seconds;
//        long long head, tail, freq;
//
//        // Get frequency
//        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//
//        // 测试 f_parallel_semaphore
//        QueryPerformanceCounter((LARGE_INTEGER*)&head);
//        f_parallel_semaphore();
//        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//        seconds = (tail - head) * 1000.0 / freq;
//        cout << "Data size: " << dataSize[i] << ", f_parallel_semaphore() took " << seconds << " ms\n";
//
//        // Reset arrays for the next test
//        memset(Act, 0, sizeof(unsigned int) * (Num + 1));
//        memset(Pas, 0, sizeof(unsigned int) * (Num + 1));
//    }
//
//    return 0;
//}
//









//////////////////////////////////////////////静态线程 + 信号量同步版本 + 三重循环全部纳入线程函数

//#include <pmmintrin.h>
//#include <xmmintrin.h>
//#include <emmintrin.h>
//#include <smmintrin.h>
//#include <tmmintrin.h>
//#include <nmmintrin.h>
//#include <iostream>
//#include <sstream>
//#include <fstream>
//#include <immintrin.h> 
//#include <windows.h>
//#include <cstring> 
//#include <vector>
//#include <semaphore.h>
//using namespace std;
//
//// 数据规模和数组将在 main 中动态定义
//unsigned int Act[37960][1188]; // 最大规模的数组预设大小
//unsigned int Pas[37960][1188];
//
//int Num;
//int pasNum;
//int lieNum;
//
//void init_A(const string& filename)
//{
//    ifstream infile(filename);
//    char fin[10000] = { 0 };
//    int index;
//
//    while (infile.getline(fin, sizeof(fin)))
//    {
//        stringstream line(fin);
//        int biaoji = 0;
//        unsigned int a;
//
//        while (line >> a)
//        {
//            if (biaoji == 0)
//            {
//                index = a;
//                biaoji = 1;
//            }
//            int k = a % 32;
//            int j = a / 32;
//
//            int temp = 1 << k;
//            Act[index][Num - 1 - j] |= temp; // Use bitwise OR to set the bit
//        }
//        Act[index][Num] = 1; // Mark this row as non-empty
//    }
//}
//
//void init_P(const string& filename)
//{
//    ifstream infile(filename);
//    char fin[10000] = { 0 };
//    int index = 0;
//
//    while (infile.getline(fin, sizeof(fin)))
//    {
//        stringstream line(fin);
//        int biaoji = 0;
//        unsigned int a;
//
//        while (line >> a)
//        {
//            if (biaoji == 0)
//            {
//                Pas[index][Num] = a;
//                biaoji = 1;
//            }
//            int k = a % 32;
//            int j = a / 32;
//
//            int temp = 1 << k;
//            Pas[index][Num - 1 - j] |= temp;
//        }
//        index++;
//    }
//}
//
//
//void f_serial() {
//    int i;
//    for (i = lieNum - 1; i >= 0; i -= 8) {
//        for (int j = 0; j < pasNum; j++) {
//            int idxPas = Pas[j][Num];
//            if (idxPas <= i && idxPas >= i - 7) {
//                int index = idxPas;
//                if (Act[index][Num] == 1) {
//                    for (int k = 0; k < Num; k++) {
//                        Pas[j][k] ^= Act[index][k];
//                    }
//                    int S_num = 0;
//                    for (int num = 0; num < Num; num++) {
//                        if (Pas[j][num] != 0) {
//                            for (int bit = 0; bit < 32; bit++) {
//                                if ((Pas[j][num] & (1 << bit)) != 0) {
//                                    S_num = bit + num * 32;
//                                    break;
//                                }
//                            }
//                            if (S_num != 0) break;
//                        }
//                    }
//                    Pas[j][Num] = S_num;
//                }
//                else {
//                    memcpy(Act[index], Pas[j], Num * sizeof(Pas[0][0]));
//                    Act[index][Num] = 1;
//                }
//            }
//        }
//    }
//}
//
//
//
//// 定义线程参数结构体
//struct ThreadData {
//    int startIdx; // 开始的索引
//    int endIdx;   // 结束的索引
//};
//
//// 全局信号量
//HANDLE sem_main;
//vector<HANDLE> sem_workerstart;
//vector<HANDLE> sem_workerend;
//
//DWORD WINAPI ThreadFuncSemaphore(LPVOID lpParam) {
//    ThreadData* pData = (ThreadData*)lpParam;
//    int startIdx = pData->startIdx;  // 线程处理的起始索引
//    int endIdx = pData->endIdx;      // 线程处理的结束索引
//
//    // 示例处理逻辑
//    for (int i = startIdx; i <= endIdx; i++) {
//        // 处理逻辑...
//        // 在示例中，我们将使用 f_serial 的逻辑作为处理逻辑的基础
//        for (int j = 0; j < pasNum; j++) {
//            unsigned int idxPas = Pas[j][Num];
//            if (idxPas <= i && idxPas >= i - 7) {
//                unsigned int index = idxPas;
//                unsigned int* pasPtr = Pas[j];
//                unsigned int* actPtr = Act[index];
//                if (actPtr[Num] == 1) {
//                    for (int k = 0; k < Num; k++) {
//                        pasPtr[k] ^= actPtr[k];
//                    }
//                    // 更新S_num等逻辑...
//                    int S_num = 0;
//                    for (int num = 0; num < Num; num++) {
//                        if (pasPtr[num] != 0) {
//                            for (int bit = 0; bit < 32; bit++) {
//                                if ((pasPtr[num] & (1 << bit)) != 0) {
//                                    S_num = bit + num * 32;
//                                    break;
//                                }
//                            }
//                            if (S_num != 0) break;
//                        }
//                    }
//                    pasPtr[Num] = S_num;
//                }
//                else {
//                    memcpy(actPtr, pasPtr, Num * sizeof(unsigned int));
//                    actPtr[Num] = 1;
//                }
//            }
//        }
//    }
//
//    // 信号量同步结束，通知主线程
//    ReleaseSemaphore(sem_main, 1, NULL);
//
//    return 0;
//}
//
//void f_parallel_semaphore_full() {
//    int threadCount = 7; // 或根据需要设置线程数量
//    std::vector<HANDLE> threads(threadCount);
//    std::vector<ThreadData> threadData(threadCount);
//
//    // 创建信号量
//    sem_main = CreateSemaphore(NULL, 0, threadCount, NULL);
//    sem_workerstart.resize(threadCount);
//    sem_workerend.resize(threadCount);
//    for (int i = 0; i < threadCount; ++i) {
//        sem_workerstart[i] = CreateSemaphore(NULL, 0, 1, NULL);
//        sem_workerend[i] = CreateSemaphore(NULL, 0, 1, NULL);
//    }
//
//    // 分配工作范围并创建线程
//    int partSize = lieNum / threadCount;
//    for (int i = 0; i < threadCount; ++i) {
//        threadData[i].startIdx = i * partSize;
//        threadData[i].endIdx = (i + 1) * partSize - 1;
//        if (i == threadCount - 1) threadData[i].endIdx = lieNum - 1;
//
//        threads[i] = CreateThread(NULL, 0, ThreadFuncSemaphore, &threadData[i], 0, NULL);
//    }
//
//    // 等待所有线程完成
//    WaitForMultipleObjects(threadCount, threads.data(), TRUE, INFINITE);
//
//    // 清理资源
//    CloseHandle(sem_main);
//    for (int i = 0; i < threadCount; ++i) {
//        CloseHandle(sem_workerstart[i]);
//        CloseHandle(sem_workerend[i]);
//    }
//}
//
//
//int main()
//{
//    // Data sizes, corresponding pasNum and Num values
//    const int dataSize[] = { 130, 254, 562, 1011, 2362, 3799, 8399, 23045, 37960 };
//    const int pasNums[] = { 8, 53, 53, 263, 453, 1953, 4535, 14325, 14921 };
//    const int Nums[] = { 5, 8, 18, 32, 74, 119, 263, 721, 1187 };
//    const int dataCount = sizeof(dataSize) / sizeof(dataSize[0]);
//
//    // Loop over each data size for testing
//    for (int i = 0; i < dataCount; ++i)
//    {
//        lieNum = dataSize[i];
//        pasNum = pasNums[i];
//        Num = Nums[i];
//
//        string xFilename = "x" + to_string(dataSize[i]) + ".txt";
//        string bFilename = "b" + to_string(dataSize[i]) + ".txt";
//
//        init_A(xFilename);
//        init_P(bFilename);
//
//        // Performance measurement variables
//        double seconds;
//        long long head, tail, freq;
//
//        // Get frequency
//        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//
//        // 测试 f_parallel_semaphore
//        QueryPerformanceCounter((LARGE_INTEGER*)&head);
//        f_parallel_semaphore_full();
//        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//        seconds = (tail - head) * 1000.0 / freq;
//        cout << "Data size: " << dataSize[i] << ", f_parallel_semaphore() took " << seconds << " ms\n";
//
//        // Reset arrays for the next test
//        memset(Act, 0, sizeof(unsigned int) * (Num + 1));
//        memset(Pas, 0, sizeof(unsigned int) * (Num + 1));
//    }
//
//    return 0;
//}
//













//#include <pthread.h>
//#include <semaphore.h>
//#include <cstring>
//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <vector>
//#include <windows.h>
//#include <immintrin.h>
//using namespace std;
//
//// 数据规模和数组将在 main 中动态定义
//unsigned int Act[37960][1188]; // 最大规模的数组预设大小
//unsigned int Pas[37960][1188];
//
//int Num;
//int pasNum;
//int lieNum;
//
//void init_A(const string& filename)
//{
//    ifstream infile(filename);
//    char fin[10000] = { 0 };
//    int index;
//
//    while (infile.getline(fin, sizeof(fin)))
//    {
//        stringstream line(fin);
//        int biaoji = 0;
//        unsigned int a;
//
//        while (line >> a)
//        {
//            if (biaoji == 0)
//            {
//                index = a;
//                biaoji = 1;
//            }
//            int k = a % 32;
//            int j = a / 32;
//
//            int temp = 1 << k;
//            Act[index][Num - 1 - j] |= temp; // Use bitwise OR to set the bit
//        }
//        Act[index][Num] = 1; // Mark this row as non-empty
//    }
//}
//
//void init_P(const string& filename)
//{
//    ifstream infile(filename);
//    char fin[10000] = { 0 };
//    int index = 0;
//
//    while (infile.getline(fin, sizeof(fin)))
//    {
//        stringstream line(fin);
//        int biaoji = 0;
//        unsigned int a;
//
//        while (line >> a)
//        {
//            if (biaoji == 0)
//            {
//                Pas[index][Num] = a;
//                biaoji = 1;
//            }
//            int k = a % 32;
//            int j = a / 32;
//
//            int temp = 1 << k;
//            Pas[index][Num - 1 - j] |= temp;
//        }
//        index++;
//    }
//}
//
//
//void f_serial() {
//    int i;
//    for (i = lieNum - 1; i >= 0; i -= 8) {
//        for (int j = 0; j < pasNum; j++) {
//            int idxPas = Pas[j][Num];
//            if (idxPas <= i && idxPas >= i - 7) {
//                int index = idxPas;
//                if (Act[index][Num] == 1) {
//                    for (int k = 0; k < Num; k++) {
//                        Pas[j][k] ^= Act[index][k];
//                    }
//                    int S_num = 0;
//                    for (int num = 0; num < Num; num++) {
//                        if (Pas[j][num] != 0) {
//                            for (int bit = 0; bit < 32; bit++) {
//                                if ((Pas[j][num] & (1 << bit)) != 0) {
//                                    S_num = bit + num * 32;
//                                    break;
//                                }
//                            }
//                            if (S_num != 0) break;
//                        }
//                    }
//                    Pas[j][Num] = S_num;
//                }
//                else {
//                    memcpy(Act[index], Pas[j], Num * sizeof(Pas[0][0]));
//                    Act[index][Num] = 1;
//                }
//            }
//        }
//    }
//}
//
//
//
//struct ThreadDataBarrier {
//    int startIdx; 
//    int endIdx; 
//};
//pthread_barrier_t computation_barrier;
//void* ThreadFuncBarrier(void* lpParam) {
//    ThreadDataBarrier* pData = (ThreadDataBarrier*)lpParam;
//    int startIdx = pData->startIdx;  // 线程处理的起始索引
//    int endIdx = pData->endIdx;      // 线程处理的结束索引
//    for (int i = startIdx; i <= endIdx; i++) {
//        // 对Pas数组的每一列进行操作，这里的逻辑与f_serial中相同，但是是分段的
//        for (int j = 0; j < pasNum; j++) {
//            unsigned int idxPas = Pas[j][Num];
//            if (idxPas <= i && idxPas >= i - 7) {
//                unsigned int index = idxPas;
//                unsigned int* pasPtr = Pas[j];
//                unsigned int* actPtr = Act[index];
//                if (actPtr[Num] == 1) {
//                    for (int k = 0; k < Num; k++) {
//                        pasPtr[k] ^= actPtr[k];
//                    }
//                    int S_num = 0;
//                    for (int num = 0; num < Num; num++) {
//                        if (pasPtr[num] != 0) {
//                            for (int bit = 0; bit < 32; bit++) {
//                                if ((pasPtr[num] & (1 << bit)) != 0) {
//                                    S_num = bit + num * 32;
//                                    break;
//                                }
//                            }
//                            if (S_num != 0) break;
//                        }
//                    }
//                    pasPtr[Num] = S_num;
//                }
//                else {
//                    memcpy(actPtr, pasPtr, Num * sizeof(unsigned int));
//                    actPtr[Num] = 1;
//                }
//            }
//        }
//    }
//    pthread_barrier_wait(&computation_barrier);
//    return 0;
//}
//void f_parallel_barrier() {
//    const int threadCount = 7;
//    pthread_t threads[threadCount];
//    ThreadDataBarrier threadData[threadCount];
//    pthread_barrier_init(&computation_barrier, NULL, threadCount);
//    int partSize = lieNum / threadCount;
//    for (int i = 0; i < threadCount; ++i) {
//        threadData[i].startIdx = i * partSize;
//        threadData[i].endIdx = (i + 1) * partSize - 1;
//        if (i == threadCount - 1) threadData[i].endIdx = lieNum - 1; 
//
//        pthread_create(&threads[i], NULL, ThreadFuncBarrier, &threadData[i]);
//    }
//    for (int i = 0; i < threadCount; ++i) {
//        pthread_join(threads[i], NULL);
//    }
//    pthread_barrier_destroy(&computation_barrier);
//}
//
//int main()
//{
//    // Data sizes, corresponding pasNum and Num values
//    const int dataSize[] = { 130, 254, 562, 1011, 2362, 3799, 8399, 23045, 37960 };
//    const int pasNums[] = { 8, 53, 53, 263, 453, 1953, 4535, 14325, 14921 };
//    const int Nums[] = { 5, 8, 18, 32, 74, 119, 263, 721, 1187 };
//    const int dataCount = sizeof(dataSize) / sizeof(dataSize[0]);
//
//    // Loop over each data size for testing
//    for (int i = 0; i < dataCount; ++i)
//    {
//        lieNum = dataSize[i];
//        pasNum = pasNums[i];
//        Num = Nums[i];
//
//        string xFilename = "x" + to_string(dataSize[i]) + ".txt";
//        string bFilename = "b" + to_string(dataSize[i]) + ".txt";
//
//        init_A(xFilename);
//        init_P(bFilename);
//
//        // Performance measurement variables
//        double seconds;
//        long long head, tail, freq;
//
//        // Get frequency
//        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//
//        // 测试 f_parallel_semaphore
//        QueryPerformanceCounter((LARGE_INTEGER*)&head);
//        f_parallel_barrier();
//        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//        seconds = (tail - head) * 1000.0 / freq;
//        cout << "Data size: " << dataSize[i] << ", f_parallel_semaphore() took " << seconds << " ms\n";
//
//        // Reset arrays for the next test
//        memset(Act, 0, sizeof(unsigned int) * (Num + 1));
//        memset(Pas, 0, sizeof(unsigned int) * (Num + 1));
//    }
//
//    return 0;
//}