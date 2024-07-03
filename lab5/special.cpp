////#include <pmmintrin.h>
////#include <xmmintrin.h>
////#include <emmintrin.h>
////#include <smmintrin.h>
////#include <tmmintrin.h>
////#include <nmmintrin.h>
////#include <iostream>
////#include <sstream>
////#include <fstream>
////#include <immintrin.h> 
////#include <windows.h>
////#include <cstring> 
////using namespace std;
////
////// 数据规模和数组将在 main 中动态定义
////unsigned int Act[37960][1188]; // 最大规模的数组预设大小
////unsigned int Pas[37960][1188];
////
////int Num;
////int pasNum;
////int lieNum;
////
////void init_A(const string& filename)
////{
////    ifstream infile(filename);
////    char fin[10000] = { 0 };
////    int index;
////
////    while (infile.getline(fin, sizeof(fin)))
////    {
////        stringstream line(fin);
////        int biaoji = 0;
////        unsigned int a;
////
////        while (line >> a)
////        {
////            if (biaoji == 0)
////            {
////                index = a;
////                biaoji = 1;
////            }
////            int k = a % 32;
////            int j = a / 32;
////
////            int temp = 1 << k;
////            Act[index][Num - 1 - j] |= temp; // Use bitwise OR to set the bit
////        }
////        Act[index][Num] = 1; // Mark this row as non-empty
////    }
////}
////
////void init_P(const string& filename)
////{
////    ifstream infile(filename);
////    char fin[10000] = { 0 };
////    int index = 0;
////
////    while (infile.getline(fin, sizeof(fin)))
////    {
////        stringstream line(fin);
////        int biaoji = 0;
////        unsigned int a;
////
////        while (line >> a)
////        {
////            if (biaoji == 0)
////            {
////                Pas[index][Num] = a;
////                biaoji = 1;
////            }
////            int k = a % 32;
////            int j = a / 32;
////
////            int temp = 1 << k;
////            Pas[index][Num - 1 - j] |= temp;
////        }
////        index++;
////    }
////}
////
////
////void f_serial() {
////    int i;
////    for (i = lieNum - 1; i >= 0; i -= 8) {
////        for (int j = 0; j < pasNum; j++) {
////            int idxPas = Pas[j][Num];
////            if (idxPas <= i && idxPas >= i - 7) {
////                int index = idxPas;
////                if (Act[index][Num] == 1) {
////                    for (int k = 0; k < Num; k++) {
////                        Pas[j][k] ^= Act[index][k];
////                    }
////                    int S_num = 0;
////                    for (int num = 0; num < Num; num++) {
////                        if (Pas[j][num] != 0) {
////                            for (int bit = 0; bit < 32; bit++) {
////                                if ((Pas[j][num] & (1 << bit)) != 0) {
////                                    S_num = bit + num * 32;
////                                    break;
////                                }
////                            }
////                            if (S_num != 0) break;
////                        }
////                    }
////                    Pas[j][Num] = S_num;
////                }
////                else {
////                    memcpy(Act[index], Pas[j], Num * sizeof(Pas[0][0]));
////                    Act[index][Num] = 1;
////                }
////            }
////        }
////    }
////}
////
////
////int main()
////{
////    // Data sizes, corresponding pasNum and Num values
////    const int dataSize[] = { 130, 254, 562, 1011, 2362, 3799, 8399, 23045, 37960 };
////    const int pasNums[] = { 8, 53, 53, 263, 453, 1953, 4535, 14325, 14921 };
////    const int Nums[] = { 5, 8, 18, 32, 74, 119, 263, 721, 1187 };
////    const int dataCount = sizeof(dataSize) / sizeof(dataSize[0]);
////
////    // Loop over each data size for testing
////    for (int i = 0; i < dataCount; ++i)
////    {
////        lieNum = dataSize[i];
////        pasNum = pasNums[i];
////        Num = Nums[i];
////
////        string xFilename = "x" + to_string(dataSize[i]) + ".txt";
////        string bFilename = "b" + to_string(dataSize[i]) + ".txt";
////
////        init_A(xFilename);
////        init_P(bFilename);
////
////        // Performance measurement variables
////        double seconds;
////        long long head, tail, freq;
////
////        // Get frequency
////        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
////
////        // Test f_ordinary
////        QueryPerformanceCounter((LARGE_INTEGER*)&head);
////        f_serial();
////        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
////        seconds = (tail - head) * 1000.0 / freq;
////        cout << "Data size: " << dataSize[i] << ", f_serial() took " << seconds << " ms\n";
////
////        // Reset arrays for the next test
////        memset(Act, 0, sizeof(unsigned int) * (Num + 1));
////        memset(Pas, 0, sizeof(unsigned int) * (Num + 1));
////    }
////
////    return 0;
////}
////#include <mpi.h>
////#include <pmmintrin.h>
////#include <xmmintrin.h>
////#include <emmintrin.h>
////#include <smmintrin.h>
////#include <tmmintrin.h>
////#include <nmmintrin.h>
////#include <iostream>
////#include <sstream>
////#include <fstream>
////#include <immintrin.h> 
////#include <windows.h>
////#include <cstring> 
////using namespace std;
////
////// 数据规模和数组将在 main 中动态定义
////unsigned int Act[37960][1188]; // 最大规模的数组预设大小
////unsigned int Pas[37960][1188];
////
////int Num;
////int pasNum;
////int lieNum;
////
////void init_A(const string& filename)
////{
////    ifstream infile(filename);
////    char fin[10000] = { 0 };
////    int index;
////
////    while (infile.getline(fin, sizeof(fin)))
////    {
////        stringstream line(fin);
////        int biaoji = 0;
////        unsigned int a;
////
////        while (line >> a)
////        {
////            if (biaoji == 0)
////            {
////                index = a;
////                biaoji = 1;
////            }
////            int k = a % 32;
////            int j = a / 32;
////
////            int temp = 1 << k;
////            Act[index][Num - 1 - j] |= temp; // Use bitwise OR to set the bit
////        }
////        Act[index][Num] = 1; // Mark this row as non-empty
////    }
////}
////
////void init_P(const string& filename)
////{
////    ifstream infile(filename);
////    char fin[10000] = { 0 };
////    int index = 0;
////
////    while (infile.getline(fin, sizeof(fin)))
////    {
////        stringstream line(fin);
////        int biaoji = 0;
////        unsigned int a;
////
////        while (line >> a)
////        {
////            if (biaoji == 0)
////            {
////                Pas[index][Num] = a;
////                biaoji = 1;
////            }
////            int k = a % 32;
////            int j = a / 32;
////
////            int temp = 1 << k;
////            Pas[index][Num - 1 - j] |= temp;
////        }
////        index++;
////    }
////}
////
////
////void f_parallel(int rank, int size) {
////    int i;
////    for (i = lieNum - 1; i >= 0; i -= 8) {
////        for (int j = rank; j < pasNum; j += size) {
////            int idxPas = Pas[j][Num];
////            if (idxPas <= i && idxPas >= i - 7) {
////                int index = idxPas;
////                if (Act[index][Num] == 1) {
////                    for (int k = 0; k < Num; k++) {
////                        Pas[j][k] ^= Act[index][k];
////                    }
////                    int S_num = 0;
////                    for (int num = 0; num < Num; num++) {
////                        if (Pas[j][num] != 0) {
////                            for (int bit = 0; bit < 32; bit++) {
////                                if ((Pas[j][num] & (1 << bit)) != 0) {
////                                    S_num = bit + num * 32;
////                                    break;
////                                }
////                            }
////                            if (S_num != 0) break;
////                        }
////                    }
////                    Pas[j][Num] = S_num;
////                }
////                else {
////                    memcpy(Act[index], Pas[j], Num * sizeof(Pas[0][0]));
////                    Act[index][Num] = 1;
////                }
////            }
////        }
////    }
////}
////
////
////int main(int argc, char* argv[])
////{
////    MPI_Init(&argc, &argv);
////
////    int rank, size;
////    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
////    MPI_Comm_size(MPI_COMM_WORLD, &size);
////
////    if (rank == 0) {
////        // Data sizes, corresponding pasNum and Num values
////        const int dataSize[] = { 130, 254, 562, 1011, 2362, 3799, 8399, 23045, 37960 };
////        const int pasNums[] = { 8, 53, 53, 263, 453, 1953, 4535, 14325, 14921 };
////        const int Nums[] = { 5, 8, 18, 32, 74, 119, 263, 721, 1187 };
////        const int dataCount = sizeof(dataSize) / sizeof(dataSize[0]);
////
////        // Loop over each data size for testing
////        for (int i = 0; i < dataCount; ++i)
////        {
////            lieNum = dataSize[i];
////            pasNum = pasNums[i];
////            Num = Nums[i];
////
////            string xFilename = "x" + to_string(dataSize[i]) + ".txt";
////            string bFilename = "b" + to_string(dataSize[i]) + ".txt";
////
////            init_A(xFilename);
////            init_P(bFilename);
////
////            // Performance measurement variables
////            double seconds;
////            long long head, tail, freq;
////
////            // Get frequency
////            QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
////
////            // Test f_parallel
////            QueryPerformanceCounter((LARGE_INTEGER*)&head);
////            MPI_Bcast(&lieNum, 1, MPI_INT, 0, MPI_COMM_WORLD);
////            MPI_Bcast(&pasNum, 1, MPI_INT, 0, MPI_COMM_WORLD);
////            MPI_Bcast(&Num, 1, MPI_INT, 0, MPI_COMM_WORLD);
////            MPI_Bcast(&Act, sizeof(Act), MPI_BYTE, 0, MPI_COMM_WORLD);
////            MPI_Bcast(&Pas, sizeof(Pas), MPI_BYTE, 0, MPI_COMM_WORLD);
////
////            f_parallel(rank, size);
////
////            MPI_Reduce(MPI_IN_PLACE, Pas, pasNum * (Num + 1), MPI_UNSIGNED, MPI_BOR, 0, MPI_COMM_WORLD);
////            QueryPerformanceCounter((LARGE_INTEGER*)&tail);
////            if (rank == 0) {
////                seconds = (tail - head) * 1000.0 / freq;
////                cout << "Data size: " << dataSize[i] << ", f_parallel() took " << seconds << " ms\n";
////            }
////
////            // Reset arrays for the next test
////            memset(Act, 0, sizeof(unsigned int) * (Num + 1));
////            memset(Pas, 0, sizeof(unsigned int) * (Num + 1));
////        }
////    }
////    else {
////        MPI_Bcast(&lieNum, 1, MPI_INT, 0, MPI_COMM_WORLD);
////        MPI_Bcast(&pasNum, 1, MPI_INT, 0, MPI_COMM_WORLD);
////        MPI_Bcast(&Num, 1, MPI_INT, 0, MPI_COMM_WORLD);
////        MPI_Bcast(&Act, sizeof(Act), MPI_BYTE, 0, MPI_COMM_WORLD);
////        MPI_Bcast(&Pas, sizeof(Pas), MPI_BYTE, 0, MPI_COMM_WORLD);
////
////        f_parallel(rank, size);
////
////        MPI_Reduce(Pas, nullptr, pasNum * (Num + 1), MPI_UNSIGNED, MPI_BOR, 0, MPI_COMM_WORLD);
////    }
////
////    MPI_Finalize();
////    return 0;
////}
////
////#include <pmmintrin.h>
////#include <xmmintrin.h>
////#include <emmintrin.h>
////#include <smmintrin.h>
////#include <tmmintrin.h>
////#include <nmmintrin.h>
////#include <iostream>
////#include <sstream>
////#include <fstream>
////#include <immintrin.h>
////#include <windows.h>
////#include <cstring>
////#include <mpi.h> // 添加MPI头文件
////
////using namespace std;
////
////// 数据规模和数组将在 main 中动态定义
////unsigned int Act[37960][1188]; // 最大规模的数组预设大小
////unsigned int Pas[37960][1188];
////
////int Num;
////int pasNum;
////int lieNum;
////
////void init_A(const string& filename)
////{
////    ifstream infile(filename);
////    char fin[10000] = { 0 };
////    int index;
////
////    while (infile.getline(fin, sizeof(fin)))
////    {
////        stringstream line(fin);
////        int biaoji = 0;
////        unsigned int a;
////
////        while (line >> a)
////        {
////            if (biaoji == 0)
////            {
////                index = a;
////                biaoji = 1;
////            }
////            int k = a % 32;
////            int j = a / 32;
////
////            int temp = 1 << k;
////            Act[index][Num - 1 - j] |= temp; // Use bitwise OR to set the bit
////        }
////        Act[index][Num] = 1; // Mark this row as non-empty
////    }
////}
////
////void init_P(const string& filename)
////{
////    ifstream infile(filename);
////    char fin[10000] = { 0 };
////    int index = 0;
////
////    while (infile.getline(fin, sizeof(fin)))
////    {
////        stringstream line(fin);
////        int biaoji = 0;
////        unsigned int a;
////
////        while (line >> a)
////        {
////            if (biaoji == 0)
////            {
////                Pas[index][Num] = a;
////                biaoji = 1;
////            }
////            int k = a % 32;
////            int j = a / 32;
////
////            int temp = 1 << k;
////            Pas[index][Num - 1 - j] |= temp;
////        }
////        index++;
////    }
////}
////
////void f_parallel(int rank, int size) {
////    int i;
////    for (i = lieNum - 1 - rank; i >= 0; i -= size) {
////        for (int j = 0; j < pasNum; j++) {
////            int idxPas = Pas[j][Num];
////            if (idxPas <= i && idxPas >= i - 7) {
////                int index = idxPas;
////                if (Act[index][Num] == 1) {
////                    for (int k = 0; k < Num; k++) {
////                        Pas[j][k] ^= Act[index][k];
////                    }
////                    int S_num = 0;
////                    for (int num = 0; num < Num; num++) {
////                        if (Pas[j][num] != 0) {
////                            for (int bit = 0; bit < 32; bit++) {
////                                if ((Pas[j][num] & (1 << bit)) != 0) {
////                                    S_num = bit + num * 32;
////                                    break;
////                                }
////                            }
////                            if (S_num != 0) break;
////                        }
////                    }
////                    Pas[j][Num] = S_num;
////                }
////                else {
////                    memcpy(Act[index], Pas[j], Num * sizeof(Pas[0][0]));
////                    Act[index][Num] = 1;
////                }
////            }
////        }
////    }
////}
////
////int main(int argc, char* argv[])
////{
////    MPI_Init(&argc, &argv); // 初始化MPI环境
////    int rank, size;
////    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // 获取当前进程的秩
////    MPI_Comm_size(MPI_COMM_WORLD, &size); // 获取进程总数
////
////    // Data sizes, corresponding pasNum and Num values
////    const int dataSize[] = { 130, 254, 562, 1011, 2362, 3799, 8399, 23045, 37960 };
////    const int pasNums[] = { 8, 53, 53, 263, 453, 1953, 4535, 14325, 14921 };
////    const int Nums[] = { 5, 8, 18, 32, 74, 119, 263, 721, 1187 };
////    const int dataCount = sizeof(dataSize) / sizeof(dataSize[0]);
////
////    // Loop over each data size for testing
////    for (int i = 0; i < dataCount; ++i)
////    {
////        if (rank == 0) { // 只有主进程进行初始化
////            lieNum = dataSize[i];
////            pasNum = pasNums[i];
////            Num = Nums[i];
////
////            string xFilename = "x" + to_string(dataSize[i]) + ".txt";
////            string bFilename = "b" + to_string(dataSize[i]) + ".txt";
////
////            init_A(xFilename);
////            init_P(bFilename);
////        }
////
////        // 广播变量给所有进程
////        MPI_Bcast(&lieNum, 1, MPI_INT, 0, MPI_COMM_WORLD);
////        MPI_Bcast(&pasNum, 1, MPI_INT, 0, MPI_COMM_WORLD);
////        MPI_Bcast(&Num, 1, MPI_INT, 0, MPI_COMM_WORLD);
////        MPI_Bcast(Act, 37960 * 1188, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
////        MPI_Bcast(Pas, 37960 * 1188, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
////
////        // Performance measurement variables
////        double seconds;
////        long long head, tail, freq;
////
////        // Get frequency
////        if (rank == 0) {
////            QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
////            QueryPerformanceCounter((LARGE_INTEGER*)&head);
////        }
////
////        // Test f_parallel
////        f_parallel(rank, size);
////
////        // 创建用于收集结果的缓冲区
////        unsigned int* recvbuf = nullptr;
////        if (rank == 0) {
////            recvbuf = new unsigned int[pasNum * (Num + 1) * size];
////        }
////
////        // Gather results to the root process
////        MPI_Gather(Pas, pasNum * (Num + 1), MPI_UNSIGNED, recvbuf, pasNum * (Num + 1), MPI_UNSIGNED, 0, MPI_COMM_WORLD);
////
////        if (rank == 0) {
////            QueryPerformanceCounter((LARGE_INTEGER*)&tail);
////            seconds = (tail - head) * 1000.0 / freq;
////            cout << "Data size: " << dataSize[i] << ", f_parallel() took " << seconds << " ms\n";
////
////            // 释放接收缓冲区
////            delete[] recvbuf;
////        }
////
////        // Reset arrays for the next test
////        memset(Act, 0, sizeof(unsigned int) * (Num + 1));
////        memset(Pas, 0, sizeof(unsigned int) * (Num + 1));
////    }
////
////    MPI_Finalize(); // 结束MPI环境
////    return 0;
////}
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
//#include <mpi.h> // 添加MPI头文件
//
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
//void f_parallel(int rank, int size) {
//    int i;
//    for (i = lieNum - 1 - rank; i >= 0; i -= size) {
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
//int main(int argc, char* argv[])
//{
//    MPI_Init(&argc, &argv); // 初始化MPI环境
//    int rank, size;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // 获取当前进程的秩
//    MPI_Comm_size(MPI_COMM_WORLD, &size); // 获取进程总数
//
//    // Data sizes, corresponding pasNum and Num values
//    const int dataSize[] = { 130, 254, 562, 1011, 2362, 3799, 8399, 23045, 37960 };
//    const int pasNums[] = { 8, 53, 53, 263, 453, 1953, 4535, 14325, 14921 };
//    const int Nums[] = { 5, 8, 18, 32, 74, 119, 263, 721, 1187 };
//    const int dataCount = sizeof(dataSize) / sizeof(dataSize[0]);
//
//    // Loop over each data size for testing
//    for (int i = 0; i < dataCount; ++i)
//    {
//        if (rank == 0) { // 只有主进程进行初始化
//            lieNum = dataSize[i];
//            pasNum = pasNums[i];
//            Num = Nums[i];
//
//            string xFilename = "x" + to_string(dataSize[i]) + ".txt";
//            string bFilename = "b" + to_string(dataSize[i]) + ".txt";
//
//            init_A(xFilename);
//            init_P(bFilename);
//        }
//
//        // 广播变量给所有进程
//        MPI_Bcast(&lieNum, 1, MPI_INT, 0, MPI_COMM_WORLD);
//        MPI_Bcast(&pasNum, 1, MPI_INT, 0, MPI_COMM_WORLD);
//        MPI_Bcast(&Num, 1, MPI_INT, 0, MPI_COMM_WORLD);
//        MPI_Bcast(Act, 37960 * 1188, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//        MPI_Bcast(Pas, 37960 * 1188, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//
//        // Performance measurement variables
//        double seconds;
//        long long head, tail, freq;
//
//        // Get frequency
//        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
//
//        // Test f_serial
//        if (rank == 0) {
//            QueryPerformanceCounter((LARGE_INTEGER*)&head);
//            f_serial();
//            QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//            seconds = (tail - head) * 1000.0 / freq;
//            cout << "Data size: " << dataSize[i] << ", f_serial() took " << seconds << " ms\n";
//        }
//
//        // Reset arrays for the next test
//        memset(Act, 0, sizeof(unsigned int) * (Num + 1));
//        memset(Pas, 0, sizeof(unsigned int) * (Num + 1));
//
//        if (rank == 0) { // 只有主进程重新初始化数据
//            init_A("x" + to_string(dataSize[i]) + ".txt");
//            init_P("b" + to_string(dataSize[i]) + ".txt");
//        }
//
//        // 广播新的初始化数据给所有进程
//        MPI_Bcast(Act, 37960 * 1188, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//        MPI_Bcast(Pas, 37960 * 1188, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//
//        // Test f_parallel
//        if (rank == 0) {
//            QueryPerformanceCounter((LARGE_INTEGER*)&head);
//        }
//
//        f_parallel(rank, size);
//
//        // 创建用于收集结果的缓冲区
//        unsigned int* recvbuf = nullptr;
//        if (rank == 0) {
//            recvbuf = new unsigned int[pasNum * (Num + 1) * size];
//        }
//
//        // Gather results to the root process
//        MPI_Gather(Pas, pasNum * (Num + 1), MPI_UNSIGNED, recvbuf, pasNum * (Num + 1), MPI_UNSIGNED, 0, MPI_COMM_WORLD);
//
//        if (rank == 0) {
//            QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//            seconds = (tail - head) * 1000.0 / freq;
//                cout << "Data size: " << dataSize[i] << ", f_parallel() took " << seconds << " ms\n";
//        }
//
//        // 释放用于收集结果的缓冲区
//        if (rank == 0) {
//            delete[] recvbuf;
//        }
//
//        // 重新初始化数组以备下次测试
//        memset(Act, 0, sizeof(unsigned int) * (Num + 1));
//        memset(Pas, 0, sizeof(unsigned int) * (Num + 1));
//    }
//
//    MPI_Finalize(); // 结束MPI环境
//    return 0;
//}








/////////////////////////////////////////////////////////////001
//#include<iostream>
//#include<algorithm>
//#include<fstream>
//#include<sstream>
//#include<omp.h>
//#include"mpi.h"
//using namespace std;
//const int col = 23100, elinenum = 14325, num_thread = 8, mpisize = 8;//列数、被消元行数
//bool parallel = 1;
//int isupgrade;
//int tmp = 0;
//const int bytenum = (col - 1) / 32 + 1;   //每个实例中的byte型数组数
//
//
//class bitmatrix {
//public:
//	int mycol;    //首项
//	int* mybyte;
//	bitmatrix() {    //初始化
//		mycol = -1;
//		mybyte = new int[bytenum];
//		for (int i = 0; i < bytenum; i++)
//			mybyte[i] = 0;
//	}
//	bool isnull() {  //判断当前行是否为空行
//		if (mycol == -1)return 1;
//		return 0;
//	}
//	void insert(int x) { //数据读入
//		if (mycol == -1)mycol = x;
//		int a = x / 32, b = x % 32;
//		mybyte[a] |= (1 << b);
//	}
//	void doxor(bitmatrix b) {  //两行做异或操作，由于结果留在本实例中，只有被消元行能执行这一操作,且异或操作后要更新首项
//		for (int i = 0; i < bytenum; i++)
//			mybyte[i] ^= b.mybyte[i];
//		for (int i = bytenum - 1; i >= 0; i--)
//			for (int j = 31; j >= 0; j--)
//				if ((mybyte[i] & (1 << j)) != 0) {
//					mycol = i * 32 + j;
//					return;
//				}
//		mycol = -1;
//	}
//};
//bitmatrix* eliminer = new bitmatrix[col], * eline = new bitmatrix[elinenum], * pass = new bitmatrix[mpisize];
//void readdata() {
//	ifstream ifs;
//	ifs.open("x23045.txt");  //消元子
//	string temp;
//	while (getline(ifs, temp)) {
//		istringstream ss(temp);
//		int x;
//		int trow = 0;
//		while (ss >> x) {
//			if (!trow)trow = x;    //第一个读入元素代表行号
//			eliminer[trow].insert(x);
//		}
//	}
//	ifs.close();
//	ifstream ifs2;
//	ifs2.open("b23045.txt");     //被消元行,读入方式与消元子不同
//	int trow = 0;
//	while (getline(ifs2, temp)) {
//		istringstream ss(temp);
//		int x;
//		while (ss >> x) {
//			eline[trow].insert(x);
//		}
//		trow++;
//	}
//	ifs2.close();
//}
//void printres() { //打印结果
//	for (int i = 0; i < elinenum; i++) {
//		if (eline[i].isnull()) { puts(""); continue; }   //空行的特殊情况
//		for (int j = bytenum - 1; j >= 0; j--) {
//			for (int k = 31; k >= 0; k--)
//				if ((eline[i].mybyte[j] & (1 << k)) != 0) {     //一个错误调了半小时，谨记当首位为1时>>不等于除法！
//					printf("%d ", j * 32 + k);
//				}
//		}
//		puts("");
//	}
//}
//void dowork() {  
//	int rank;
//	MPI_Status status;
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	//获取当前进程号
//	int r1 = rank * (elinenum / mpisize), r2 = (rank == mpisize - 1) ? elinenum - 1 : (rank + 1) * (elinenum / mpisize) - 1;
//	double st = MPI_Wtime();	 //计时开始
//	int i, j, k;
//	for (i = col - 1; i >= 0; i--) {
//		if (!eliminer[i].isnull()) {
//			for (j = r1; j <= r2; j++) {
//				if (eline[j].mycol == i)
//					eline[j].doxor(eliminer[i]);
//			}
//		}
//		else {
//			//cout << rank<<" "<<2 << endl;
//			isupgrade = -1;
//			int t = -1;
//			if (rank != 0) {
//				for (k = r1; k <= r2; k++)
//					if (eline[k].mycol == i) {
//						pass[rank] = eline[k];
//						t = k;
//						MPI_Send(&t, 1, MPI_INT, 0, k + elinenum * 4 + 3, MPI_COMM_WORLD);
//						MPI_Send(&pass[rank].mybyte[0], bytenum, MPI_INT, 0, k + 3, MPI_COMM_WORLD);	//各进程内消元完毕的被消元行传回0号进程
//						MPI_Send(&pass[rank].mycol, 1, MPI_INT, 0, k + elinenum + 3, MPI_COMM_WORLD);
//						break;
//					}
//				if (k > r2) {
//					MPI_Send(&t, 1, MPI_INT, 0, k + elinenum * 4 + 3, MPI_COMM_WORLD);
//					MPI_Send(&pass[rank].mybyte[0], bytenum, MPI_INT, 0, k + 3, MPI_COMM_WORLD);	//各进程内消元完毕的被消元行传回0号进程
//					MPI_Send(&pass[rank].mycol, 1, MPI_INT, 0, k + elinenum + 3, MPI_COMM_WORLD);
//				}
//			}
//			else {
//				for (k = mpisize - 1; k > 0; k--) {
//					MPI_Recv(&t, 1, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//					MPI_Recv(&pass[k].mybyte[0], bytenum, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//					MPI_Recv(&pass[k].mycol, 1, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//					if (t != -1) {
//						eliminer[i] = pass[k];
//						isupgrade = t;
//					}
//				}
//				for (k = r1; k <= r2; k++)
//					if (eline[k].mycol == i) {
//						eliminer[i] = eline[k];
//						isupgrade = k;
//						break;
//					}
//				for (k = 1; k < mpisize; k++) {
//					int t1 = k * (elinenum / mpisize), t2 = (k == mpisize - 1) ? elinenum - 1 : (k + 1) * (elinenum / mpisize) - 1;
//					MPI_Send(&isupgrade, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
//					if (isupgrade != -1 && t2 >= isupgrade) {
//						MPI_Send(&eliminer[i].mybyte[0], bytenum, MPI_INT, k, 1, MPI_COMM_WORLD);
//						MPI_Send(&eliminer[i].mycol, 1, MPI_INT, k, 2, MPI_COMM_WORLD);
//					}
//				}
//			}
//
//			//MPI_Barrier(MPI_COMM_WORLD);
//			if (rank != 0) {
//				MPI_Recv(&isupgrade, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
//				if (isupgrade != -1 && r2 >= isupgrade) {
//					MPI_Recv(&eliminer[i].mybyte[0], bytenum, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
//					MPI_Recv(&eliminer[i].mycol, bytenum, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
//				}
//			}
//			if (isupgrade != -1 && r2 >= isupgrade)
//				for (j = r1; j <= r2; j++) {
//					if (eline[j].mycol == i && j != isupgrade)
//						eline[j].doxor(eliminer[i]);
//				}
//		}
//	}
//
//	if (rank != 0)
//		for (k = r1; k <= r2; k++) {
//			MPI_Send(&eline[k].mybyte[0], bytenum, MPI_INT, 0, k + 3 + elinenum * 2, MPI_COMM_WORLD);	//各进程内消元完毕的被消元行传回0号进程
//			MPI_Send(&eline[k].mycol, 1, MPI_INT, 0, k + 3 + elinenum * 3, MPI_COMM_WORLD);
//		}
//	else
//		for (k = 1; k < mpisize; k++) {
//			int t1 = k * (elinenum / mpisize), t2 = (k == mpisize - 1) ? elinenum - 1 : (k + 1) * (elinenum / mpisize) - 1;
//			for (int q = t1; q <= t2; q++) {
//				MPI_Recv(&eline[q].mybyte[0], bytenum, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//				MPI_Recv(&eline[q].mycol, 1, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
//			}
//		}
//	double ed = MPI_Wtime();	 //计时结束
//	if (rank == 0) {	//只有0号进程中有最终结果
//		printf("cost time:%.4lf\n", ed - st);
//		//printres();
//	}
//}
//
//int main() {
//	MPI_Init(0, 0);
//	readdata();
//	dowork();
//	MPI_Finalize();
//	return 0;
//}




#include<iostream>
#include<algorithm>
#include<fstream>
#include<sstream>
#include<omp.h>
#include"mpi.h"
using namespace std;
const int col = 23100, elinenum = 14325, num_thread = 8, mpisize = 8;//列数、被消元行数
bool parallel = 1;
int isupgrade;
int tmp = 0;
const int bytenum = (col - 1) / 32 + 1;   //每个实例中的byte型数组数


class bitmatrix {
public:
	int mycol;    //首项
	int* mybyte;
	bitmatrix() {    //初始化
		mycol = -1;
		mybyte = new int[bytenum];
		for (int i = 0; i < bytenum; i++)
			mybyte[i] = 0;
	}
	bool isnull() {  //判断当前行是否为空行
		if (mycol == -1)return 1;
		return 0;
	}
	void insert(int x) { //数据读入
		if (mycol == -1)mycol = x;
		int a = x / 32, b = x % 32;
		mybyte[a] |= (1 << b);
	}
	void doxor(bitmatrix b) {  //两行做异或操作，由于结果留在本实例中，只有被消元行能执行这一操作,且异或操作后要更新首项
		for (int i = 0; i < bytenum; i++)
			mybyte[i] ^= b.mybyte[i];
		for (int i = bytenum - 1; i >= 0; i--)
			for (int j = 31; j >= 0; j--)
				if ((mybyte[i] & (1 << j)) != 0) {
					mycol = i * 32 + j;
					return;
				}
		mycol = -1;
	}
};
bitmatrix* eliminer = new bitmatrix[col], * eline = new bitmatrix[elinenum], * pass = new bitmatrix[mpisize];
void readdata() {
	ifstream ifs;
	ifs.open("x23045.txt");  //消元子
	string temp;
	while (getline(ifs, temp)) {
		istringstream ss(temp);
		int x;
		int trow = 0;
		while (ss >> x) {
			if (!trow)trow = x;    //第一个读入元素代表行号
			eliminer[trow].insert(x);
		}
	}
	ifs.close();
	ifstream ifs2;
	ifs2.open("b23045.txt");     //被消元行,读入方式与消元子不同
	int trow = 0;
	while (getline(ifs2, temp)) {
		istringstream ss(temp);
		int x;
		while (ss >> x) {
			eline[trow].insert(x);
		}
		trow++;
	}
	ifs2.close();
}
void printres() { //打印结果
	for (int i = 0; i < elinenum; i++) {
		if (eline[i].isnull()) { puts(""); continue; }   //空行的特殊情况
		for (int j = bytenum - 1; j >= 0; j--) {
			for (int k = 31; k >= 0; k--)
				if ((eline[i].mybyte[j] & (1 << k)) != 0) {     //一个错误调了半小时，谨记当首位为1时>>不等于除法！
					printf("%d ", j * 32 + k);
				}
		}
		puts("");
	}
}

void dowork() {  //串行消元--消元子->被消元行
	int rank;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);	//获取当前进程号
	int r1 = rank * (elinenum / mpisize), r2 = (rank == mpisize - 1) ? elinenum - 1 : (rank + 1) * (elinenum / mpisize) - 1;
	double st = MPI_Wtime();	 //计时开始
	int i, j, k;
#pragma omp parallel if(parallel),num_threads(num_thread),private(i,j)
	for (i = col - 1; i >= 0; i--) {
		if (!eliminer[i].isnull()) {
#pragma omp for 
			for (j = r1; j <= r2; j++) {
				if (eline[j].mycol == i)
					eline[j].doxor(eliminer[i]);
			}
		}
		else {
			isupgrade = -1;
			int t = -1;
			if (rank != 0) {
#pragma omp single
				for (k = r1; k <= r2; k++)
					if (eline[k].mycol == i) {
						pass[rank] = eline[k];
						t = k;
						MPI_Send(&t, 1, MPI_INT, 0, k + elinenum * 4 + 3, MPI_COMM_WORLD);
						MPI_Send(&pass[rank].mybyte[0], bytenum, MPI_INT, 0, k + 3, MPI_COMM_WORLD);
						MPI_Send(&pass[rank].mycol, 1, MPI_INT, 0, k + elinenum + 3, MPI_COMM_WORLD);
						break;
					}
				if (k > r2) {
					MPI_Send(&t, 1, MPI_INT, 0, k + elinenum * 4 + 3, MPI_COMM_WORLD);
					MPI_Send(&pass[rank].mybyte[0], bytenum, MPI_INT, 0, k + 3, MPI_COMM_WORLD);	//各进程内消元完毕的被消元行传回0号进程
					MPI_Send(&pass[rank].mycol, 1, MPI_INT, 0, k + elinenum + 3, MPI_COMM_WORLD);
				}
			}
			else {
#pragma omp single
				for (k = mpisize - 1; k > 0; k--) {
					MPI_Recv(&t, 1, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
					MPI_Recv(&pass[k].mybyte[0], bytenum, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
					MPI_Recv(&pass[k].mycol, 1, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
					if (t != -1) {
						eliminer[i] = pass[k];
						isupgrade = t;

					}
				}
#pragma omp single
				for (k = r1; k <= r2; k++)
					if (eline[k].mycol == i) {
						eliminer[i] = eline[k];
						isupgrade = k;
						break;
					}
#pragma omp for
				for (k = 1; k < mpisize; k++) {
					int t1 = k * (elinenum / mpisize), t2 = (k == mpisize - 1) ? elinenum - 1 : (k + 1) * (elinenum / mpisize) - 1;
					MPI_Send(&isupgrade, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
					if (isupgrade != -1 && t2 >= isupgrade) {
						MPI_Send(&eliminer[i].mybyte[0], bytenum, MPI_INT, k, 1, MPI_COMM_WORLD);
						MPI_Send(&eliminer[i].mycol, 1, MPI_INT, k, 2, MPI_COMM_WORLD);
					}
				}
			}

			//MPI_Barrier(MPI_COMM_WORLD);
			if (rank != 0) {
				MPI_Recv(&isupgrade, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
				if (isupgrade != -1 && r2 >= isupgrade) {
					MPI_Recv(&eliminer[i].mybyte[0], bytenum, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
					MPI_Recv(&eliminer[i].mycol, bytenum, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
				}
			}
			if (isupgrade != -1 && r2 >= isupgrade) {
#pragma omp for
				for (j = r1; j <= r2; j++) {
					if (eline[j].mycol == i && j != isupgrade)
						eline[j].doxor(eliminer[i]);
				}
			}
		}
	}
	if (rank != 0)
		for (k = r1; k <= r2; k++) {
			MPI_Send(&eline[k].mybyte[0], bytenum, MPI_INT, 0, k + 3 + elinenum * 2, MPI_COMM_WORLD);	//各进程内消元完毕的被消元行传回0号进程
			MPI_Send(&eline[k].mycol, 1, MPI_INT, 0, k + 3 + elinenum * 3, MPI_COMM_WORLD);
		}
	else
		for (k = 1; k < mpisize; k++) {
			int t1 = k * (elinenum / mpisize), t2 = (k == mpisize - 1) ? elinenum - 1 : (k + 1) * (elinenum / mpisize) - 1;
			for (int q = t1; q <= t2; q++) {
				MPI_Recv(&eline[q].mybyte[0], bytenum, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				MPI_Recv(&eline[q].mycol, 1, MPI_INT, k, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
		}
	double ed = MPI_Wtime();	 //计时结束
	if (rank == 0) {	//只有0号进程中有最终结果
		printf("cost time:%.4lf\n", ed - st);
		//printres();
	}
}

int main() {
	MPI_Init(0, 0);
	readdata();
	dowork();
	MPI_Finalize();
	return 0;
}