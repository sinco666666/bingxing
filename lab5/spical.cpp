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
//void f_ordinary()
//{
//    int i;
//    for (i = lieNum - 1; i - 8 >= -1; i -= 8)
//    {
//        //每轮处理8个消元子，范围：首项在 i-7 到 i
//        for (int j = 0; j < pasNum; j++)
//        {
//
//            //看4535个被消元行有没有首项在此范围内的
//            while (Pas[j][Num] <= i && Pas[j][Num] >= i - 7)
//            {
//                int index = Pas[j][Num];
//                if (Act[index][Num] == 1)//消元子不为空
//                {
//                    //Pas[j][]和Act[（Pas[j][18]）][]做异或
//                    for (int k = 0; k < Num; k++)
//                    {
//                        Pas[j][k] = Pas[j][k] ^ Act[index][k];
//                    }
//
//                    //更新Pas[j][18]存的首项值
//                    //做完异或之后继续找这个数的首项，存到Pas[j][18]，若还在范围里会继续while循环
//                    //找异或之后Pas[j][ ]的首项
//                    int num = 0, S_num = 0;
//                    for (num = 0; num < Num; num++)
//                    {
//                        if (Pas[j][num] != 0)
//                        {
//                            unsigned int temp = Pas[j][num];
//                            while (temp != 0)
//                            {
//                                temp = temp >> 1;
//                                S_num++;
//                            }
//                            S_num += num * 32;
//                            break;
//                        }
//                    }
//                    Pas[j][Num] = S_num - 1;
//
//                }
//                else//消元子为空
//                {
//                    //Pas[j][]来补齐消元子
//                    for (int k = 0; k < Num; k++)
//                        Act[index][k] = Pas[j][k];
//
//                    Act[index][Num] = 1;//设置消元子非空
//                    break;
//                }
//
//            }
//        }
//    }
//
//    for (i = i + 8; i >= 0; i--)
//    {
//        //每轮处理1个消元子，范围：首项等于i
//
//        for (int j = 0; j < pasNum; j++)
//        {
//            //看53个被消元行有没有首项等于i的
//            while (Pas[j][Num] == i)
//            {
//                if (Act[i][Num] == 1)//消元子不为空
//                {
//                    //Pas[j][]和Act[i][]做异或
//                    for (int k = 0; k < Num; k++)
//                    {
//                        Pas[j][k] = Pas[j][k] ^ Act[i][k];
//                    }
//
//                    //更新Pas[j][18]存的首项值
//                    //做完异或之后继续找这个数的首项，存到Pas[j][18]，若还在范围里会继续while循环
//                    //找异或之后Pas[j][ ]的首项
//                    int num = 0, S_num = 0;
//                    for (num = 0; num < Num; num++)
//                    {
//                        if (Pas[j][num] != 0)
//                        {
//                            unsigned int temp = Pas[j][num];
//                            while (temp != 0)
//                            {
//                                temp = temp >> 1;
//                                S_num++;
//                            }
//                            S_num += num * 32;
//                            break;
//                        }
//                    }
//                    Pas[j][Num] = S_num - 1;
//
//                }
//                else//消元子为空
//                {
//                    //Pas[j][]来补齐消元子
//                    for (int k = 0; k < Num; k++)
//                        Act[i][k] = Pas[j][k];
//
//                    Act[i][Num] = 1;//设置消元子非空
//                    break;
//                }
//            }
//        }
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
//void f_cache_optimized() {
//    for (int i = lieNum - 1; i >= 0; i -= 8) {
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
//}
//
//__m128 vakj;
//__m128 vaij;
//
//void f_sse() {
//    int i;
//    for (i = lieNum - 1; i >= 0; i -= 8) {
//        for (int j = 0; j < pasNum; j++) {
//            int idxPas = Pas[j][Num];
//            if (idxPas <= i && idxPas >= i - 7) {
//                int index = idxPas;
//                if (Act[index][Num] == 1) {
//                    for (int k = 0; k + 4 <= Num; k += 4) {
//                        vakj = _mm_loadu_ps((float*)&Pas[j][k]);
//                        vaij = _mm_loadu_ps((float*)&Act[index][k]);
//                        vakj = _mm_xor_ps(vakj, vaij);
//                        _mm_storeu_ps((float*)&Pas[j][k], vakj);
//                    }
//                    for (int k = (Num / 4) * 4; k < Num; k++) {
//                        Pas[j][k] ^= Act[index][k];
//                    }
//                    int S_num = 0;
//                    for (int num = 0; num < Num; num++) {
//                        unsigned int temp = Pas[j][num];
//                        if (temp != 0) {
//                            unsigned long index;
//                            _BitScanForward(&index, temp);
//                            S_num = index + num * 32;
//                            break;
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
//void f_sse_a() {
//    int i;
//    for (i = lieNum - 1; i >= 0; i -= 8) {
//        for (int j = 0; j < pasNum; j++) {
//            int idxPas = Pas[j][Num];
//            if (idxPas <= i && idxPas >= i - 7) {
//                int index = idxPas;
//                if (Act[index][Num] == 1) {
//                    for (int k = 0; k + 4 <= Num; k += 4) {
//                        vakj = _mm_load_ps((float*)&Pas[j][k]); // Aligned load
//                        vaij = _mm_load_ps((float*)&Act[index][k]); // Aligned load
//                        vakj = _mm_xor_ps(vakj, vaij);
//                        _mm_storeu_ps((float*)&Pas[j][k], vakj);
//                    }
//                    for (int k = (Num / 4) * 4; k < Num; k++) {
//                        Pas[j][k] ^= Act[index][k];
//                    }
//                    int S_num = 0;
//                    for (int num = 0; num < Num; num++) {
//                        unsigned int temp = Pas[j][num];
//                        if (temp != 0) {
//                            unsigned long index;
//                            _BitScanForward(&index, temp);
//                            S_num = index + num * 32;
//                            break;
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
//__m256 va_Pas2;
//__m256 va_Act2;
//
//void f_avx256() {
//    int i;
//    for (i = lieNum - 1; i >= 0; i -= 8) {
//        for (int j = 0; j < pasNum; j++) {
//            int idxPas = Pas[j][Num];
//            if (idxPas <= i && idxPas >= i - 7) {
//                int index = idxPas;
//                if (Act[index][Num] == 1) {
//                    for (int k = 0; k + 8 <= Num; k += 8) {
//                        va_Pas2 = _mm256_loadu_ps((float*)&Pas[j][k]);
//                        va_Act2 = _mm256_loadu_ps((float*)&Act[index][k]);
//                        va_Pas2 = _mm256_xor_ps(va_Pas2, va_Act2);
//                        _mm256_storeu_ps((float*)&Pas[j][k], va_Pas2);
//                    }
//                    for (int k = (Num / 8) * 8; k < Num; k++) {
//                        Pas[j][k] ^= Act[index][k];
//                    }
//
//                    int S_num = 0;
//                    for (int num = 0; num < Num; num++) {
//                        unsigned int temp = Pas[j][num];
//                        if (temp != 0) {
//                            unsigned long index;
//                            _BitScanForward(&index, temp);
//                            S_num = index + num * 32;
//                            break;
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
//void f_avx256_a() {
//    int i;
//    for (i = lieNum - 1; i >= 0; i -= 8) {
//        for (int j = 0; j < pasNum; j++) {
//            int idxPas = Pas[j][Num];
//            if (idxPas <= i && idxPas >= i - 7) {
//                int index = idxPas;
//                if (Act[index][Num] == 1) {
//                    for (int k = 0; k + 8 <= Num; k += 8) {
//                        va_Pas2 = _mm256_load_ps((float*)&Pas[j][k]); // Aligned load
//                        va_Act2 = _mm256_load_ps((float*)&Act[index][k]); // Aligned load
//                        va_Pas2 = _mm256_xor_ps(va_Pas2, va_Act2);
//                        _mm256_storeu_ps((float*)&Pas[j][k], va_Pas2);
//                    }
//                    for (int k = (Num / 8) * 8; k < Num; k++) {
//                        Pas[j][k] ^= Act[index][k];
//                    }
//
//                    int S_num = 0;
//                    for (int num = 0; num < Num; num++) {
//                        unsigned int temp = Pas[j][num];
//                        if (temp != 0) {
//                            unsigned long index;
//                            _BitScanForward(&index, temp);
//                            S_num = index + num * 32;
//                            break;
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
//        // Test f_ordinary
//        QueryPerformanceCounter((LARGE_INTEGER*)&head);
//        //f_cache_optimized();
//        f_ordinary();
//        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//        seconds = (tail - head) * 1000.0 / freq;
//        cout << "Data size: " << dataSize[i] << ",f_ordinary() took " << seconds << " ms\n";
//
//        // Reset arrays for the next test
//        memset(Act, 0, sizeof(unsigned int) * (Num + 1));
//        memset(Pas, 0, sizeof(unsigned int) * (Num + 1));
//    }
//
//    return 0;
//}




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
//void f_ordinary()
//{
//    int i;
//    for (i = lieNum - 1; i - 8 >= -1; i -= 8)
//    {
//        //每轮处理8个消元子，范围：首项在 i-7 到 i
//        for (int j = 0; j < pasNum; j++)
//        {
//
//            //看4535个被消元行有没有首项在此范围内的
//            while (Pas[j][Num] <= i && Pas[j][Num] >= i - 7)
//            {
//                int index = Pas[j][Num];
//                if (Act[index][Num] == 1)//消元子不为空
//                {
//                    //Pas[j][]和Act[（Pas[j][18]）][]做异或
//                    for (int k = 0; k < Num; k++)
//                    {
//                        Pas[j][k] = Pas[j][k] ^ Act[index][k];
//                    }
//
//                    //更新Pas[j][18]存的首项值
//                    //做完异或之后继续找这个数的首项，存到Pas[j][18]，若还在范围里会继续while循环
//                    //找异或之后Pas[j][ ]的首项
//                    int num = 0, S_num = 0;
//                    for (num = 0; num < Num; num++)
//                    {
//                        if (Pas[j][num] != 0)
//                        {
//                            unsigned int temp = Pas[j][num];
//                            while (temp != 0)
//                            {
//                                temp = temp >> 1;
//                                S_num++;
//                            }
//                            S_num += num * 32;
//                            break;
//                        }
//                    }
//                    Pas[j][Num] = S_num - 1;
//
//                }
//                else//消元子为空
//                {
//                    //Pas[j][]来补齐消元子
//                    for (int k = 0; k < Num; k++)
//                        Act[index][k] = Pas[j][k];
//
//                    Act[index][Num] = 1;//设置消元子非空
//                    break;
//                }
//
//            }
//        }
//    }
//
//    for (i = i + 8; i >= 0; i--)
//    {
//        //每轮处理1个消元子，范围：首项等于i
//
//        for (int j = 0; j < pasNum; j++)
//        {
//            //看53个被消元行有没有首项等于i的
//            while (Pas[j][Num] == i)
//            {
//                if (Act[i][Num] == 1)//消元子不为空
//                {
//                    //Pas[j][]和Act[i][]做异或
//                    for (int k = 0; k < Num; k++)
//                    {
//                        Pas[j][k] = Pas[j][k] ^ Act[i][k];
//                    }
//
//                    //更新Pas[j][18]存的首项值
//                    //做完异或之后继续找这个数的首项，存到Pas[j][18]，若还在范围里会继续while循环
//                    //找异或之后Pas[j][ ]的首项
//                    int num = 0, S_num = 0;
//                    for (num = 0; num < Num; num++)
//                    {
//                        if (Pas[j][num] != 0)
//                        {
//                            unsigned int temp = Pas[j][num];
//                            while (temp != 0)
//                            {
//                                temp = temp >> 1;
//                                S_num++;
//                            }
//                            S_num += num * 32;
//                            break;
//                        }
//                    }
//                    Pas[j][Num] = S_num - 1;
//
//                }
//                else//消元子为空
//                {
//                    //Pas[j][]来补齐消元子
//                    for (int k = 0; k < Num; k++)
//                        Act[i][k] = Pas[j][k];
//
//                    Act[i][Num] = 1;//设置消元子非空
//                    break;
//                }
//            }
//        }
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
//        // Test f_ordinary
//        QueryPerformanceCounter((LARGE_INTEGER*)&head);
//        //f_cache_optimized();
//        f_ordinary();
//        QueryPerformanceCounter((LARGE_INTEGER*)&tail);
//        seconds = (tail - head) * 1000.0 / freq;
//        cout << "Data size: " << dataSize[i] << ",f_ordinary() took " << seconds << " ms\n";
//
//        // Reset arrays for the next test
//        memset(Act, 0, sizeof(unsigned int) * (Num + 1));
//        memset(Pas, 0, sizeof(unsigned int) * (Num + 1));
//    }
//
//    return 0;
//}
