#include <pmmintrin.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <tmmintrin.h>
#include <nmmintrin.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <immintrin.h>
#include <windows.h>
#include <cstring>
#include <mpi.h> // ���MPIͷ�ļ�

using namespace std;

// ���ݹ�ģ�����齫�� main �ж�̬����
unsigned int Act[37960][1188]; // ����ģ������Ԥ���С
unsigned int Pas[37960][1188];

int Num;
int pasNum;
int lieNum;

void init_A(const string& filename)
{
    ifstream infile(filename);
    char fin[10000] = { 0 };
    int index;

    while (infile.getline(fin, sizeof(fin)))
    {
        stringstream line(fin);
        int biaoji = 0;
        unsigned int a;

        while (line >> a)
        {
            if (biaoji == 0)
            {
                index = a;
                biaoji = 1;
            }
            int k = a % 32;
            int j = a / 32;

            int temp = 1 << k;
            Act[index][Num - 1 - j] |= temp; // Use bitwise OR to set the bit
        }
        Act[index][Num] = 1; // Mark this row as non-empty
    }
}

void init_P(const string& filename)
{
    ifstream infile(filename);
    char fin[10000] = { 0 };
    int index = 0;

    while (infile.getline(fin, sizeof(fin)))
    {
        stringstream line(fin);
        int biaoji = 0;
        unsigned int a;

        while (line >> a)
        {
            if (biaoji == 0)
            {
                Pas[index][Num] = a;
                biaoji = 1;
            }
            int k = a % 32;
            int j = a / 32;

            int temp = 1 << k;
            Pas[index][Num - 1 - j] |= temp;
        }
        index++;
    }
}

void f_serial() {
    int i;
    for (i = lieNum - 1; i >= 0; i -= 8) {
        for (int j = 0; j < pasNum; j++) {
            int idxPas = Pas[j][Num];
            if (idxPas <= i && idxPas >= i - 7) {
                int index = idxPas;
                if (Act[index][Num] == 1) {
                    for (int k = 0; k < Num; k++) {
                        Pas[j][k] ^= Act[index][k];
                    }
                    int S_num = 0;
                    for (int num = 0; num < Num; num++) {
                        if (Pas[j][num] != 0) {
                            for (int bit = 0; bit < 32; bit++) {
                                if ((Pas[j][num] & (1 << bit)) != 0) {
                                    S_num = bit + num * 32;
                                    break;
                                }
                            }
                            if (S_num != 0) break;
                        }
                    }
                    Pas[j][Num] = S_num;
                }
                else {
                    memcpy(Act[index], Pas[j], Num * sizeof(Pas[0][0]));
                    Act[index][Num] = 1;
                }
            }
        }
    }
}

void f_parallel(int rank, int size) {
    int i;
    for (i = lieNum - 1 - rank; i >= 0; i -= size) {
        for (int j = 0; j < pasNum; j++) {
            int idxPas = Pas[j][Num];
            if (idxPas <= i && idxPas >= i - 7) {
                int index = idxPas;
                if (Act[index][Num] == 1) {
                    for (int k = 0; k < Num; k++) {
                        Pas[j][k] ^= Act[index][k];
                    }
                    int S_num = 0;
                    for (int num = 0; num < Num; num++) {
                        if (Pas[j][num] != 0) {
                            for (int bit = 0; bit < 32; bit++) {
                                if ((Pas[j][num] & (1 << bit)) != 0) {
                                    S_num = bit + num * 32;
                                    break;
                                }
                            }
                            if (S_num != 0) break;
                        }
                    }
                    Pas[j][Num] = S_num;
                }
                else {
                    memcpy(Act[index], Pas[j], Num * sizeof(Pas[0][0]));
                    Act[index][Num] = 1;
                }
            }
        }
    }
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv); // ��ʼ��MPI����
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // ��ȡ��ǰ���̵���
    MPI_Comm_size(MPI_COMM_WORLD, &size); // ��ȡ��������

    // Data sizes, corresponding pasNum and Num values
    const int dataSize[] = { 130, 254, 562, 1011, 2362, 3799, 8399, 23045, 37960 };
    const int pasNums[] = { 8, 53, 53, 263, 453, 1953, 4535, 14325, 14921 };
    const int Nums[] = { 5, 8, 18, 32, 74, 119, 263, 721, 1187 };
    const int dataCount = sizeof(dataSize) / sizeof(dataSize[0]);

    // Loop over each data size for testing
    for (int i = 0; i < dataCount; ++i)
    {
        if (rank == 0) { // ֻ�������̽��г�ʼ��
            lieNum = dataSize[i];
            pasNum = pasNums[i];
            Num = Nums[i];

            string xFilename = "x" + to_string(dataSize[i]) + ".txt";
            string bFilename = "b" + to_string(dataSize[i]) + ".txt";

            init_A(xFilename);
            init_P(bFilename);
        }

        // �㲥���������н���
        MPI_Bcast(&lieNum, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&pasNum, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&Num, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(Act, 37960 * 1188, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI_Bcast(Pas, 37960 * 1188, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        // Performance measurement variables
        double seconds;
        long long head, tail, freq;

        // Get frequency
        QueryPerformanceFrequency((LARGE_INTEGER*)&freq);

        // Test f_serial
        if (rank == 0) {
            QueryPerformanceCounter((LARGE_INTEGER*)&head);
            f_serial();
            QueryPerformanceCounter((LARGE_INTEGER*)&tail);
            seconds = (tail - head) * 1000.0 / freq;
            cout << "Data size: " << dataSize[i] << ", f_serial() took " << seconds << " ms\n";
        }

        // Reset arrays for the next test
        memset(Act, 0, sizeof(unsigned int) * (Num + 1));
        memset(Pas, 0, sizeof(unsigned int) * (Num + 1));

        if (rank == 0) { // ֻ�����������³�ʼ������
            init_A("x" + to_string(dataSize[i]) + ".txt");
            init_P("b" + to_string(dataSize[i]) + ".txt");
        }

        // �㲥�µĳ�ʼ�����ݸ����н���
        MPI_Bcast(Act, 37960 * 1188, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI_Bcast(Pas, 37960 * 1188, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        // Test f_parallel
        if (rank == 0) {
            QueryPerformanceCounter((LARGE_INTEGER*)&head);
        }

        f_parallel(rank, size);

        // ���������ռ�����Ļ�����
        unsigned int* recvbuf = nullptr;
        if (rank == 0) {
            recvbuf = new unsigned int[pasNum * (Num + 1) * size];
        }

        // Gather results to the root process
        MPI_Gather(Pas, pasNum * (Num + 1), MPI_UNSIGNED, recvbuf, pasNum * (Num + 1), MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            QueryPerformanceCounter((LARGE_INTEGER*)&tail);
            seconds = (tail - head) * 1000.0 / freq;
                cout << "Data size: " << dataSize[i] << ", f_parallel() took " << seconds << " ms\n";
        }

        // �ͷ������ռ�����Ļ�����
        if (rank == 0) {
            delete[] recvbuf;
        }

        // ���³�ʼ�������Ա��´β���
        memset(Act, 0, sizeof(unsigned int) * (Num + 1));
        memset(Pas, 0, sizeof(unsigned int) * (Num + 1));
    }

    MPI_Finalize(); // ����MPI����
    return 0;
}
