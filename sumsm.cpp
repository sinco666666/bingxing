#include <iostream>
#include <windows.h>
using namespace std;

#define ull unsigned long long int

const ull N = 4096 * 128;
ull a[N];
int LOOP = 100;

void init()
{
    for (ull i = 0; i < N; i++)
        a[i] = i;
}

void Trivial()
{
    long long int begin, end, freq;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&begin);
    for (int l = 0; l < LOOP; l++)
    {
        // init();
        ull sum = 0;
        for (int i = 0; i < N; i++)
            sum += a[i];
    }
    QueryPerformanceCounter((LARGE_INTEGER*)&end);
    cout << "Trivial:" << (end - begin) * 1000.0 / freq / LOOP << "ms" << endl;
}

void optimization()
{
    long long int begin, end, freq;
    QueryPerformanceFrequency((LARGE_INTEGER*)&freq);
    QueryPerformanceCounter((LARGE_INTEGER*)&begin);
    for (int l = 0; l < LOOP; l++)
    {
        ull sum1 = 0, sum2 = 0;
        for (int i = 0; i < N - 1; i += 2)
            sum1 += a[i], sum2 += a[i + 1];
        ull sum = sum1 + sum2;
    }
    QueryPerformanceCounter((LARGE_INTEGER*)&end);
    cout << "optimization:" << (end - begin) * 1000.0 / freq / LOOP << "ms" << endl;
}

int main()
{
    init(); cout << "n=" << N << endl;
    //Trivial();
    optimization();
}

