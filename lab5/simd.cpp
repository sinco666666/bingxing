//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <chrono>
//#include <emmintrin.h>  // Include SSE2 intrinsics
//#include <immintrin.h>  // Include AVX intrinsics
//#include <malloc.h>     // For _aligned_malloc and _aligned_free
//
//using namespace std;
//
//const int col = 23145, elinenum = 14325;
//int bytenum = ((col - 1) / 32 + 1);
//
//class bitmatrix {
//public:
//    int mycol;
//    int* mybyte;
//
//    bitmatrix() : mycol(-1), mybyte(reinterpret_cast<int*>(_aligned_malloc(bytenum * sizeof(int), 32))) {
//        for (int i = 0; i < bytenum; i++) {
//            mybyte[i] = 0;
//        }
//    }
//
//    bitmatrix(const bitmatrix& other) : mycol(other.mycol), mybyte(reinterpret_cast<int*>(_aligned_malloc(bytenum * sizeof(int), 32))) {
//        for (int i = 0; i < bytenum; i++) {
//            mybyte[i] = other.mybyte[i];
//        }
//    }
//
//    bitmatrix& operator=(const bitmatrix& other) {
//        if (this != &other) {
//            mycol = other.mycol;
//            for (int i = 0; i < bytenum; i++) {
//                mybyte[i] = other.mybyte[i];
//            }
//        }
//        return *this;
//    }
//
//    ~bitmatrix() {
//        _aligned_free(mybyte);
//    }
//
//    bool isnull() const {
//        return mycol == -1;
//    }
//
//    void insert(int x) {
//        if (mycol == -1) mycol = x;
//        int a = x / 32, b = x % 32;
//        mybyte[a] |= (1 << b);
//    }
//
//    void doxor(const bitmatrix& b) {
//        for (int i = 0; i < bytenum; i++) {
//            mybyte[i] ^= b.mybyte[i];
//        }
//        updateMyCol();
//    }
//
//    void simd_doxor(const bitmatrix& b) {
//        int i = 0;
//        for (; i <= bytenum - 4; i += 4) {
//            __m128i mb = _mm_loadu_si128(reinterpret_cast<__m128i*>(mybyte + i));
//            __m128i bb = _mm_loadu_si128(reinterpret_cast<__m128i*>(b.mybyte + i));
//            __m128i result = _mm_xor_si128(mb, bb);
//            _mm_storeu_si128(reinterpret_cast<__m128i*>(mybyte + i), result);
//        }
//        for (; i < bytenum; i++) {
//            mybyte[i] ^= b.mybyte[i];
//        }
//        updateMyCol();
//    }
//
    //void avx_doxor(const bitmatrix& b) {
    //    int i = 0;
    //    for (; i <= bytenum - 8; i += 8) {
    //        __m256i mb = _mm256_loadu_si256(reinterpret_cast<__m256i*>(mybyte + i));
    //        __m256i bb = _mm256_loadu_si256(reinterpret_cast<__m256i*>(b.mybyte + i));
    //        __m256i result = _mm256_xor_si256(mb, bb);
    //        _mm256_storeu_si256(reinterpret_cast<__m256i*>(mybyte + i), result);
    //    }
    //    for (; i < bytenum; i++) {
    //        mybyte[i] ^= b.mybyte[i];
    //    }
    //    updateMyCol();
    //}
//
//    void simd_aligned_doxor(const bitmatrix& b) {
//        int i = 0;
//        for (; i <= bytenum - 4; i += 4) {
//            __m128i mb = _mm_load_si128(reinterpret_cast<__m128i*>(mybyte + i));
//            __m128i bb = _mm_load_si128(reinterpret_cast<__m128i*>(b.mybyte + i));
//            __m128i result = _mm_xor_si128(mb, bb);
//            _mm_store_si128(reinterpret_cast<__m128i*>(mybyte + i), result);
//        }
//        for (; i < bytenum; i++) {
//            mybyte[i] ^= b.mybyte[i];
//        }
//        updateMyCol();
//    }
//
//    void avx_aligned_doxor(const bitmatrix& b) {
//        int i = 0;
//        for (; i <= bytenum - 8; i += 8) {
//            __m256i mb = _mm256_load_si256(reinterpret_cast<__m256i*>(mybyte + i));
//            __m256i bb = _mm256_load_si256(reinterpret_cast<__m256i*>(b.mybyte + i));
//            __m256i result = _mm256_xor_si256(mb, bb);
//            _mm256_store_si256(reinterpret_cast<__m256i*>(mybyte + i), result);
//        }
//        for (; i < bytenum; i++) {
//            mybyte[i] ^= b.mybyte[i];
//        }
//        updateMyCol();
//    }
//
//private:
//    void updateMyCol() {
//        mycol = -1;
//        for (int i = bytenum - 1; i >= 0; i--) {
//            for (int j = 31; j >= 0; j--) {
//                if ((mybyte[i] & (1 << j)) != 0) {
//                    mycol = i * 32 + j;
//                    return;
//                }
//            }
//        }
//    }
//};
//
//bitmatrix* eliminator = nullptr;
//bitmatrix* eline = nullptr;
//
//void readdata(const string& eliminatorFile, const string& elineFile) {
//    delete[] eliminator;
//    delete[] eline;
//    eliminator = new bitmatrix[col];
//    eline = new bitmatrix[elinenum];
//
//    ifstream ifs(eliminatorFile);
//    string temp;
//    while (getline(ifs, temp)) {
//        istringstream ss(temp);
//        int x, trow = 0;
//        while (ss >> x) {
//            if (!trow) trow = x;
//            eliminator[trow].insert(x);
//        }
//    }
//    ifs.close();
//
//    ifstream ifs2(elineFile);
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
//void dowork() {
//    for (int i = 0; i < elinenum; i++) {
//        while (!eline[i].isnull()) {
//            int tcol = eline[i].mycol;
//            if (!eliminator[tcol].isnull()) {
//                eline[i].doxor(eliminator[tcol]);
//            }
//            else {
//                eliminator[tcol] = eline[i];
//                break;
//            }
//        }
//    }
//}
//
//void simd_dowork() {
//    for (int i = 0; i < elinenum; i++) {
//        while (!eline[i].isnull()) {
//            int tcol = eline[i].mycol;
//            if (!eliminator[tcol].isnull()) {
//                eline[i].simd_doxor(eliminator[tcol]);
//            }
//            else {
//                eliminator[tcol] = eline[i];
//                break;
//            }
//        }
//    }
//}
//
//void avx_dowork() {
//    for (int i = 0; i < elinenum; i++) {
//        while (!eline[i].isnull()) {
//            int tcol = eline[i].mycol;
//            if (!eliminator[tcol].isnull()) {
//                eline[i].avx_doxor(eliminator[tcol]);
//            }
//            else {
//                eliminator[tcol] = eline[i];
//                break;
//            }
//        }
//    }
//}
//
//void simd_aligned_dowork() {
//    for (int i = 0; i < elinenum; i++) {
//        while (!eline[i].isnull()) {
//            int tcol = eline[i].mycol;
//            if (!eliminator[tcol].isnull()) {
//                eline[i].simd_aligned_doxor(eliminator[tcol]);
//            }
//            else {
//                eliminator[tcol] = eline[i];
//                break;
//            }
//        }
//    }
//}
//
//void avx_aligned_dowork() {
//    for (int i = 0; i < elinenum; i++) {
//        while (!eline[i].isnull()) {
//            int tcol = eline[i].mycol;
//            if (!eliminator[tcol].isnull()) {
//                eline[i].avx_aligned_doxor(eliminator[tcol]);
//            }
//            else {
//                eliminator[tcol] = eline[i];
//                break;
//            }
//        }
//    }
//}
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
//
//int main() {
//    const string eliminatorFile = "x23045.txt";
//    const string elineFile = "b23045.txt";
//
//    readdata(eliminatorFile, elineFile);
//
//    auto start = chrono::high_resolution_clock::now();
//    dowork();
//    auto end = chrono::high_resolution_clock::now();
//    chrono::duration<double> elapsed = end - start;
//    cout << "Normal dowork time: " << elapsed.count() << " seconds." << endl;
//
//    readdata(eliminatorFile, elineFile);
//
//    start = chrono::high_resolution_clock::now();
//    simd_dowork();
//    end = chrono::high_resolution_clock::now();
//    elapsed = end - start;
//    cout << "SIMD dowork time: " << elapsed.count() << " seconds." << endl;
//
//    readdata(eliminatorFile, elineFile);
//
//    start = chrono::high_resolution_clock::now();
//    avx_dowork();
//    end = chrono::high_resolution_clock::now();
//    elapsed = end - start;
//    cout << "AVX dowork time: " << elapsed.count() << " seconds." << endl;
//
//    readdata(eliminatorFile, elineFile);
//
//    start = chrono::high_resolution_clock::now();
//    simd_aligned_dowork();
//    end = chrono::high_resolution_clock::now();
//    elapsed = end - start;
//    cout << "SIMD aligned dowork time: " << elapsed.count() << " seconds." << endl;
//
//    readdata(eliminatorFile, elineFile);
//
//    start = chrono::high_resolution_clock::now();
//    avx_aligned_dowork();
//    end = chrono::high_resolution_clock::now();
//    elapsed = end - start;
//    cout << "AVX aligned dowork time: " << elapsed.count() << " seconds." << endl;
//    //printres();
//    return 0;
//}
//
