//#include<iostream>
//#include<fstream>
//#include<sstream>
//#include <ctime>
//#include <ratio>
//#include <chrono>
//using namespace std;
//const int col = 3799, elinenum = 1953; //����������Ԫ����
//int bytenum = (col - 1) / 32 + 1;   //ÿ��ʵ���е�byte��������
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
//    void doxor(bitmatrix b) {  //�����������������ڽ�����ڱ�ʵ���У�ֻ�б���Ԫ����ִ����һ����,����������Ҫ��������
//        for (int i = 0; i < bytenum; i++)
//            mybyte[i] ^= b.mybyte[i];
//        for (int i = bytenum - 1; i >= 0; i--)
//            for (int j = 31; j >= 0; j--)
//                if ((mybyte[i] & (1 << j)) != 0) {
//                    mycol = i * 32 + j;
//                    return;
//                }
//        mycol = -1;
//    }
//};
//bitmatrix* eliminer = new bitmatrix[col], * eline = new bitmatrix[elinenum];
//void readdata() {
//    ifstream ifs;
//    ifs.open("x3799.txt");  //��Ԫ��
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
//    ifs2.open("b3799.txt");     //����Ԫ��,���뷽ʽ����Ԫ�Ӳ�ͬ
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
//void dowork() {  //��Ԫ
//    for (int i = 0; i < elinenum; i++) {
//        while (!eline[i].isnull()) {  //ֻҪ����Ԫ�зǿգ�ѭ������
//            int tcol = eline[i].mycol;  //����Ԫ�е�����
//            if (!eliminer[tcol].isnull())    //������ڶ�Ӧ��Ԫ��
//                eline[i].doxor(eliminer[tcol]);
//            else {
//                eliminer[tcol] = eline[i];    //���ڱ���Ԫ������Ϊ��Ԫ�Ӻ󲻲��������������ֱ����=��ǳ����
//                break;
//            }
//        }
//    }
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
//    //printres();
//    system("pause");
//    return 0;
//}
