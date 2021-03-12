#include <iostream>
#include <omp.h>
//#include <typeinfo>
#include <regex>
#include <vector>
#include <vector>
int main() {
    omp_set_num_threads(1);
//#pragma omp parallel
//    {
//        int tid = omp_get_thread_num();
//        printf("%d\n", tid);
//    }
//    int a[5][6];
//    std::cout << sizeof(a) << std::endl;
    std::vector<int> v1 = {3, 6, 5, 7, 4, 7, 2};
    int *v3 = v1.data();
    uint64_t c = 5;
    v3 = v3 + c;
    int v2[4];
//    v2 = std::vector<int>(v1.begin() + 2, v1.end() - 1);
    for (int i=0;i<5; i++) {
        std::cout << v3[i] << std::endl;
    }
}