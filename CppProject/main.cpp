//#include <iostream>
//#include <stdio.h>
#include <cstdio>
#include <list>
#include <omp.h>
#include <vector>
using namespace std;

//__global__ void VecAdd(float* A, float* B, float* C)
//{
//    int i = threadIdx.x;
//    C[i] = A[i] + B[i];
//}

void compute(int p){
    printf("Haha%d\n", p);
}

void test_vector(const std::vector<int>& v){
    printf("%d", v[0]);
}

int main() {
//    omp_set_num_threads(6);
//    int maxID = omp_get_max_threads();
    std::list<int> second (4,100);
    auto p = second.begin();

    int a = 1;
    std::vector<int> vec1;
    int arr1[] = {10};
    vec1.push_back(1);
    test_vector(vec1);
    return 0;
}
