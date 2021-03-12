#include <iostream>
#include <stdio.h>
//#include <typeinfo>
#include <regex>
#include <vector>
//#include <string>
//#include <cstdio>
//using namespace std;
/*
void sum_up(int& a, int b){
    a = a + b;
}

void operation_printer(std::string& s1, std::string&& s2, const char& operation){
    double result;
    double a1 = std::stod(s1), a2 = std::stod(s2);

    switch(operation){
        case '+':
            result = a1 + a2;
            break;
        case '-':
            result = a1 - a2;
            break;
        case '*':
            result = a1 * a2;
            break;
        case '/':
            result = a1 / a2;
            break;
        default:
            result = 0;
    }
    if (std::regex_search(s1, std::regex(".")) || std::regex_search(s2, std::regex("."))){
        std::cout << "Float output: " << s1 << " " << operation << " " << s2 << " = " << result << "\n";
    }
    else{
        std::cout << s1 << " " << operation << " " << s2 << " = " << (int)result << "\n";
    }
}

*/
void update(double &&a, double &b) {
    b += a;
    std::cout << b << a << std::endl;
}

class Cool{
private:
    int ax = 0;
public:
    explicit Cool(int a){
        ax = a;
    }

    int get_ax(){
        return ax;
    }
    int b = ax;
};

int main() {
    double a = 1;
    std::vector<int> b = {2,5,1,7,2,5};
    int c = b.size();
    Cool cool = Cool(5);
    std::cout << cool.get_ax() << std::endl;
}

auto begin = std::chrono::steady_clock::now();

auto end = std::chrono::steady_clock::now();
std::cout << "elapsed time: " << (end - begin).count() / 1000000000.0 << " sec"<< std::endl;
std::cout << "Performance: " << (end - begin).count() /
(num_of_iterations * (1 + 8 + num_of_asteroids)) << " nanosecs / particle update"<< std::endl;

