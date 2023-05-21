#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <optional>
#include <functional>

using namespace std;

//Прозводная в точке многочлена лежандра
double Pn(double x, int n)
{
    if(n == 0) return 1;
    if(n == 1) return x;
    return ( (2* n - 1) * x * Pn(x, n-1) - (n - 1) * Pn(x, n - 2) ) / n;
}

double Pn_derivative(double x, int n)
{
    return n / (1 - x*x) * ( Pn(x, n-1) - x * Pn(x, n));
}

/*vector<vector<double>> xx = { {0}, {0}, {-0.774597, 0 , 0.774597}, {-0.861136, -0.339981, 0.339981, 0.861136 }, {-0.90618, -0.538469, 0, 0.538469, 0.90618} }; */

vector<vector<double>> xx = { {0}, {0}, { -sqrt( 3.0 / 5.0 ), 0 , sqrt( 3.0 / 5.0 ) },
    {-0.861136, -0.339981, 0.339981, 0.861136 },
    {double((-1)) / 3 * sqrt(5 + 2 * sqrt(double(10) / 7)),
    double((-1)) / 3 * sqrt(5 - 2 * sqrt(double(10) / 7)), 0,
    double(1) / 3 * sqrt(5 - 2 * sqrt(double(10) / 7)),
    double(1) / 3 * sqrt(5 + 2 * sqrt(double(10) / 7))} };

template<unsigned int n>
struct Data {
    array<double, n> X;
    array<double, n> W;
};

template<unsigned int n>
Data<n> Make_Xi_and_Wi() {
    Data<n> Xi_and_Wi;
    for(double i = 0; i < n; i++){
        Xi_and_Wi.X[i] = xx[n-1][i];
        Xi_and_Wi.W[i] = 2.0 / ( 1 - Xi_and_Wi.X[i] * Xi_and_Wi.X[i] ) / pow( Pn_derivative(Xi_and_Wi.X[i], n), 2 );
    };
    
    return Xi_and_Wi;
};


double func(double x)
{
    return cos(x);
};
/*** Функция для интегрирования на одном отрезке
* a - начало отезка
* b - конец отрезка
* n - количество точек для квадратуры
* func - функция для интегрирования
***/
template<unsigned n>
[[nodiscard]] double integrateOneSeg(double a, double b, const function<double(double)> & func) noexcept{
    double S = 0;
    Data<n> Xi_and_Wi = Make_Xi_and_Wi<n>();
    
    for (int i = 0; i < n; ++i) {
        S += (b - a) / 2.0 * (func(((b + a) / 2.0) + ((b - a) / 2.0) * Xi_and_Wi.X[i])) * Xi_and_Wi.W[i];
    };
    return S;
    
};

/*** Функция для интегрирования на всем отрезке
* a - начало отезка
* b - конец отрезка
* n - количество точек для квадратуры
* s - количество отрезков
* func - функция для интегрирования
***/
template<unsigned n>
[[nodiscard]] double integrate(double a, double b, unsigned s, const function<double(double)> & func) noexcept{
    double S = 0;
    double step = (b - a) / s;
    for (unsigned i = 0; i < s; ++i)
    {
        S += integrateOneSeg<n>(a + step * i, a + step * (i+1), func);
    };
    return S;
};


template <const unsigned int... n>
void test()
{
    ofstream file1, file2;  // окрываем файлы для записи
    file1.open("/Users/paveldibrovenko/Desktop/deductions/Integration/Integration/Integration/Integration_1.txt");
    file2.open("/Users/paveldibrovenko/Desktop/deductions/Integration/Integration/Integration/Integration_2.txt");
    double a = 0;
    double b = 10;
    double f;
    ([&]{
        for(unsigned s = 1; s < 11; s++){
            f = integrate<n>(a, b, s, func);
            file1 << log( abs( f + sin(a) - sin(b) ) ) << endl;
            file2 << log(s) << endl;
        };
    }(), ...);
};


int main() {
    
    test<3, 4, 5>();
    return 0;
}

