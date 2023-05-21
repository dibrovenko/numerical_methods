#define _USE_MATH_DEFINES
#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <fstream>
#include <iterator>
#include <string>
#include <optional>

using namespace std;

double eq1(double x) {
    return x - 0.1 * sin(x) - M_PI / 4;
}

double eq2(double x) {
    return tan(x) - 4.0 * x / M_PI;
}

double eq3(double x) {
    return log( ( exp(x) + exp(-x) ) / 2.0 );
}
/**
* Функция для решения уравнения func(x) = 0 при помощи метода половинного деления
* a - левая граница отрезка локализации корня
* b - правая гранича отрезка локализации корня
* func - функция, корень которой нужно найти
* numberOfIterations - количество итераций, которое должен освершить метод
* ДОСТАТОЧНЫЕ УСЛОВИЯ СХОДИМОСТИ ПРОВЕРЯТЬ НЕ НУЖНО!
**/
double bisectionMethod(double a, double b,  const function<double(double)> & func, unsigned numberOfIterations) noexcept {
    bool grow = (func(a) <= 0);

    double ans = (a + b) / 2;
    if (!numberOfIterations)
        return ans;

    if ((grow && func(ans) > 0) || (!grow && func(ans) < 0))
        return bisectionMethod(a, ans, func, numberOfIterations - 1);
    else
        if ((grow && func(ans) < 0) || (!grow && func(ans) > 0))
            return bisectionMethod(ans, b, func, numberOfIterations - 1);
        else
            return ans;
}

/**
* Функция для решения уравнения func(x) = 0 при помощи метода простой итерации с релаксацией
* inital - начальное приближение для решения
* func - функция, корень которой нужно найти
* tau - параметр в методе простой итерации
* numberOfIterations - количество итераций, которое должен освершить метод
**/
[[nodiscard]] double simpleIterationMethod(double inital, const function<double(double)>& func, double tau, unsigned numberOfIterations) noexcept {
    double x = inital;
    for (int i = 0; i < numberOfIterations; i++) {
        x = x + tau * func(x);
    };
    return x;
}


double deriv(double x, const function<double(double)>& func) {
    std::array<double, 5> koef = { 0.08333333333 , -0.66666666667 , 0, 0.66666666667, -0.08333333333 };
    std::array<double, 5> set = {-2, -1, 0, 1, 2 };

    double h = 1e-4;
    double ans = 0;
    for (int i = 0; i < set.size(); i++)
        ans += (func(x + h * set[i]) * koef[i]);
        ans /= h;
    
    return ans;
};

/**
* Функция для решения уравнения func(x) = 0 при помощи метода Ньютона
* inital - начальное приближение для решения
* func - функция, корень которой нужно найти
* numberOfIterations - количество итераций, которое должен освершить метод
* Производную стоит оценить численно (предложите схему не ниже 2 порядка самостоятельно)
**/
double newtonMethod(double inital, const function<double(double)>& func, unsigned numberOfIterations) noexcept {
    if (!numberOfIterations)
        return inital;

    return newtonMethod(inital - func(inital) / deriv(inital, func), func, numberOfIterations - 1);
}


int main() {
    
    
    double f;
    double a;
    for (int i = 100; i < 100000; i+=100) {
        f = simpleIterationMethod(1.00, eq2, -0.5, i);
        if ( abs(0.78 - f) > 0.1){
            cout << f << endl;}
    };
    
    f = simpleIterationMethod(1.00, eq2, -0.5, 50);
    cout << f << endl;
    
    double d;
    d = bisectionMethod(0.6, 0.9, eq2, 1000);
    cout << d << endl;
    
    double m;
    m = bisectionMethod(-0.2, 0.2, eq3, 1000);
    cout << m << endl;
    m = simpleIterationMethod(0.2, eq3, -10.0, 1000);
    cout << m << endl;
    
    m = newtonMethod(0.9, eq1, 4);
    cout << m << endl;
    
    return 0;
}
