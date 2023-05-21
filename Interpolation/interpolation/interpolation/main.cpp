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

template<unsigned int NumPoints>
struct Data {
    array<double, NumPoints> X;
    array<double, NumPoints> Y;
};

template<unsigned int n>
Data<n> CosValueCheb(const double& EndSegment) {
    Data<n> Points_Cheb;
    for (int k = 1; k <= n; k++) {
        Points_Cheb.X[k-1] = (EndSegment) / 2 + (EndSegment) / 2 * cos(((2 * k - 1) * M_PI) / (2 * n));
        Points_Cheb.Y[k-1] = cos(Points_Cheb.X[k-1]);
    };
    return Points_Cheb;
};

template<unsigned int n>
Data<n> CosValue(const double& EndSegment) {
    
    double step = EndSegment / (n - 1);
    Data<n> Points;
    
    for (int i = 0; i < n; i++) {
        Points.X[i] = (i * step);
        Points.Y[i] = cos( i * step );
    };
    return Points;
};



template<typename xType, typename yType, unsigned int size>
class NewtonianInterpolator {
    
private:
    array<xType, size> coeffs;
    array<xType, size> x;
    
public:
    /*** Конструктор, который принимает поля класса и копирует их в поля.
    Для создания по массиву x и y использовать фабрику ***/
    explicit NewtonianInterpolator(const array<xType, size>& coef, const array<xType, size>& X) noexcept{
        this->coeffs = coef;
        this-> x = X;
    };

    /*** Это фабрика - она создает объект класса, а если не получается создать, то возвращает std::nullopt ***/
    [[nodiscard]] static optional<NewtonianInterpolator> Build_Interpolar(const array<xType, size>& x, const array<yType, size>& y) noexcept{
        array<yType, size> coef;
        int m = 1;
        
        for (int i = 0; i < size; ++i) {
            coef[i] = y[i];
        };
        
        for (int i = size - 1; i >= 0; i--) {
            for (int j = 1; j <= i; j++) {
                coef[size - j] = (coef[size - j] - coef[size - 1 - j]) / (x[size - j] - x[size - j - m]);
            }
            m++;
        };

        return NewtonianInterpolator(coef, x);
        
    };
    
    /*** Метод, выполняющий подсчет интерполяционного полинома в точке  [[nodiscard]] yType */
    [[nodiscard]] yType interpolate(xType x0) const noexcept{
        yType sum = coeffs[0];
        yType mult = 1;
        for (int i = 0; i < size - 1; ++i) {
            mult *= x0 - x[i];
            sum += mult * coeffs[i+1];
        };
        
        return sum;
    };
};
//template<unsigned int ... numbers> // template parameter pack

template <const unsigned int... n>
void test()
{
    ofstream file1, file2;
    file1.open("/Users/paveldibrovenko/Desktop/deductions/Derivative/Derivative/Empty.txt"); // окрываем файл для записи
    file2.open("/Users/paveldibrovenko/Desktop/deductions/Derivative/Derivative/Empty2.txt");
    
    double f;
    double max = 0;
    array<double, 7> EndSegments = {8, 4, 2, 1, 0.5, 0.25, 0.125};
    
    ([&]{
        for ( auto &number : EndSegments ){
            Data<n> CosPoints = CosValue<n>(number);
            optional<NewtonianInterpolator<double, double, n>> I = NewtonianInterpolator<double, double, n> :: Build_Interpolar(CosPoints.X, CosPoints.Y);
            for(double i = 0; i < number; i += number/1000){
                f = I->interpolate(i);
                if ( abs( f - cos(i) ) > max){ max = abs( f - cos(i) ); };
            };
            file1 << log( max ) << endl;
            file2 << log(number) << endl;
            max = 0;
        };
    }(), ...);
};

int main() {
    
    test<3, 5, 10>();
    
    //Это для чебышева
    ofstream file1, file2;
    file1.open("/Users/paveldibrovenko/Desktop/deductions/Derivative/Derivative/Empty.txt", std::ios_base::app); // окрываем файл для записи
    file2.open("/Users/paveldibrovenko/Desktop/deductions/Derivative/Derivative/Empty2.txt", std::ios_base::app);
    
    double f;
    double max = 0;
    array<double, 7> EndSegments = {8, 4, 2 , 1, 0.5 ,0.25 ,0.125};
    
    const unsigned int t = 5;
    for ( auto &number : EndSegments ){
        Data<t> CosPoints = CosValueCheb<t>(number);
        optional<NewtonianInterpolator<double, double, t>> I = NewtonianInterpolator<double, double, t> :: Build_Interpolar(CosPoints.X, CosPoints.Y);
        for(double i = 0; i < number; i += number/1000){
            f = I->interpolate(i);
            if ( abs( f - cos(i) ) > max){ max = abs( f - cos(i) ); };
        };
        file1 << log( max ) << endl;
        file2 << log(number) << endl;
        max = 0;
    };
    
    return 0;
}


