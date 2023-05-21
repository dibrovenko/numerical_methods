#include <iostream>
#include <sstream>
#include <vector>
#include <array>
#include <string>
#include <math.h>
#include <optional>
#include <fstream>
#include <iterator>

using namespace std;


template<unsigned int n>
[[nodiscard]] optional<array<double, n>> firstDer5(array<double, n>& REFpoints) noexcept{
   
    unsigned long count = 0;
    int count_zero = 0;
    int count_L = 0;
    
    // Нам нужно решить систему: A_matrix * x_array = b_array
    array<double, n> b_array;
    b_array.fill(0);
    array<double, n> x_array;
    array<std::array<double, n+1>, n> A_matrix;
    
    // Проверка входных данных
    for(int i = 0; i < n-1; ++i){
        if ( REFpoints[i+1] - REFpoints[i] != 1 ){count++;}
        if ( REFpoints[i] < 0 ){count_L++;}
        if ( REFpoints[i] == 0 ){count_zero++;}
    };
    if ( REFpoints[n-1] == 0 ){count_zero++;};
    
    if (count > 0 || count_zero != 1){
        cout << "пошла херня " << endl;
        return nullopt;}
    else
        //заполняем вектор
        b_array[1] = 1.0;
    //заполняем расширенную матрицу A
    for(int i = 0; i < n; ++i){
        for(int j = 0; j <= n; ++j){
            if(j == n)
                A_matrix[i][j] = b_array[i];
            else
                A_matrix[i][j] = pow(j - count_L, i );
        };
    };
    //приводим матрицу к виду, где диагональные элементы больше чем элементы ниже (на одну строчку)
    for (int i = 0; i < n; i++){
        for (int k = i+1; k < n; k++){
            if ( abs(A_matrix[i][i]) < abs(A_matrix[k][i]) ){
                for (int j = 0; j <= n; j++){
                    double temp = A_matrix[i][j];
                    A_matrix[i][j] = A_matrix[k][j];
                    A_matrix[k][j] = temp;
                }
            }
        }
    };
    //приводим матрицу к гаусу
    for (int i = 0; i < n-1; i++){
        for (int k = i+1; k < n; k++){
            
            const double t = A_matrix[k][i] / A_matrix[i][i];
            for (int j = 0; j <= n; j++) {
                A_matrix[k][j] = A_matrix[k][j] - t * A_matrix[i][j];
            }
        }
    };
    //из гауса записываем ответы в x_vect
    for (int i = n-1; i >= 0; i--){
        x_array[i] = A_matrix[i][n];
        
        for (int j = i+1; j < n; j++)
            if (j != i) {x_array[i] = x_array[i] - A_matrix[i][j] * x_array[j];}
        
        x_array[i] = x_array[i] / A_matrix[i][i];
    };
    
    return x_array;
};



[[nodiscard]] optional<vector<double>> firstDer(const vector<double>& REFpoints) noexcept {
    
    unsigned long n = REFpoints.size();
    unsigned long count = 0;
    int count_zero = 0;
    int count_L = 0;

    // Нам нужно решить систему: A_matrix * x_vect = b_vect
    vector<double> b_vect(n, 0);
    vector<double> x_vect(n);
    vector<vector<double>> A_matrix(n, vector<double>(n + 1));
   
    // Проверка входных данных
    for(int i = 0; i < n-1; ++i){
        if ( REFpoints[i+1] - REFpoints[i] != 1 ){count++;}
        if ( REFpoints[i] < 0 ){count_L++;}
        if ( REFpoints[i] == 0 ){count_zero++;}
        };
    if ( REFpoints[n-1] == 0 ){count_zero++;};
   
    if (count > 0 || count_zero != 1){
        cout << "пошла херня " << endl;
        return nullopt;}
    else
        //заполняем вектор
        b_vect[1] = 1.0;
        //заполняем расширенную матрицу A
        for(int i = 0; i < n; ++i){
            for(int j = 0; j <= n; ++j){
                if(j == n)
                    A_matrix[i][j] = b_vect[i];
                else
                    A_matrix[i][j] = pow(j - count_L, i );
            };
        };
        //приводим матрицу к виду, где диагональные элементы больше чем элементы ниже (на одну строчку)
        for (int i = 0; i < n; i++){
            for (int k = i+1; k < n; k++){
                if ( abs(A_matrix[i][i]) < abs(A_matrix[k][i]) ){
                    for (int j = 0; j <= n; j++){
                        double temp = A_matrix[i][j];
                        A_matrix[i][j] = A_matrix[k][j];
                        A_matrix[k][j] = temp;
                    }
                }
            }
        };
        //приводим матрицу к гаусу
        for (int i = 0; i < n-1; i++){
            for (int k = i+1; k < n; k++){
            
                double t = A_matrix[k][i] / A_matrix[i][i];
                for (int j = 0; j <= n; j++) {
                    A_matrix[k][j] = A_matrix[k][j] - t * A_matrix[i][j];
                }
            }
        };
        //из гауса записываем ответы в x_vect
        for (int i = n-1; i >= 0; i--){
            x_vect[i] = A_matrix[i][n];
            
            for (int j = i+1; j < n; j++)
                if (j != i) {x_vect[i] = x_vect[i] - A_matrix[i][j] * x_vect[j];}
            
            x_vect[i] = x_vect[i] / A_matrix[i][i];
        };
    
        return x_vect;
}

double f(double x) {
    return pow(x, x);
}

int main()
{
    string s;
    vector<double> points;
    getline(cin, s);
    stringstream ss(s);
    string i;
    
    while (getline(ss, i, ' ')) {
        points.push_back(stod(i));
    }

    optional<vector<double>> coefficient = firstDer(points);
    
    if (coefficient.has_value()){
        for (int i = 0; i < coefficient.value().size(); i++){
            cout << coefficient.value()[i] << "   " ;
        };
        cout << "- коэфф vector"<< endl;
    }
    
    //НА ОТЛ:
    const unsigned int size = 9;
    array<double, size> points5 = {-4, -3, -2, -1, 0, 1, 2, 3, 4};
    optional<array<double, size>> coefficient5 = firstDer5<size>(points5);
    
    if ( coefficient5.has_value() ){
        for ( auto &number : coefficient5.value() ){
           cout << number << "   " ;
        };
        cout << "- коэфф array"<< endl;
    };
    
    double h = 0.01;
    double dfdx = 0;
    double pi = 3.14159265358979311600;
    
    ofstream file1, file2;
    file1.open("/Users/paveldibrovenko/Desktop/deductions/Derivative/Derivative/Empty.txt"); // окрываем файл для записи
    file2.open("/Users/paveldibrovenko/Desktop/deductions/Derivative/Derivative/Empty2.txt");
  
    
    if (coefficient.has_value()){
        for (double h=1; h > 1e-10; h/=5){
            for (int i = 0; i < coefficient.value().size(); i++){
                dfdx += f(pi + points[i] * h ) * coefficient.value()[i] / h;
            };
            cout <<  dfdx - f(pi) - f(pi) * log(pi) << "  " << h << endl;
            cout << h << " j " << endl;
            file1 << log(abs(dfdx - f(pi) - f(pi) * log(pi)))  << endl;
            file2 << log(h) << endl;
            dfdx = 0;
        };
    }
    
    
    
    return 0;
}
