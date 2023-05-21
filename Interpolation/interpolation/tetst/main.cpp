#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <optional>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace std;


template <unsigned int size>
optional<unsigned int> BinarySearch(const array<double, size>& arr, double X0) noexcept
{
    unsigned int left = 0;
    unsigned int right = size - 1;
    unsigned int middle;

    while (left <= right)
    {
        middle = ( left + right ) / 2;
        if ( X0 <= arr[middle] && X0 >= arr[middle-1] )
            return middle - 1;
        if (arr[middle] < X0)
            left = middle + 1;
        else
            right = middle -1;
    }
    return nullopt;
}
int main()
{
    const unsigned int size = 4;
    array<double, size> m = {0, 1, 2, 3.1};
    //optional<array<double, size>> coefficient5 = firstDer5<size>(points5);
    optional<unsigned int> t = BinarySearch<size>(m, 2.2);
    if (t.has_value()) cout << t.value() ;

}
