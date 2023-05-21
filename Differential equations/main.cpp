#include "third_party/Eigen/Core"
#include "iostream"
#include <vector>
#include <fstream>
#define _USE_MATH_DEFINES


using namespace std;

template<unsigned dim>
using Vector = Eigen::Vector<double, dim>;

class RightPart{
public:
    constexpr static unsigned dim = 1;
    using Time = double;
    using State = Vector<dim>;
    using RP = Vector<dim>;
    struct StateAndTime {
        Time time;
        State state;
    };

    struct Parameters{};

    static RP calcRightPart(const StateAndTime& stateAndTime, const Parameters& params) noexcept {
        return - stateAndTime.state;
    }
};



template<typename RightPartAbstract>
typename RightPartAbstract::StateAndTime oneStepRK4(const typename RightPartAbstract::StateAndTime& stateAndTime,
                                                    const typename RightPartAbstract::Parameters& params,
                                                    const double h) noexcept {
    //k1
    const typename RightPartAbstract::RP k1 = RightPartAbstract::calcRightPart(stateAndTime, params);

    //k2
    typename RightPartAbstract::StateAndTime change_result;
    change_result.time = stateAndTime.time + h/2;
    change_result.state = stateAndTime.state + h/2 * k1;
    const typename RightPartAbstract::RP k2 = RightPartAbstract::calcRightPart(change_result, params);

    //k3
    change_result.state = stateAndTime.state + h/2 * k2;
    const typename RightPartAbstract::RP k3 = RightPartAbstract::calcRightPart(change_result, params);

    //k4
    change_result.time = stateAndTime.time + h;
    change_result.state = stateAndTime.state + h * k3;
    const typename RightPartAbstract::RP k4 = RightPartAbstract::calcRightPart(change_result, params);

    typename RightPartAbstract::StateAndTime result;
    result.state = stateAndTime.state + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    result.time = stateAndTime.time + h;

    return result;
}

//цикл
template<typename RightPartAbstract>
std::vector<typename RightPartAbstract::StateAndTime> RK4_cycle(int iterations, double h,
                                                                const typename RightPartAbstract::StateAndTime& Y0,
                                                                const typename RightPartAbstract::Parameters& t){

     std::vector<typename RightPartAbstract::StateAndTime> result = {Y0};

    for(int i = 0; i < iterations; i++){
        result.push_back(oneStepRK4<RightPartAbstract>(result.back(), t, h));
    };

    return result;
};

int main(){

    const typename RightPart::Parameters t = {};
    typename RightPart::State b;
    b(0) = 1;
    const typename RightPart::StateAndTime Y_T = {0, b};

    ofstream file1, file2;
    file1.open("/Users/paveldibrovenko/Desktop/deductions/Derivative/Derivative/Empty.txt", std::ios_base::app); // окрываем файл для записи
    file2.open("/Users/paveldibrovenko/Desktop/deductions/Derivative/Derivative/Empty2.txt", std::ios_base::app);

    for (double h = 0.0001; h < 4; h*=2){

        const auto f = RK4_cycle<RightPart>(100, h, Y_T, t);
        double max = 0;

        for ( auto &number : f ){
            double err = abs(number.state(0) - exp(-number.time));
            max = std::max(err, max);
        }
        std::cout << log(max) << std::endl;
        file1 << log( max ) << endl;
        file2 << log(h) << endl;

    };
    file1.close();
    file2.close();


}
