#include <pybind11/pybind11.h>
#include <vector>
#include <pybind11/numpy.h>
#include <string>

namespace py = pybind11;

template<class T>
std::vector<T> make_vector_1d_numpy(py::array_t<T> py_array){
    return std::vector<T>(py_array.data(), py_array.data()+py_array.size());
};
namespace parameter_aggregate{
struct WallParameters{
    WallParameters(float reflection_coef) : pw(reflection_coef) {}
    float pw;
};

struct TransmitterParameters{
    TransmitterParameters(py::array_t<float> coordinate, float beta, float alpha, float ang_rad, float m) :
    coordinate(make_vector_1d_numpy(coordinate)), beta(beta), alpha(alpha), ang_rad(ang_rad), m(m){};

    std::vector<float> coordinate;
    float alpha;
    float beta;
    float ang_rad;
    float m;
};

struct TransmitterAggregate{
    void pushTransmitter(TransmitterParameters *transmitter){transmitters.push_back(transmitter);};
    std::string get_head(){return std::to_string(transmitters[0]->coordinate[0]);}
    std::vector<TransmitterParameters*> transmitters;
};

struct TransmitterConfigurations{
    void pushTransmitter(TransmitterAggregate *transmitter){transmitters.push_back(transmitter);};
    std::string get_head(){return std::to_string(transmitters[0]->transmitters[0]->coordinate[0]);}
    std::vector<TransmitterAggregate*> transmitters;
};

struct ReceiverParameters{
    ReceiverParameters(float Ap, float eta, float fov, float angle, float ele, py::array_t<float> center, py::array_t<float> coordinate) :
    Ap(Ap), eta(eta), fov(fov), angle(angle), ele(ele), center(make_vector_1d_numpy(center)), coordinate(make_vector_1d_numpy(coordinate)) {}
    float r = 0.1;
    float Ap;
    float eta;
    float angle;
    float ele;
    float fov;
    std::vector<float> center;
    std::vector<float> coordinate;
};

struct ReceiverAggregate{
    void pushReceiver(ReceiverParameters *recv){receivers.push_back(recv);};
    std::string get_head(){return std::to_string(receivers[0]->coordinate[0]);}
    std::vector<ReceiverParameters*> receivers;

};

struct ReceiverConfigurations{
    void pushReceiver(ReceiverAggregate *recv){receivers.push_back(recv);};
    std::string get_head(){return std::to_string(receivers[0]->receivers[0]->coordinate[0]);}
    std::vector<ReceiverAggregate*> receivers;
};

struct TunnelParameters{
    TunnelParameters(float x, float y, float z) : x(x), y(y), z(z){
    }
    float x;
    float y;
    float z;

};

struct SimulationParameters{
    SimulationParameters(float t, float c, py::array_t<float> time, py::array_t<float> h_led) :
    t(t), c(c), time(make_vector_1d_numpy(time)), h_led(make_vector_1d_numpy(h_led)) {}
    float t;
    float c;
    std::vector<float> time;
    std::vector<float> h_led;
};
};