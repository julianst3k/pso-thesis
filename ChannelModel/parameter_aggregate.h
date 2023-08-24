#pragma once
#include <pybind11/pybind11.h>
#include <vector>
#include <pybind11/numpy.h>
#include <string>
#include <iostream>

namespace py = pybind11;

template<class T>
std::vector<T> make_vector_1d_numpy(py::array_t<T> py_array){
    auto ret = std::vector<T>(py_array.data(), py_array.data()+py_array.size());
    return ret;
};
namespace parameter_aggregate{
struct WallParameters{
    WallParameters(float reflection_coef) : pw(reflection_coef) {}
    float get_reflection_param(){return pw;};
    float pw;
};

struct TransmitterParameters{
    TransmitterParameters(py::array_t<float> coordinate, float beta, float alpha, float ang_rad, float m) :
    coordinate(make_vector_1d_numpy(coordinate)), beta(beta), alpha(alpha), ang_rad(ang_rad), m(m){};
    py::array_t<float> get_coord(){
        float* ret = (float*)malloc(3*sizeof(float));
        size_t  it = 0;
        for(float e: this->coordinate){
            ret[it] = e;
            it++;
        }

        return py::array_t<float>({3}, {sizeof(float)}, ret);
        };
    std::vector<float> coordinate;
    float alpha;
    float beta;
    float ang_rad;
    float m;
};

struct TransmitterAggregate{
    void pushTransmitter(TransmitterParameters *transmitter){transmitters.push_back(transmitter);};
    int get_head(){
        std::cout << (this->transmitters[0])->coordinate[1];
        std::cout << (this->transmitters[0])->coordinate[2];
        return (this->transmitters[0])->coordinate[0];}
    std::vector<TransmitterParameters*> transmitters;
};

struct TransmitterConfigurations{
    void pushTransmitter(TransmitterAggregate *transmitter){transmitters.push_back(transmitter);};
    std::string get_head(){return std::to_string(this->transmitters[0]->transmitters[0]->coordinate[0]);}

    std::vector<TransmitterAggregate*> transmitters;
};

struct ReceiverParameters{
    ReceiverParameters(float Ap, float eta, float fov, float alpha, float ele, py::array_t<float> center, py::array_t<float> coordinate) :
    Ap(Ap), eta(eta), fov(fov), alpha(alpha), ele(ele), center(make_vector_1d_numpy(center)), coordinate(make_vector_1d_numpy(coordinate)) {}
        py::array_t<float> get_coord(){
        float* ret = (float*)malloc(3*sizeof(float));
        size_t  it = 0;
        for(float e: this->coordinate){
            ret[it] = e;
            it++;
        }

        return py::array_t<float>({3}, {sizeof(float)}, ret);
        };
    float r = 0.1;
    float Ap;
    float eta;
    float alpha;
    float ele;
    float fov;
    std::vector<float> center;
    std::vector<float> coordinate;
};

struct ReceiverAggregate{
    void pushReceiver(ReceiverParameters *recv){receivers.push_back(recv);};
    py::array_t<float> get_head(int i){return this->receivers[i]->get_coord();}
    std::vector<ReceiverParameters*> receivers;

};

struct ReceiverConfigurations{
    void pushReceiver(ReceiverAggregate *recv){receivers.push_back(recv);};
    py::array_t<float> get_head(int i, int j){return this->receivers[i]->receivers[j]->get_coord();}
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