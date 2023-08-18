#include <pybind11/pybind11.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <parameter_aggregate.h>

#define pi 3.1415926535
namespace pa = parameter_aggregate;
namespace py = pybind11;


float gain(float eta, float inct, float incr, float fov){
    float fov_radian = pi*fov/180;
    if(incr <= fov && incr>=0){
        return std::pow(eta,2)/std::sin(std::pow(fov_radian,2));
    }
    else{
        return 0;
    }
}


std::vector<float> vector_distance(float tx, float ty, float tz, float rx, float ry, float rz){
    std::vector<float> vect{rx-tx, ry-ty, rz-tz};
    return vect;
}
std::vector<float> norm_vector_trans(float alpha, float beta){
    std::vector<float> vect{std::cos(alpha)*std::sin(beta), std::sin(alpha)*std::sin(beta), -std::cos(beta)};
    return vect;
}
std::vector<float> norm_vector_recv(float alpha, float beta){
    std::vector<float> vect{-std::cos(alpha)*std::sin(beta), -std::sin(alpha)*std::sin(beta), std::cos(beta)};
    return vect;
}
float dot_product(std::vector<float> v1, std::vector<float> v2){
    float sum = 0;
    for(int i=0; i < v1.size(); i++){
        sum += v1[i]*v2[i];
    }
    return sum;
}


float HLos(float tx, float ty, float tz, float rx, float ry, float rz, float Ap, float eta, float alphat, float alphar,
    float betat, float betar, float incr, float incj, int m, float fov, float X, float Y, float t, float es, float c){
    std::vector<float> dist_vector; std::vector<float> norm_vec_t; std::vector<float> norm_vec_r;
    std::vector<float> dist_vector_neg;
    float dist_vector_norm; float p1; float p2; float dm; float g;
    dist_vector = vector_distance(tx, ty, tz, rx, ry, rz);
    dist_vector_norm = std::sqrt(std::pow(tx-rx,2)+std::pow(ty-ry,2)+std::pow(tz-rz,2));
    norm_vec_t = norm_vector_trans(alphat, betat);
    p1 = dot_product(dist_vector, norm_vec_t);
    std::transform(dist_vector.begin(), dist_vector.end(), dist_vector_neg.begin(), [](float i) -> float { return -i;});
    norm_vec_r = norm_vector_recv(alphar, betar);
    p2 = dot_product(dist_vector_neg, norm_vec_r);
    g = gain(eta, incr, incj, fov);
    dm = dist_vector_norm/c;
    return std::abs((m+1)*Ap/(2*pi*std::pow(dist_vector_norm,2))*(std::pow(p1,m)/dist_vector_norm)*(p2/dist_vector_norm)*g);
}

float HNLos(float tx, float ty, float tz, float rx, float ry, float rz, float wx, float wy, float wz, float Aw, float pw, float Ap, float eta, float alphat, float alphar,
    float alphaw, float betat, float betar, float betaw, float incr, float incj, int m, float fov, float X, float Y, float t, float es, float c){
    std::vector<float> dist_vector_tw; std::vector<float> norm_vec_t; std::vector<float> norm_vec_r;
    std::vector<float> dist_vector_tw_neg; std::vector<float> dist_vector_wr; std::vector<float> dist_vector_wr_neg;
    float dist_vector_norm_tw; float dist_vector_norm_wr; float p1; float p2; float p3; float p4; float dm; float g; std::vector<float> norm_vec_w;
    dist_vector_tw = vector_distance(tx, ty, tz, wx, wy, wz);
    std::transform(dist_vector_tw.begin(), dist_vector_tw.end(), dist_vector_tw_neg.begin(), [](float i) -> float { return -i;});
    dist_vector_norm_tw = std::sqrt(std::pow(tx-wx,2)+std::pow(ty-wy,2)+std::pow(tz-wz,2));
    dist_vector_wr = vector_distance(wx, wy, wz, rx, ry, rz);
    std::transform(dist_vector_wr.begin(), dist_vector_wr.end(), dist_vector_wr_neg.begin(), [](float i) -> float { return -i;});
    dist_vector_norm_wr = std::sqrt(std::pow(rx-wx,2)+std::pow(ry-wy,2)+std::pow(rz-wz,2));
    norm_vec_t = norm_vector_trans(alphat, betat);
    norm_vec_r = norm_vector_recv(alphar, betar);
    norm_vec_w = norm_vector_recv(alphaw, betaw);
    p1 = dot_product(dist_vector_tw, norm_vec_t);
    p4 = dot_product(dist_vector_wr_neg, norm_vec_r);
    p2 = dot_product(dist_vector_tw_neg, norm_vec_w);
    p3 = dot_product(dist_vector_wr, norm_vec_w);
    g = gain(eta, incr, incj, fov);
    dm = (dist_vector_norm_wr+dist_vector_norm_tw)/c;
    return std::abs((m+1)*Ap*Aw*pw*p1*p2*p3*p4*g/(std::pow(dist_vector_norm_tw,4)*
    std::pow(dist_vector_norm_wr,4)*pi*pi));
}




PYBIND11_MODULE(example, m){
    py::class_<pa::WallParameters>(m, "WallParameters")
    .def(py::init<float> ());
    py::class_<pa::TransmitterParameters>(m, "TransmitterParameters")
    .def(py::init<py::array_t<float>,float,float,float,float> ());
    py::class_<pa::TransmitterAggregate>(m, "TransmitterAggregate")
    .def("pushTransmitter", &pa::TransmitterAggregate::pushTransmitter)
    .def("getHead", &pa::TransmitterAggregate::get_head);
    py::class_<pa::TransmitterConfigurations>(m, "TransmitterConfigurations")
    .def("pushTransmitter", &pa::TransmitterConfigurations::pushTransmitter)
    .def("getHead", &pa::TransmitterConfigurations::get_head);
    py::class_<pa::ReceiverParameters>(m, "ReceiverParameters")
    .def(py::init<float,float,float,float,float,py::array_t<float>,py::array_t<float>> ());
    py::class_<pa::ReceiverAggregate>(m, "ReceiverAggregate")
    .def("pushReceiver", &pa::ReceiverAggregate::pushReceiver)
    .def("getHead", &pa::ReceiverAggregate::get_head);
    py::class_<pa::ReceiverConfigurations>(m, "ReceiverConfigurations")
    .def("pushReceiver", &pa::ReceiverConfigurations::pushReceiver)
    .def("getHead", &pa::ReceiverConfigurations::get_head);;
    py::class_<pa::TunnelParameters>(m, "TunnelParameters")
    .def(py::init<float,float,float>());
    py::class_<pa::SimulationParameters>(m, "SimulationParameters")
    .def(py::init<float,float,py::array_t<float>,py::array_t<float>> ());

}