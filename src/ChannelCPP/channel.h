#include <vector>
#include <cmath>
#include <algorithm>
#include "parameter_aggregate.h"
#include <pybind11/numpy.h>
#include <random>
#include "shadowing.h"

namespace sh = shadowing;
namespace pa = parameter_aggregate;
namespace py = pybind11;
namespace ch{
struct loss_and_time{
    float loss;
    float time;
};
}
py::array_t<float> channel_calculation(pa::WallParameters* wall, pa::TransmitterAggregate* tconfig, pa::ReceiverConfigurations* rconfig, pa::TunnelParameters* tunnel, pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll);
py::array_t<float> response_calculation(pa::WallParameters* wall, pa::TransmitterAggregate* tconfig, pa::ReceiverConfigurations* rconfig, pa::TunnelParameters* tunnel, pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll);
py::array_t<float> response_calculation_scattering(pa::WallParameters* wall, pa::TransmitterAggregate* tconfig, pa::ReceiverConfigurations* rconfig, pa::TunnelParameters* tunnel, pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll);
py::array_t<float> response_calculation_los(pa::WallParameters* wall, pa::TransmitterAggregate* tconfig, pa::ReceiverConfigurations* rconfig, pa::TunnelParameters* tunnel, pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll);
py::array_t<float> response_calculation_nlos(pa::WallParameters* wall, pa::TransmitterAggregate* tconfig, pa::ReceiverConfigurations* rconfig, pa::TunnelParameters* tunnel, pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll);

int find_index_t(py::array_t<float> arr_t, float x);
void initialize_3d_matrix(float*** out_matrix, int size, int h, int w);

void channel_subroutine(pa::WallParameters* wall, pa::TransmitterAggregate* trans, pa::ReceiverAggregate* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll, std::vector<std::vector<float>>* matrix);

void HNLos_Vector(pa::WallParameters* wall, pa::TransmitterParameters* trans, pa::ReceiverParameters* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, std::vector<float>* final_response, float wall_pos, sh::WH_Probabilities* wh_coll);
void HNLos_Vector_Out(pa::WallParameters* wall, pa::TransmitterParameters* trans, pa::ReceiverParameters* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, std::vector<float>* final_response, float wall_pos, sh::WH_Probabilities* wh_coll);

void HLos_Vector(pa::WallParameters* wall, pa::TransmitterParameters* trans, pa::ReceiverParameters* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, std::vector<float>* final_response, sh::WH_Probabilities* wh_coll);
ch::loss_and_time* HNLos(float tx, float ty, float tz, float rx, float ry, float rz, float wx, float wy, float wz, float Aw, float pw, float Ap, float eta, float alphat, float alphar,
    float alphaw, float betat, float betar, float betaw, float incr, float incj, int m, float fov, float X, float Y, float t, float c, sh::WH_Probabilities* wh_coll);

ch::loss_and_time* HLos(float tx, float ty, float tz, float rx, float ry, float rz, float Ap, float eta, float alphat, float alphar,
    float betat, float betar, float incr, float incj, int m, float fov, float X, float Y, float t, float c, sh::WH_Probabilities* wh_coll);

ch::loss_and_time* HScatter(float tx, float ty, float tz, float rx, float ry, float rz, float Ap, float eta, float alphat, float alphar,
    float betat, float betar, int m, float fov, float X, float Y, float t, float c, sh::WH_Probabilities* wh_coll, pa::ScatterParameters* scatters);

void HScatter_Vector(pa::WallParameters* wall, pa::TransmitterParameters* trans, pa::ReceiverParameters* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, std::vector<float>* final_response, float wall_pos, sh::WH_Probabilities* wh_coll);
void HScatter_Vector_Out(pa::WallParameters* wall, pa::TransmitterParameters* trans, pa::ReceiverParameters* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, std::vector<float>* final_response, float wall_pos, sh::WH_Probabilities* wh_coll);



float dot_product(std::vector<float> v1, std::vector<float> v2);

std::vector<float> norm_vector_recv(float alpha, float beta);

std::vector<float> norm_vector_trans(float alpha, float beta);

std::vector<float> vector_distance(float tx, float ty, float tz, float rx, float ry, float rz);

float gain(float eta, float inct, float incr, float fov);

float incline(float xt, float yt, float zt, float xr, float yr, float zr, float alpha, float beta, bool print);

int find_indexs(std::vector<float> arr, float x);

int find_binary(std::vector<float> arr, int low, int high, float x);

void conv(std::vector<float>* final_response, float** h_vector, std::vector<float>* h_led, int sz);

float led_pd_channel(pa::WallParameters* wall, pa::TransmitterParameters* trans, pa::ReceiverParameters* recv, pa::TunnelParameters* tunnel,
pa::SimulationParameters* simulation, sh::WH_Probabilities* wh_coll);

float gain_scattering(pa::ScatterParameters* scatters, float inc_sr);