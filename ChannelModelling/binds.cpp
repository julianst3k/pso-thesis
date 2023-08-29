#include <pybind11/pybind11.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <parameter_aggregate.h>
#include <channel.h>

namespace pa = parameter_aggregate;
namespace py = pybind11;

PYBIND11_MODULE(model, m){
    py::class_<pa::WallParameters>(m, "WallParameters")
    .def(py::init<float> ())
    .def("get_reflection_param", &pa::WallParameters::get_reflection_param);
    py::class_<pa::TransmitterParameters>(m, "TransmitterParameters")
    .def(py::init<py::array_t<float>,float,float,float,float> ())
    .def("getCoordinate", &pa::TransmitterParameters::get_coord, py::return_value_policy::copy);
    py::class_<pa::TransmitterAggregate>(m, "TransmitterAggregate")
    .def(py::init<> ())
    .def("pushTransmitter", &pa::TransmitterAggregate::pushTransmitter)
    .def("getHead", &pa::TransmitterAggregate::get_head);
    py::class_<pa::TransmitterConfigurations>(m, "TransmitterConfigurations")
    .def(py::init<> ())
    .def("pushTransmitter", &pa::TransmitterConfigurations::pushTransmitter)
    .def("getHead", &pa::TransmitterConfigurations::get_head);
    py::class_<pa::ReceiverParameters>(m, "ReceiverParameters")
    .def(py::init<float,float,float,float,float,py::array_t<float>,py::array_t<float>> ());
    py::class_<pa::ReceiverAggregate>(m, "ReceiverAggregate")
    .def(py::init<> ())
    .def("pushReceiver", &pa::ReceiverAggregate::pushReceiver)
    .def("getHead", &pa::ReceiverAggregate::get_head);
    py::class_<pa::ReceiverConfigurations>(m, "ReceiverConfigurations")
    .def(py::init<> ())
    .def("pushReceiver", &pa::ReceiverConfigurations::pushReceiver)
    .def("getHead", &pa::ReceiverConfigurations::get_head);;
    py::class_<pa::TunnelParameters>(m, "TunnelParameters")
    .def(py::init<float,float,float>());
    py::class_<pa::SimulationParameters>(m, "SimulationParameters")
    .def(py::init<float,float,py::array_t<float>,py::array_t<float>> ());
    m.def("channel_calculation", &channel_calculation);
    m.def("channel_calculation", &HLos);
    m.def("channel_calculation", &HNLos);
    m.def("channel_calculation", &incline);
    m.def("channel_calculation", &find_index_t);
    m.def("channel_calculation", &led_pd_channel);
}