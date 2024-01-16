#include <pybind11/pybind11.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include "parameter_aggregate.h"
#include "channel.h"
#include "shadowing.h"

namespace sh = shadowing;
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
    .def(py::init<float,float,py::array_t<float>,py::array_t<float>,pa::ScatterParameters*> ());
    py::class_<sh::WH_Probabilities>(m, "WH_Probabilities")
    .def(py::init<> ())
    .def("push_probability", &sh::WH_Probabilities::push_probability);
    py::class_<sh::Shadowing_Parameters_Coll>(m, "Shadow_Coll")
    .def(py::init<> ())
    .def("create_collection", &sh::Shadowing_Parameters_Coll::create_collection);
    py::class_<pa::ScatterParameters>(m, "ScatteringParameters")
    .def(py::init<int, float, float , float , float , float , float , float>());

    m.def("channel_calculation", &channel_calculation);
    m.def("response_calculation", &response_calculation);
    m.def("response_calculation_scattering", &response_calculation_scattering);
    m.def("response_calculation_los", &response_calculation_los);
    m.def("response_calculation_nlos", &response_calculation_nlos);

    // m.def("channel_calculation", &HLos);
    // m.def("channel_calculation", &HNLos);
    // m.def("channel_calculation", &incline);
    // m.def("channel_calculation", &find_index_t);
    // m.def("channel_calculation", &led_pd_channel);
    m.def("calculate_expectancy", &sh::calculate_expectancy);
}