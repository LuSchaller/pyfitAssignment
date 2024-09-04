#include <pybind11/pybind11.h>
namespace py = pybind11;

#include "fitAssignment.cc"

PYBIND11_MODULE(pyfitAssignment, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("setbestcombi", &setbestcombi, "A function that sets the best combination of jets based on the chi2 value of the kinematic fit", py::arg("inputdata"), py::arg("gendata"), py::arg("dRLimit"), py::arg("option"));
    py::class_<Selection>(m, "Selection")
        .def(py::init<>())
        .def_readwrite("chi2", &Selection::chi2)
        .def_readwrite("bestPermutation", &Selection::bestPermutation)
        .def_readwrite("fitJets", &Selection::fitJets)
        .def_readwrite("combinationType", &Selection::combinationType);
}
