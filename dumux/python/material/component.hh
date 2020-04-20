#ifndef DUMUX_PYTHON_MATERIAL_COMPONENT_HH
#define DUMUX_PYTHON_MATERIAL_COMPONENT_HH

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/stl.h>

#include <dumux/material/components/liquid.hh>
#include <dumux/material/components/gas.hh>
#include <dumux/material/components/solid.hh>
#include <dumux/material/components/ion.hh>

namespace Dumux::Python {



template <class Comp, class... options>
void registerComponent(pybind11::handle scope,
                          pybind11::class_<Comp, options...> cls)
{
    using pybind11::operator""_a;

    cls.def(pybind11::init([](){
        return new Comp();
    }));

    using Scalar = typename Comp::Scalar;

    static constexpr bool isLiquid = std::is_base_of<Dumux::Components::Liquid<Scalar, Comp>, Comp>::value;
    static constexpr bool isGas = std::is_base_of<Dumux::Components::Gas<Scalar, Comp>, Comp>::value;
    static constexpr bool isSolid = std::is_base_of<Dumux::Components::Solid<Scalar, Comp>, Comp>::value;
    static constexpr bool isIon = std::is_base_of<Dumux::Components::Ion<Scalar, Comp>, Comp>::value;

    cls.def_static("name", &Comp::name);
    cls.def_static("molarMass", &Comp::molarMass);


    if constexpr (isLiquid)
    {
        cls.def_static("liquidDensity", &Comp::liquidDensity, "temperature"_a, "pressure"_a);
        cls.def_static("liquidMolarDensity", &Comp::liquidMolarDensity, "temperature"_a, "pressure"_a);
        cls.def_static("liquidIsCompressible", &Comp::liquidIsCompressible);
        cls.def_static("liquidViscosity", &Comp::liquidViscosity, "temperature"_a, "pressure"_a);
        cls.def_static("liquidEnthalpy", &Comp::liquidEnthalpy, "temperature"_a, "pressure"_a);
        cls.def_static("liquidInternalEnergy", &Comp::liquidInternalEnergy, "temperature"_a, "pressure"_a);
        cls.def_static("liquidHeatCapacity", &Comp::liquidHeatCapacity, "temperature"_a, "pressure"_a);
        cls.def_static("liquidThermalConductivity", &Comp::liquidThermalConductivity, "temperature"_a, "pressure"_a);
        cls.def_static("vaporPressure", &Comp::vaporPressure, "temperature"_a);
    }

    if constexpr (isGas)
    {
        cls.def_static("gasDensity", &Comp::gasDensity, "temperature"_a, "pressure"_a);
        cls.def_static("gasMolarDensity", &Comp::gasMolarDensity, "temperature"_a, "pressure"_a);
        cls.def_static("gasIsCompressible", &Comp::gasIsCompressible);
        cls.def_static("gasIsIdeal", &Comp::gasIsIdeal);
        cls.def_static("gasPressure", &Comp::gasPressure, "temperature"_a, "pressure"_a);
        cls.def_static("gasViscosity", &Comp::gasViscosity, "temperature"_a, "pressure"_a);
        cls.def_static("gasEnthalpy", &Comp::gasEnthalpy, "temperature"_a, "pressure"_a);
        cls.def_static("gasInternalEnergy", &Comp::gasInternalEnergy, "temperature"_a, "pressure"_a);
        cls.def_static("gasHeatCapacity", &Comp::gasHeatCapacity, "temperature"_a, "pressure"_a);
        cls.def_static("gasThermalConductivity", &Comp::gasThermalConductivity, "temperature"_a, "pressure"_a);
    }

    if constexpr (isSolid)
    {
        cls.def_static("solidIsCompressible", &Comp::solidIsCompressible);
        cls.def_static("solidDensity", &Comp::solidDensity, "temperature_a");
        cls.def_static("solidThermalConductivity", &Comp::solidThermalConductivity, "temperature_a");
        cls.def_static("solidHeatCapacity", &Comp::solidHeatCapacity, "temperature_a");
    }

    if constexpr (isIon)
    {
        cls.def_static("charge", &Comp::charge);
    }

    if constexpr (isLiquid || isGas)
    {
        cls.def_static("criticalTemperature", &Comp::criticalTemperature);
        cls.def_static("criticalPressure", &Comp::criticalPressure);
    }
}


} // namespace Dumux::Python

#endif
