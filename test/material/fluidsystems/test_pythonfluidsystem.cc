#include <config.h>

#include <dumux/python/material/fluidsystems/1ppython.hh>
#include <dumux/material/components/constant.hh>
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/embed.h>

int main()
{
    using FS = Dumux::FluidSystems::Python::OnePLiquid<double, Dumux::Components::Constant<0, double>>;
    std::cout << "density is " << FS::density(300.0, 1e5) << std::endl;

    return 0;
}
