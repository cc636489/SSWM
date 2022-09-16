

#from fenics import Constant, interpolate, project, compile_cpp_code, Expression, CompiledExpression, as_vector, sqrt
from dolfin import *
from make_sto_modes import make_sto_modes
from netCDF4 import Dataset
import numpy as np


def make_bath(bath_dict, pc_basis_str, bath_function_space):

    # check if the key are entered right.
    key = list(bath_dict.keys())
    n_modes = pc_basis_str.get("n_modes")
    
    if len(key) != 1:
        raise TypeError("enter more than one domain! not allowed! Exit program...")
    if key[0] != "flat" and key[0] != "vary" and key[0] != "class":
        raise TypeError("don't have this kind of bathymetry implemented! Exit program...")

    # check if the value are entered right.
    value = bath_dict.get(key[0])
    if key[0] == "flat":
        if not isinstance(value, float):
            raise TypeError("wrong bathymetry type entered! Could only input float! Exit program...")
    elif key[0] == "vary":
        if not isinstance(value, str):
            raise TypeError("wrong bathymetry type entered! Could only input string ! Exit program...")
    elif key[0] == "class":
        if not isinstance(value, list):
            raise TypeError("wrong bathymetry list entered! Could only be ['type*', ...] ! Exit program...")
        if value[0] != "type1" and value[0] != "type2" and value[0] != "type3":
            raise TypeError("other types of bathymetry class has not been implemented yet! Exit program...")
        if value[0] == "type1" and len(value) != 5:
            raise TypeError("type1 class should read in a list with length 5 Exit program...")
        if value[0] == "type2" and len(value) != 4:
            raise TypeError("type2 class should read in a list with length 4 Exit program...")
        if value[0] == "type3" and len(value) != 2:
            raise TypeError("type3 class should read in a generalized bathymetry file. Exit program...")
    else:
        raise TypeError("other type of bathymetry has not been implemented yet! Exit program...")

    # make stochastic version of bathymetry
    if key[0] == "flat":
        bath_list = make_sto_modes(pc_basis_str, str(value))
        if n_modes == 1:
            bath_expression = Constant(eval(bath_list[0]))
        else:
            bath_expression = Constant(bath_list)
        bath_function = interpolate(bath_expression, bath_function_space)
    elif key[0] == "vary":
        bath_list = make_sto_modes(pc_basis_str, value)
        if n_modes == 1:
            bath_expression = Expression(bath_list[0], element=bath_function_space.ufl_element())
        else:
            bath_expression = Expression(bath_list, element=bath_function_space.ufl_element())
        bath_function = interpolate(bath_expression, bath_function_space)
    elif key[0] == "class":
        if value[0] == "type1":
            bath_list_at3, bath_list_at4 = make_sto_modes(pc_basis_str, value[3], value[4])

            class BottomExpression(UserExpression):
                def eval(self, mode, x):
                    if x[0] < value[1] or x[0] > value[2]:
                        for ii in range(n_modes):
                            mode[ii] = eval(bath_list_at3[ii])
                    else:
                        for ii in range(n_modes):
                            mode[ii] = eval(bath_list_at4[ii])

                def value_shape(self):
                    return (n_modes,)

            bath_expression = BottomExpression(element=bath_function_space.ufl_element())
            bath_function = interpolate(bath_expression, bath_function_space)
        elif value[0] == "type2":
            bath_list_at3, bath_list_at4 = make_sto_modes(pc_basis_str, value[2], value[3])

            class BottomExpression2(UserExpression):
                def eval(self, mode, x):
                    if x[0] <= value[1]:
                        for ii in range(n_modes):
                            mode[ii] = eval(bath_list_at3[ii])
                    else:
                        for ii in range(n_modes):
                            mode[ii] = eval(bath_list_at4[ii])
                    hump_center_x = 1775.0
                    hump_center_y = 1500.0
                    hump_radius = 375.0
                    hump_height = 5.44186
                    if (x[0]-hump_center_x)**2 + (x[1]-hump_center_y)**2 <= hump_radius**2:
                        temp = sqrt((x[0]-hump_center_x)**2 + (x[1]-hump_center_y)**2)
                        mode[0] = mode[0] - (2.0*hump_height*temp**3/hump_radius**3 -
                                             3.0*hump_height*temp**2/hump_radius**2 + hump_height)

                def value_shape(self):
                    return (n_modes,)

            bath_expression = BottomExpression2(element=bath_function_space.ufl_element())
            bath_function = interpolate(bath_expression, bath_function_space)
        elif value[0] == "type3":
            code = """

#include <dolfin/function/Expression.h>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

class MyFun : public dolfin::Expression
{
    public:
        
        Eigen::VectorXd xcoordinate;
        Eigen::VectorXd ycoordinate;
        Eigen::VectorXd bathymetryValue;

    public:

        MyFun(): dolfin::Expression() {}

        void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
        {
            double vv = findval(x);
            values[0] = vv;
        }

        double findval(Eigen::Ref<const Eigen::VectorXd> x) const
        {
            for(int index=0; index < xcoordinate.size(); index++)
            {
                if ( (fabs(x[0] - xcoordinate[index])<1.0e-4) && (fabs(x[1] - ycoordinate[index])<1.0e-4) )
                {
                    return bathymetryValue[index];
                }
            }
            return 0;
        }
};

PYBIND11_MODULE(SIGNATURE, m)
{
  pybind11::class_<MyFun, std::shared_ptr<MyFun>, dolfin::Expression>
    (m, "MyFun")
    .def(pybind11::init<>())
    .def_readwrite("xcoordinate", &MyFun::xcoordinate)
    .def_readwrite("ycoordinate", &MyFun::ycoordinate)
    .def_readwrite("bathymetryValue", &MyFun::bathymetryValue)
    ;
}

"""
            nc = Dataset(value[1], 'r', format='NETCDF4')
            y_coord = nc.variables['y_coord'][:].data.tolist()
            x_coord = nc.variables['x_coord'][:].data.tolist()
            bath_depth = nc.variables['bathy'][:].data.tolist()
            if n_modes == 1:
                bath_expression = CompiledExpression(compile_cpp_code(code).MyFun(), element=bath_function_space.ufl_element())
                
                bath_expression.xcoordinate = np.array(x_coord, dtype=float)
                bath_expression.ycoordinate = np.array(y_coord, dtype=float)
                bath_expression.bathymetryValue = np.array(bath_depth, dtype=float)
                
                
                bath_function = interpolate(bath_expression, bath_function_space)
            else:
                bath_expression = CompiledExpression(compile_cpp_code(code).MyFun(), element=bath_function_space.sub(0).ufl_element())
                
                bath_expression.xcoordinate = np.array(x_coord, dtype=float)
                bath_expression.ycoordinate = np.array(y_coord, dtype=float)
                bath_expression.bathymetryValue = np.array(bath_depth, dtype=float)
                
                bath_string = "[bath_expression"
                for _ in range(n_modes-1):
                    bath_string += ", 0"
                bath_string += "]"
                bath_form = as_vector(eval(bath_string))
                bath_function = project(bath_form, bath_function_space)
        else:
            raise TypeError("don't implement this class bathymetry type.")
    else:
        raise TypeError("don't implement this key bathymetry type.")

    return bath_function
