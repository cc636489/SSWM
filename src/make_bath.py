

#from fenics import Constant, interpolate, project, compile_cpp_code, Expression, CompiledExpression, as_vector, sqrt
from dolfin import *
from make_sto_modes import make_sto_modes
from netCDF4 import Dataset


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

            class BottomExpression(Expression):
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

            class BottomExpression2(Expression):
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
    struct node {
        double x;
        double y;
        double v;
        node *next;
    };
    
    private:
        node *root;

    public:

        MyFun(): dolfin::Expression()
        {
            root = new node;
            root->next = 0;
        }

        void update(double _x, double _y, double _v)
        {
            node *newNode = new node;
            newNode->x = _x;
            newNode->y = _y;
            newNode->v = _v;
            newNode->next = root;
            root = newNode;
        }

        void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
        {
            double vv = findval(x[0], x[1]);
            values[0] = vv;
            
            //cout << x[0] << "  " << x[1] << "  " << vv << endl;
            //cout << "  " << endl;
        }

        double findval(double _x, double _y) const
        {
            node *p = root;
            while (p->next != 0)
            {

            // Assume that root node has biggest x-value.
            // Traverse down the list until p->x = _x and p->y = _y, then assign p->v to _v, and return value, break.
                
                if ( (fabs(p->x - _x)<1.0e-4) && (fabs(p->y - _y)<1.0e-4) )
                {
                    // cout << fabs(p->x-_x) << "  " << fabs(p->y-_y) << " " << p->x << " " << _x << " " 
                    // << p->y << " " << _y << endl;
                
                    double find = p->v;
                    return find;
                }
                else
                {
                    p = p->next;
                }
            }
            return 0;

        }

};

PYBIND11_MODULE(SIGNATURE, m)
{
  pybind11::class_<MyFun, std::shared_ptr<MyFun>, dolfin::Expression>(m, "MyFun").def(pybind11::init<>());
}

"""
            nc = Dataset(value[1], 'r', format='NETCDF4')
            y_coord = nc.variables['y_coord'][:].data.tolist()
            x_coord = nc.variables['x_coord'][:].data.tolist()
            bath_depth = nc.variables['bathy'][:].data.tolist()
            if n_modes == 1:
                #bath = Expression(code, element=bath_function_space.ufl_element())
                #bath_expression = CompiledExpression(compile_cpp_code(code).MyFun(), element=bath_function_space.ufl_element())
                bath_expression = CompiledExpression(compile_cpp_code(code).MyFun(), degree=1)
                type(bath_expression)
                print("Chen I am okay here...", type(bath_expression))
                for jj in range(len(x_coord)):
                    bath_expression.update(x_coord[jj], y_coord[jj], bath_depth[jj])
                bath_function = interpolate(bath_expression, bath_function_space)
            else:
                bath = Expression(code, element=bath_function_space.sub(0).ufl_element())
                for jj in range(len(x_coord)):
                    bath.update(x_coord[jj], y_coord[jj], bath_depth[jj])
                bath_string = "[bath"
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
