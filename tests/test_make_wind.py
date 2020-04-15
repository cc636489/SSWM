

import sys
sys.path.insert(0, '/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/src')
sys.path.insert(0, '/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input')

import unittest
from make_wind import MakeWind
from model_read_input import ModelReadInput
from fenics import Mesh, FunctionSpace, Expression, XDMFFile, project, Function
from netCDF4 import Dataset
import numpy as np


class MakeMeshTestCase(unittest.TestCase):

    def setUp(self):
        self.inputs = ModelReadInput()

    def test_wind(self):

        # now we can update the wind speed in every time step.
        days = 2.58 
        curr_time = 0
        final_time = 3600 * 24 * days
        dt = 447


        this_mesh = Mesh(self.inputs.input_dir + self.inputs.mesh_file)
        this_scalar_space = FunctionSpace(this_mesh, "CG", 1)
        code = '''
    
        namespace dolfin {
    
        struct node {
            double x;
            double y;
            double v;
            node *next;
        };
    
        class MyFun : public Expression
        {
            private:
                node *root;
    
            public:
    
                MyFun(): Expression()
                {
                    root = new node;
                    root->next = 0;
                };
    
                void initial(double _x, double _y, double _v)
                {
                    node *newNode = new node;
                    newNode->x = _x;
                    newNode->y = _y;
                    newNode->v = _v;
                    newNode->next = root;
                    root = newNode;
                    //cout << _v << endl;
                };
    
                void update(double _x, double _y, double _v)
                {
                    node *p = findp(_x, _y);
                    //cout << p->v << "  " << _v << endl;
                    p->v = _v;
                };
    
                void eval(Array<double>& values, const Array<double>& x) const
                {
                    double vv = findval(x[0], x[1]);
                    values[0] = vv;
                    //cout << x[0] << "  " << x[1] << "  " << vv << endl;
                    //cout << "  " << endl;
                };
    
                node * findp(double _x, double _y) const
                {
                    node *p = root;
                    while (p->next != 0)
                    {
                        if ( (fabs(p->x - _x)<1.0e-4) && (fabs(p->y - _y)<1.0e-4) )
                        {
                            return p;
                        }
                        else
                        {
                            p = p->next;
                        }
                    }
                    return 0;
                }
    
                double findval(double _x, double _y) const
                {
                    node *p = root;
                    while (p->next != 0)
                    {   
                        if ( (fabs(p->x - _x)<1.0e-4) && (fabs(p->y - _y)<1.0e-4) )
                        {
                            //cout << fabs(p->x-_x) << "  " << fabs(p->y-_y) << endl;
                            double find = p->v;
                            return find;
                        }
                        else
                        {
                            p = p->next;
                        }
                    }
                    return 0;
                };
    
        };
        };
    
        '''

        # create a list for wind_para_x, wind_para_y.
        wind_para_x = Expression(cppcode=code, element=this_scalar_space.ufl_element())
        wind_para_y = Expression(cppcode=code, element=this_scalar_space.ufl_element())
        wdrag = Expression(cppcode=code, element=this_scalar_space.ufl_element())

        # read in variables.
        nc = Dataset(self.inputs.input_dir + self.inputs.bath_file, 'r', format='NETCDF4')
        x_coord = nc.variables['x_coord'][:].data.tolist()
        y_coord = nc.variables['y_coord'][:].data.tolist()
        x_deg = nc.variables['lon'][:].data.tolist()
        y_deg = nc.variables['lat'][:].data.tolist()

        # initialize a list for wind_para_x, wind_para_y.
        for st in range(len(x_coord)):
            wind_para_x.initial(x_coord[st], y_coord[st], 0.0)
            wind_para_y.initial(x_coord[st], y_coord[st], 0.0)
            wdrag.initial(x_coord[st], y_coord[st], 0.00075)

        wind = MakeWind(self.inputs)
        wind_x_max = []
        wind_y_max = []
        wdrag_max = []

        f = XDMFFile("wind_drag_test.pvd")
        drag_func = Function(this_scalar_space, name="drag")

        while curr_time <= final_time:

            wind.get_prepared(current_time=curr_time)
            wind_para_x_list = []
            wind_para_y_list = []
            wind_drag_list = []
            for t in range(len(x_deg)):
                wind_x, wind_y, pressure = wind.get_wind_x_wind_y(x_coord_deg=x_deg[t], y_coord_deg=y_deg[t])
                wind_para_x_list.append(wind_x)
                wind_para_y_list.append(wind_y)
                wind_drag = wind.get_wind_drag(x_coord_deg=x_deg[t], y_coord_deg=y_deg[t])
                wind_drag_list.append(wind_drag)

            wind_x_max.append(np.max(wind_para_x_list))
            wind_y_max.append(np.max(wind_para_y_list))
            wdrag_max.append(np.max(wind_drag_list))

            for sm in range(len(x_coord)):
                wind_para_x.update(x_coord[sm], y_coord[sm], wind_para_x_list[sm])
                wind_para_y.update(x_coord[sm], y_coord[sm], wind_para_y_list[sm])
                wdrag.update(x_coord[sm], y_coord[sm], wind_drag_list[sm])

            # save wdrag to file.
            
            drag_func.assign(project(wdrag, this_scalar_space))
            f << (drag_func, float(curr_time))

            # update current time
            print curr_time/dt
            curr_time += dt


        self.assertTrue(1)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
