from dolfin import *

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

        void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
        {
            double vv = findval(x[0], x[1]);
            values[0] = vv;
        }

        double findval(double _x, double _y) const
        {
            node *p = root;
            while (p->next != 0)
            {
                if ( (fabs(p->x - _x)<1.0e-4) && (fabs(p->y - _y)<1.0e-4) )
                {
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

haha = compile_cpp_code(code)
 
bath_expression = CompiledExpression(compile_cpp_code(code).MyFun(), degree=1)

print("Chen Chen I am here.")
