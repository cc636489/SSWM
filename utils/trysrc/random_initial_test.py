from fenics import list_timings, TimingClear_keep, TimingType_wall, TimingType_system, Timer, RectangleMesh, FunctionSpace, Expression, Point
from random import seed, random

def main():
    m = Timer("dd")
    n_modes = 2
    mesh = RectangleMesh(Point(0, 0), Point(100, 20), 100, 20)
    V = FunctionSpace(mesh, 'CG', 1)
    m.stop()
    
    class RandomInitialConditions(Expression):
        def __init__(self, **kwargs):
            self.n_modes = kwargs["n_modes"]
            self.element = kwargs["element"]
            seed(self.n_modes)
    
        def eval(self, values, x):
            for i in range(self.n_modes):
                values[i] = - 1e-6 + (random() - 0.1) / 0.9 * 2 * 1e-6
        def value_shape(self):
            return (self.n_modes,)
    t = Timer("cc")
    a = RandomInitialConditions(n_modes=n_modes, element=V.ufl_element())
    t.stop()
    
    list_timings(TimingClear_keep, [TimingType_wall, TimingType_system])

if __name__ == "__main__":
    main()
