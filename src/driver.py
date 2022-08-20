""" main function for driving stochastic shallow water model. """
import sys
from driver_model_read_input import model_read_input_driver
from driver_model_initiate import model_initiate_driver
from driver_model_formulate import model_formulate_driver
from driver_model_run import model_run_driver
from fenics import Timer
sys.setrecursionlimit(70000)


# TODO: design a User Interface for 2D-SSWM based on this input file.


def main():

    time1 = Timer("model_read_input")
    time1.start()
    inputs = model_read_input_driver()
    time1.stop()
    time1.elapsed()

    time2 = Timer("model_initiate")
    time2.start()
    initiate = model_initiate_driver(inputs)
    time2.stop()
    time2.elapsed()

    time3 = Timer("model_formulate")
    time3.start()
    formulate = model_formulate_driver(inputs, initiate)
    time3.stop()
    time3.elapsed()

    time4 = Timer("model_run_driver")
    time4.start()
    model_run_driver(inputs, initiate, formulate)
    time4.stop()
    time4.elapsed()


if __name__ == "__main__":
    main()
