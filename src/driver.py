""" main function for driving stochastic shallow water model. """

from driver_model_read_input import model_read_input_driver
from driver_model_initiate import model_initiate_driver
from driver_model_formulate import model_formulate_driver
from driver_model_run import model_run_driver
from fenics import list_timings, TimingClear_keep, TimingType_wall, TimingType_system, Timer
import sys
sys.setrecursionlimit(70000)


# TODO: design a User Interface for 2D-SSWM based on this input file.


def main():

    time1 = Timer("model_read_input")

    inputs = model_read_input_driver()

    time1.stop()

    time2 = Timer("model_initiate")

    initiate = model_initiate_driver(inputs)

    time2.stop()

    time3 = Timer("model_formulate")

    formulate = model_formulate_driver(inputs, initiate)

    time3.stop()

    time4 = Timer("model_run_driver")

    model_run_driver(inputs, initiate, formulate)

    time4.stop()

    list_timings(TimingClear_keep, [TimingType_wall, TimingType_system])


if __name__ == "__main__":
    main()
