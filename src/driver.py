""" main function for driving stochastic shallow water model. """

from driver_model_read_input import model_read_input_driver
from driver_model_initiate import model_initiate_driver
from driver_model_formulate import model_formulate_driver
from driver_model_run import model_run_driver


import sys
sys.path.insert(0, '/workspace/Documentation/Research_Doc/SFEM_Doc/7-NS-github/input')
sys.setrecursionlimit(10000)


def main():

    inputs = model_read_input_driver()

    initiate = model_initiate_driver(inputs)

    formulate = model_formulate_driver(inputs, initiate)

    model_run_driver(inputs, initiate, formulate)


if __name__ == "__main__":
    main()
