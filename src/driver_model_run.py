""" function to solve the stochastic shallow water equation. """

from model_run import ModelRun


def model_run_driver(inputs, initiate, formulate):

    run = ModelRun(inputs, initiate, formulate)

    run.run_prep()

    print "okay to here"
    run.running()

    run.run_final()

    return 0
