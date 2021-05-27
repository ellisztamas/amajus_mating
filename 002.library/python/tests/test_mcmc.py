import pytest
import numpy as np

from amajusmating import mcmc

# Dictionary listing starting values.
initial_model = {
    'missing' : 0.3, # proportion missing fathers
    'shape'  : 1,
    'scale'  : 10,
    'mixture' : 0.8
}

# Proposed values are a Gaussian peturbation away from the previous values.
# This is controlled by the sigma of the gaussian, which is defined for each variable
proposal_sigma = {
    'missing' : 0.025,
    'shape'  : 0.05,
    'scale'  : 2,
    'mixture' : 0.025,
}

def test_check_parmeter_output():
    initial_model = {
        'missing' : 0.3, # proportion missing fathers
        'shape'  : 1,
        'scale'  : 10,
        'mixture' : 0.8
    }
    out = mcmc.check_parameters(initial_model)
    assert isinstance(out, dict)
    assert out.keys == initial_model.keys()
    assert("assortment" not in out2.keys())

    model2 = {
        'missing' : 0.3, # proportion missing fathers
        'shape'  : 1,
        'scale'  : 10,
        'mixture' : 0.8,
        'assortment' : 0.5
    }
    out2 = mcmc.check_parameters(model2)
    assert("assortment" in out2.keys())

    # Should raise ValueError because missing > 1
    model3 = {
        'missing' : 1.3, # proportion missing fathers
        'shape'  : 1,
        'scale'  : 10,
        'mixture' : 0.8,
        'assortment' : 0.5
    }
    with pytest.raises(Exception):
        mcmc.check_parameters(model3)