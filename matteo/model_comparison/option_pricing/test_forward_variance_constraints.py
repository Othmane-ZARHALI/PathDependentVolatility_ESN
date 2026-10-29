import unittest
from unittest.mock import patch
import numpy as np
from . import calibration_policy as policy
from .forward_variance_constraints import validate_forward_variance_levels


class ForwardVarianceConstraintTests(unittest.TestCase):
    def test_any_piece_count_and_nonmonotonic_levels_allowed(self):
        for count in [1, 3, 7, 14]:
            levels=np.resize([.04,.025,.06,.0036], count)
            np.testing.assert_array_equal(validate_forward_variance_levels(levels,.0036),levels)

    def test_equality_is_allowed_but_any_below_floor_is_rejected(self):
        floor=.0025
        np.testing.assert_array_equal(validate_forward_variance_levels([floor]*14,floor),[floor]*14)
        for index in range(14):
            levels=np.full(14,.04);levels[index]=np.nextafter(floor,-np.inf)
            with self.assertRaisesRegex(ValueError,'Every forward-variance level'):
                validate_forward_variance_levels(levels,floor)

    def test_invalid_inputs_rejected(self):
        for levels,floor in [([],.0025),([[.04]],.0025),([np.nan],.0025),([np.inf],.0025),
                             ([.04],0),([.04],-.0025),([.04],np.nan),([.04],[.0025])]:
            with self.subTest(levels=levels,floor=floor):
                with self.assertRaises(ValueError):validate_forward_variance_levels(levels,floor)

    def test_current_calibration_enforces_relation_independently_of_box_bounds(self):
        p=(policy.LOWER+policy.UPPER)/2;p[10]=policy.UPPER[10]
        # Current xi box bounds are stricter. Widen them temporarily to isolate
        # the hard xi>=floor check, as a future curve extension could do.
        lower=policy.LOWER.copy();lower[7:10]=policy.LOWER[10]
        with patch.object(policy,'LOWER',lower):
            for i in range(7,10):
                trial=p.copy();trial[i]=p[10]
                policy.validate(trial)
                trial[i]=np.nextafter(p[10],-np.inf)
                with self.assertRaisesRegex(ValueError,'Every forward-variance level'):
                    policy.validate(trial)


if __name__=='__main__':unittest.main()
