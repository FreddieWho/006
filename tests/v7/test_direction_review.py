import numpy as np
from v7.direction_review import fit_direction


def test_direction_does_not_use_held_out_baseline_or_labels():
    pre=np.array([[0.,1.],[2.,3.],[3.,1.],[5.,4.],[100.,200.]])
    labels=np.array([0,1,0,1,1]); train=np.array([0,1,2,3])
    first=fit_direction(pre,labels,train)
    pre[-1]=[-300.,999.];labels[-1]=0
    second=fit_direction(pre,labels,train)
    for a,b in zip(first,second):
        np.testing.assert_array_equal(a,b)
