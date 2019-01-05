from rootpy.plotting.utils import _limits_helper, Compare
from rootpy.plotting import Hist, Graph, Profile
from nose.tools import assert_equal, assert_raises
import random


def test_limits():
    assert_equal(_limits_helper(0, 1, 0, 0), (0, 1))
    assert_equal(_limits_helper(1, 1, 0, 0, snap=True), (0, 1))
    assert_equal(_limits_helper(-2, -1, 0, 0, snap=True), (-2, 0))
    assert_equal(_limits_helper(-1, 1, .1, .1, snap=True), (-1.25, 1.25))

def test_Compare():
    #Histograms
    h1 = Hist(100,0,10,name='h1')
    h2 = Hist(100,0,10,name='h2')
    h3 = Hist(100,0,10,name='h3')

    h1.fill_random('gaus',5000)
    h2.fill_random('gaus',4000)
    h3.fill_random('gaus',6000)

    h12 = h1/h2
    h13 = h1/h3

    comp = Compare({'x0':h1, 'x1':[h2,h3]}, 'x0/x1')

    comp_h12 = comp[h1,h2][0]
    comp_h13 = comp[h1,h3][0]

    diff1 = [i for i in (h12-comp_h12) if i.value != 0]
    diff2 = [i for i in (h13-comp_h13) if i.value != 0]

    assert_equal(len(diff1), 0)
    assert_equal(len(diff2), 0)
    #End of histogram checks

    #Graphs
    g1 = Graph(20,name='g1')
    for i,point in enumerate(g1):
        point.x.value = 3*i
        point.x.error_low = 1.0
        point.x.error_hi = 2.0
        point.y.value = random.gauss(50,4)
        point.y.error_hi = random.gauss(3,1)
        point.y.error_low = random.gauss(3,1)

    g2 = Graph(20,name='g2')
    for i,point in enumerate(g2):
        point.x.value = 3*i
        point.x.error_low = 1.0
        point.x.error_hi = 2.0
        point.y.value = random.gauss(50,4)
        point.y.error_hi = random.gauss(3,1)
        point.y.error_low = random.gauss(3,1)

    g12=(g1-g2)/g1

    comp = Compare({'x0':g1, 'x1':g2}, '(x0-x1)/x0')

    comp_g12 = comp.results[0]

    diff = [i for i in (g12-comp_g12) if i.y.value != 0]

    assert_equal(len(diff), 0)
    #End of graph checks



if __name__ == "__main__":
    import nose
    nose.runmodule()
