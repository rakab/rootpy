#!/usr/bin/env python
"""
===========================================
Comparing multiple histograms to each other
===========================================

This example demonstrates how to compare multiple histograms to each other in
a convenient way. Same technique can be used for TGraph and TProfile as well.
"""
print(__doc__)
import rootpy.ROOT as R
R.gStyle.SetOptStat(0)
R.gStyle.SetOptTitle(0)
from rootpy.extern.six.moves import range
from rootpy.plotting import Hist, Canvas, Pad, Legend
from rootpy.plotting.utils import draw, Compare
from rootpy.plotting.style import get_style
from rootpy.interactive import wait
import random

# create couple 1D histograms to compare
# bin number and bin edges should be the same for all of them
h_data1 = Hist(100, 0, 100, name='data1', title='Data 1')
h_data2 = Hist(100, 0, 100, name='data2', title='Data 2', linestyle='dashed',markerstyle=21)

h_mc1 = Hist(100, 0, 100, name='mc1', title='MC 1', color='red')
h_mc2 = Hist(100, 0, 100, name='mc2', title='MC 2', color='blue')
h_mc3 = Hist(100, 0, 100, name='mc3', title='MC 3', color='magenta')

h_w = Hist(100, 0, 100, name='w', title='Weights', color='green')

# fill histograms
for i in range(10000):
    h_data1.fill(random.gauss(50,20))
    h_data2.fill(random.gauss(51,24))
    h_mc1.fill(random.gauss(50,30))
    h_mc2.fill(random.gauss(47,23))
    h_mc3.fill(random.gauss(55,17))

#Uniformly distibuted weight between 0 and 1
for i in h_w:
    i.value = random.random()


#Assign histos to keys for late use as varaibles in comparisson function
#definiton
input_plotables = {'x0':[h_data1, h_data2], 'x1':[h_mc1, h_mc2, h_mc3], 'x2':h_w}

#Compare will subtract weighted mc from data and divide the result on errors of
#data
output = Compare(input_plotables,'(x0-x2*x1)/dx0', keytoclone='x1')
#Resulting histograms will have exactly the same style what x1 -
#(h_mc1, h_mc2, h_mc3) histos have

#For histos which are compared to h_data2 set markerstyle from h_data2
for h in output[h_data2]:
    h.markerstyle = h_data2.markerstyle

#To access list of all resulting histos use results attribute of Compare class
results = output.results

#To check how errors of resulting histos are being calculated use errfunc
#attribute of Compare class
print(output.errfunc)


#Draw inital and resulting histos togather
c = Canvas(width=800,height=800)
p_t = Pad(0.01,0.35,0.95,0.97)
p_b = Pad(0.01,0.0,0.95,0.35)
p_t.set_bottom_margin(0.01)
p_b.set_top_margin(0)
p_b.set_bottom_margin(0.3)
p_t.draw()
p_b.draw()

#Upper part
h_data1.yaxis.set_label_size(0.05)
draw([h_data1,h_data2,h_mc1,h_mc3], pad=p_t)
legend = Legend([h_data1,h_data2,h_mc1,h_mc2,h_mc3], pad=p_t,
                    leftmargin=0.05, rightmargin=0.5, textsize=0.02, entryheight=0.02)
legend.Draw()

#Lower part
results[0].yaxis.set_label_size(0.1)
results[0].xaxis.set_label_size(0.1)
draw(results, pad=p_b)

c.draw()

# wait for you to close all open canvases before exiting
# wait() will have no effect if ROOT is in batch mode:
# ROOT.gROOT.SetBatch(True)
wait()
