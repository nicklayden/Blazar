import adaptive
import numpy as np
import matplotlib.pyplot as plt



def ring(xy, wait=True):
    import numpy as np
    from time import sleep
    from random import random
    if wait:
        sleep(random()/10)
    x, y = xy
    a = 0.2
    return x + np.exp(-(x**2 + y**2 - 0.75**2)**2/a**4)

learner = adaptive.Learner2D(ring, bounds=[(-1, 1), (-1, 1)])
runner = adaptive.Runner(learner, goal=lambda l: l.loss() < 0.01)
def plot(learner):
    plot = learner.plot(tri_alpha=0.2)
    return (plot.Image + plot.EdgePaths.I + plot).cols(2)
#runner.ive_plot(plotter=plot, update_interval=0.1)
learner.plot(tri_alpha=0.4)