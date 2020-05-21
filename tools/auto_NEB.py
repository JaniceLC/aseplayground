#!/usr/bin/env python3
"""Analyzes the output of NEB scripts from the Atomic Simulation
Environment (ASE)."""

from ase import io
import ase.io
from matplotlib import pyplot
from tempfile import mkdtemp
from optparse import OptionParser
import shutil
import os
import matplotlib
matplotlib.use('Agg')  # Needed for headless operation.


try:  # ASE < svn 4349.
    from ase.neb import get_NEB_plot
except ImportError:  # ASE >= svn 4349.
    from ase.neb import NEBtools

    def get_NEB_plot(images):
        return NEBtools(images).plot_band()


def find_trajectories(startswith):
    """Finds NEB trajectories."""
    trajs = [file for file in os.listdir(os.getcwd()) if
             (file.startswith(startswith) and file.endswith('.traj'))]
    trajs.sort()
    if len(trajs) == 0:
        raise RuntimeError('No NEB trajectories found. Consider the '
                           '-s flag.')
    print('Found trajectories: %s' % trajs)
    return trajs


def makePDFs(trajs, n_images, output):
    """Creates PDFs of each band in the trajectory files, assuming n_images
    images per band. Then concatenates them all into a single pdf called
    <output>.pdf and deletes the source files."""
    if n_images == -1:
        n_images = guess_nimages(ase.io.read(trajs[0], index=':'))
    count = 0
    filenames = []
    tmpd = mkdtemp()  # temporary directory
    for file in trajs:
        print('Analyzing: %s' % file)
        traj = io.Trajectory(file, 'r')
        assert (len(traj) % n_images) == 0
        steps = len(traj) // n_images
        for step in range(steps):
            filenames.append(os.path.join(tmpd, 'neb_plot%04i.pdf' % count))
            indices = range(step * n_images, step * n_images + n_images)
            print('%i/%i  %s: %s' % (step + 1, steps, filenames[-1], indices))
            images = [traj[index] for index in indices]
            fig = get_NEB_plot(images)
            fig.savefig(filenames[-1])
            pyplot.close(fig)  # garbage collection
            count += 1

    command = 'pdftk ' + ('%s ' * len(filenames)) + 'cat output %s.pdf'
    command = command % tuple(filenames + [output])
    print('Combining pdfs with pdftk to %s.' % output)
    os.system(command)

    shutil.rmtree(tmpd)


def guess_nimages(images):
    """Attempts to guess the number of images per band from
    a trajectory file, based solely on the repetition of the
    potential energy of images. This might fail for symmetric
    cases."""
    e_first = images[0].get_potential_energy()
    for index, image in enumerate(images[1:], start=1):
        e = image.get_potential_energy()
        if e == e_first:
            n_images = index
            break
    # Sanity check that both first and last line up.
    e_first = images[0].get_potential_energy()
    e_nextfirst = images[n_images].get_potential_energy()
    e_last = images[n_images - 1].get_potential_energy()
    e_nextlast = images[2 * n_images - 1].get_potential_energy()
    if not (e_first == e_nextfirst) and (e_last == e_nextlast):
        raise RuntimeError('Could not guess number of images per band.')
    print('Number of images guessed to be {:d}.'.format(n_images))
    return n_images


if __name__ == '__main__':
    parser = OptionParser(usage='usage: %prog [options]')
    parser.add_option('-i', dest='images', default=-1,
                      help='number of images per band, '
                           'guessed if not supplied')
    parser.add_option('-s', dest='startswith', default='neb',
                      help='NEB traj files start with, default="neb"')
    parser.add_option('-o', dest='output', default='all-NEBs',
                      help='filename (without suffix) for output, '
                           'default="all-NEBs"')
    options, args = parser.parse_args()

    trajs = find_trajectories(options.startswith)
    makePDFs(trajs, int(options.images), options.output)
