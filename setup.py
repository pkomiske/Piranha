# Piranha Package
#
#  Questions/comments? pkomiske@mit.edu
#
#  Copyright (c) 2019-2021
#  Patrick T. Komiske III, Eric M. Metodiev, Jesse Thaler
#
#----------------------------------------------------------------------
# This file is part of FastJet contrib.
#
# It is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
#
# It is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this code. If not, see <http://www.gnu.org/licenses/>.
#----------------------------------------------------------------------

import os
import platform
import re
import subprocess
import sys

################################################################################

# Package name, with capitalization
name = 'Piranha'
lname = name.lower()

# using PyFJCore or not
use_pyfjcore = True

################################################################################

# use main fastjet library
if not use_pyfjcore:

    # function to query a config binary and get the result
    fastjet_config = os.environ.get('FASTJET_CONFIG', 'fastjet-config')
    def query_config(query):
        if not use_pyfjcore:
            return subprocess.check_output([fastjet_config, query]).decode('utf-8').strip()
        return ''

    # get fastjet info
    fj_prefix = query_config('--prefix')
    fj_cxxflags = query_config('--cxxflags')
    fj_ldflags = query_config('--libs')

if sys.argv[1] == 'swig':

    # form swig options
    if use_pyfjcore:
        opts = '-DPIRANHA_USE_PYFJCORE'
    else:
        opts = '-DFASTJET_PREFIX=' + fj_prefix + ' ' + fj_cxxflags

    command = ('swig -python -c++ -fastproxy -keyword -py3 -w509,511 -DSWIG_NUMPY {opts} '
               '-IEventGeometry -IEventGeometry/Wasserstein -IEventGeometry/PyFJCore '
               '-o {lname}/{lname}.cpp {lname}/swig/{lname}.i').format(opts=opts, lname=lname)
    print(command)
    subprocess.run(command.split())

else:

    import numpy as np
    from setuptools import setup
    from setuptools.extension import Extension

    # get contrib version
    with open(os.path.join(lname, '__init__.py'), 'r') as f:
        __version__ = re.search(r'__version__\s*=\s*[\'"]([^\'"]*)[\'"]', f.read()).group(1)

    # define containers of extension options
    sources = [os.path.join(lname, lname + '.cpp')]
    cxxflags = os.environ.get('CXXFLAGS', '').split()
    macros = []
    include_dirs = [np.get_include(), '.']
    ldflags = []
    library_dirs = []
    libraries = []

    # using main fastjet library
    if not use_pyfjcore:
        cxxflags += fj_cxxflags.split()
        macros.append(('SWIG_TYPE_TABLE', 'fastjet'))
        macros.append(('EVENTGEOMETRY_TEMPLATE_VISIBILITY', None))
        for ldflag in fj_ldflags.split():
            if ldflag.startswith('-L'):
                library_dirs.append(ldflag[2:])
            elif ldflag.startswith('-l'):
                libraries.append(ldflag[2:])
            else:
                ldflags.append(ldflag)

    # using pyfjcore
    else:
        cxxflags.append('-std=c++14')
        macros.append(('EVENTGEOMETRY_USE_PYFJCORE', None))
        macros.append(('PIRANHA_USE_PYFJCORE', None))
        macros.append(('SWIG_TYPE_TABLE', 'fjcore'))
        include_dirs.append('EventGeometry')
        include_dirs.append(os.path.join('EventGeometry', 'PyFJCore'))
        include_dirs.append(os.path.join('EventGeometry', 'Wasserstein'))

        # need to compile pyfjcore from scratch for windows
        if platform.system() == 'Windows':
            sources.append(os.path.join('PyFJCore', 'pyfjcore', 'fjcore.cc'))
            cxxflags = ['/std:c++14']

        # on non-windows we can use shared library
        else:
            macros.append(('DECLARE_EVENTGEOMETRY_TEMPLATES', None))
            library_dirs.extend(['.'])
            libraries.extend([name])

            if platform.system() == 'Darwin':
                ldflags.append('-Wl,-rpath,@loader_path/..')

            elif platform.system() == 'Linux':
                ldflags.append('-Wl,-rpath,$ORIGIN/..')

    # if not windows, further modification needed for multithreading
    if platform.system() != 'Windows':

        # no debugging
        cxxflags.append('-g0')

    module = Extension('{0}._{0}'.format(lname),
                       sources=sources,
                       define_macros=macros,
                       include_dirs=include_dirs,
                       library_dirs=library_dirs,
                       libraries=libraries,
                       extra_compile_args=cxxflags,
                       extra_link_args=ldflags
                      )

    setup(
        ext_modules=[module],
        version=__version__
    )
