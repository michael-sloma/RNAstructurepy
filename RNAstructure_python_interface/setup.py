#!/usr/bin/env python

"""
setup.py file for RNAstructure python interface
"""

from distutils.core import setup, Extension
from distutils.sysconfig import get_config_vars
import os

#this removes the annoying -Wstrict_prototypes warning
(opt,) = get_config_vars('OPT')
os.environ['OPT'] = " ".join(
    flag for flag in opt.split() if flag != '-Wstrict-prototypes'
)



RNAstructure_module = Extension('_RNAstructure',
                           sources=['RNAstructure_wrap.cxx',
                                    '../RNA_class/Dynalign_object.cpp',
                                    '../RNA_class/HybridRNA.cpp',
                                    '../RNA_class/Multilign_object.cpp',
                                    '../RNA_class/Oligowalk_object.cpp',
                                    '../RNA_class/RNA.cpp',
                                    '../RNA_class/thermodynamics.cpp',
                                    '../RNA_class/TwoRNA.cpp',
                                    '../RNA_class/ProbScan.cpp',
                                    '../src/algorithm.cpp',
                                    '../src/alltrace.cpp',
                                    '../src/alltrace_intermolecular.cpp',
                                    '../src/arrayclass.cpp',
                                    '../src/bimol.cpp',
                                    '../src/dotarray.cpp',
                                    '../src/DotPlotHandler.cpp',
                                    '../src/draw.cpp',
                                    '../src/dynalign.cpp',
                                    '../src/dynalignarray.cpp',
                                    '../src/dynalignheap.cpp',
                                    '../src/dynalignstackclass.cpp',
                                    '../src/extended_double.cpp',
                                    '../src/forceclass.cpp',
                                    '../src/intermolecular.cpp',
                                    '../src/MaxExpect.cpp',
                                    '../src/MaxExpectStack.cpp',
                                    '../src/observable.cpp',
                                    '../src/observer.cpp',
                                    '../src/OligoScreenCalc.cpp',
                                    '../src/outputconstraints.cpp',
                                    '../src/pclass.cpp',
                                    '../src/pfunction.cpp',
                                    '../src/probknot.cpp',
                                    '../src/random.cpp',
                                    '../src/rna_library.cpp',
                                    '../src/siRNAfilter.cpp',
                                    '../src/stackclass.cpp',
                                    '../src/stackstruct.cpp',
                                    '../src/stochastic.cpp',
                                    '../src/structure.cpp',
                                    '../src/StructureImageHandler.cpp',
                                    '../src/thermo.cpp',
                                    '../src/varray.cpp',
                                    '../src/wendarray.cpp',
                                    '../src/phmm/aln_env_utils.cpp',
                                    '../src/phmm/p_alignment.cpp',
                                    '../src/phmm/phmm_aln.cpp',
                                    '../src/phmm/phmm_array.cpp',
                                    '../src/phmm/phmm.cpp',
                                    '../src/phmm/phmm_ml_loops.cpp',
                                    '../src/phmm/phmm_pp_loops.cpp',
                                    '../src/phmm/structure/folding_constraints.cpp',
                                    '../src/phmm/structure/structure_object.cpp',
                                    '../src/phmm/utils/ansi_string/ansi_string.cpp',
                                    '../src/phmm/utils/file/utils.cpp',
                                    '../src/phmm/utils/xmath/linear/linear_math.cpp',
                                    '../src/phmm/utils/xmath/log/xlog_math.cpp',
                                    '../src/phmm/utils/xmath/matrix/matrix.cpp',
                                    '../TurboFold/TurboFold_object.cpp',
                                    '../TurboFold/TurboFold_thread.cpp',
                                    '../src/TProgressDialog.cpp',
                                    ],
                                    extra_compile_args=['-w','-O3'],
                           )

setup (name = 'RNAstructure',
       version = '1.0',
       author      = "Michael Sloma",
       description = """Python interface for the RNA library""",
       ext_modules = [RNAstructure_module],
       py_modules = ["RNAstructure"],
       )
