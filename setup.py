from distutils.core import setup, Extension

pys = ['orngCRS',
       'orng2Array',
       'orngCluster',       
       'orngLR',
       'orngMultiClass',
       'orngSVM']

src = ['src/clustind.i',
       'src/twins.c',
       'src/pam.c',
       'src/fanny.c',
       'src/logreg.cpp',
       'src/lsq.cpp',
       'src/svm.cpp',
       'src/interface.cpp']

docs = ['doc/orngCluster.htm',
        'doc/orngLR-SVM-MultiClass.htm',
        'doc/orngClustRef.ps']

setup(name='orngExtn',
      version='1.8',
      description='Data Mining in Python',
      author='Aleks Jakulin',
      author_email='jakulin@ieee.org',
      url='http://ai.fri.uni-lj.si/~aleks/orng/',
      py_modules=pys,
      ext_modules=[Extension('_orngCRS', src, include_dirs=["include"])],
      data_files=[('doc', docs)])

