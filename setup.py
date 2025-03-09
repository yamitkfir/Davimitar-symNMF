from setuptools import Extension, setup

module = Extension("symnmfmodule", sources=['symnmfmodule.c'])#, 'symnmfalgo.c'])
setup(name='symnmfmodule',
     version='1.0',
     description='Python wrapper for custom C extension',
     ext_modules=[module])

# install by running in terminal:
# python3 setup.py build_ext --inplace