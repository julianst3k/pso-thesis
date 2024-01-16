from setuptools import setup, find_packages
from glob import glob
from pybind11.setup_helpers import Pybind11Extension, ParallelCompile, build_ext
ext_modules = [
    Pybind11Extension(
        "model",
        sorted(glob("src/ChannelCPP/*.cpp")),
    ),
]
print(sorted(glob("src/ChannelCPP/*.cpp")))
ParallelCompile("NPY_NUM_BUILD_JOBS").install()
NAME = "PSO_Thesis"
version = '1.0'

print(find_packages())
setup(name=NAME, version=version, packages=find_packages('src'),
cmdclass={"build_ext": build_ext}, ext_modules=ext_modules, package_dir={"":"src"})