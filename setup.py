import os
from setuptools import find_packages, setup
from setuptools.command.install import install


def find_scripts():
    all_scripts = []
    for root, dirs, files in os.walk("."):
        for file in files:
            if file.endswith(".sh") and (not("test" in root or "other_classifiers" in root) or file == "HiTaxon.sh"):
                all_scripts.append(os.path.join(root, file))
    return all_scripts

class CustomInstallCommand(install):
    def run(self):
        install.run(self)
        all_scripts = find_scripts()
        for script_path in all_scripts:
            os.system(f"chmod +x {script_path}")

setup(
    name='HiTaxon',
    packages=find_packages(),
    version=1.0,
    description='HiTaxon: A hierarchical ensemble framework for taxonomic classification of short reads',
    url='https://github.com/ParkinsonLab/HiTaxon',
    author='Bhavish Verma',
    author_email='bhavish.verma@mail.utoronto.ca',
    license='MIT',

    cmdclass={'install': CustomInstallCommand},
    # Add any other configuration here, e.g., version, description, author, etc.
)
