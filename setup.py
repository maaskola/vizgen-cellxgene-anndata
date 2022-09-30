import os
import subprocess

program_name = 'vizgen'
package_dir = 'vizgen'
version = subprocess.check_output(["git", "describe", "--always", "--dirty"]).strip().decode()

with open(
    os.path.join(package_dir, "__version__.py"), "w"
) as fp:
    fp.write(f'__version__ = "{version}"\n')
    fp.write(f'__program_name__ = "{program_name}"\n')

from setuptools import setup, find_packages

setup(
    name='vizgen',
    version=version,
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click',
        'anndata',
        'pandas',
        'scipy',
        'numpy',
        'matplotlib',
    ],
    entry_points=f'''
        [console_scripts]
        {program_name}-cellxgene_to_h5ad={package_dir}.cellxgene_to_h5ad:run
    ''',
    author='Jonas Maaskola',
    author_email='jonas.maaskola@weizmann.ac.il',
    description='Analyze Vizgen MERFISH data',
    url='https://github.com/maaskola/vizgen',
)
