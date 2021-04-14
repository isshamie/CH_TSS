from setuptools import find_packages, setup

setup(
    name='tss',
    packages=find_packages(),
    version='0.1.0',
    description='Transcription Start Sites annotation pipeline',
    author='Isaac Shamie',
    license='MIT',
    include_package_data=True,
    install_requires=[
            'Click',
        ],
    entry_points='''
            [console_scripts]
            gene_centric=tss.data.gene_centric:main
        '''
)
