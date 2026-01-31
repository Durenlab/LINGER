from setuptools import setup

setup(
    name='LingerGRN',
    version='1.106',
    description='Gene regulatory network inference',
    author='Kaya Yuan',
    author_email='qyyuan33@gmail.com',
    packages=['LingerGRN'],
    license = "MIT",
    url='https://github.com/Durenlab/LINGER',
    install_requires=['torch', 'scipy==1.11.3', 'numpy==1.24.3', 'pandas==2.0.3', 'shap==0.42.0', 'scikit-learn==1.3.0', 'joblib==1.3.2','matplotlib==3.8.0','seaborn==0.13.0','statsmodels==0.14.1','umap-learn','scanpy==1.9.5','anndata==0.9.2','pybedtools==0.10.0'],
)
