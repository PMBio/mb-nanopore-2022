import setuptools

setuptools.setup(
    name="mb_analysis",
    version="0.1.0",
    author="Rene Snajder",
    author_email="r.snajder@dkfz-heidelberg.de",
    description="ASE and ASM analysis for MB project",
    url="https://github.com/snajder-r/medulloblastoma_analysis´",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "meth5",
        "tqdm",
        "nanoepitools",
        "mygene",
        "scipy",
    ],
    python_requires='>=3.7',
)
