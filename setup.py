import setuptools

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setuptools.setup(name='marker_alignments',
    version='0.4.1',
    description="Process and summarise alignments of metagenomic sequencing reads to reference databases of marker genes",
    long_description = long_description,
    long_description_content_type="text/markdown",
    url='http://github.com/wbazant/marker_alignments',
    author='wbazant',
    author_email='wojciech.bazant@gmail.com',
    license="MIT",
    entry_points={"console_scripts": [
        "marker_alignments = marker_alignments.main:main",
    ]},
    install_requires=["pysam", "scipy", "numpy", "sklearn"],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)
