import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()

with open('LICENSE') as f:
	license=f.read()

setuptools.setup(name='marker_alignments',
    version='0.1',
    description="Process and summarise alignments of metagenomic sequencing reads to reference databases of marker genes",
    long_description_content_type="text/markdown",
    url='http://github.com/wbazant/marker_alignments',
    author='wbazant',
    author_email='wojciech.bazant@gmail.com',
    license=license,
    entry_points={"console_scripts": ["process_marker_alignments = marker_alignments.main:main"]},
    packages=setuptools.find_packages(),
)
