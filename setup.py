from skbuild import setup


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()


setup(
    name="TMscore_py",
    version="0.0.1",
    author="Gesine Cauer",
    description="Python wrapper for original TMscore C++ code",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gesinecauer/TMscore_py",
    project_urls={
        "Bug Tracker": "https://github.com/gesinecauer/TMscore_py/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: C++",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=["TMscore_py"],
    python_requires=">=3.6",
)
