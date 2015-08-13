from distutils.core import setup
setup(
    name = 'SeqBuddy',
    packages = ['SeqBuddy'], # this must be the same as the name above
    version = 'v1.0',
    description = 'A python toolkit for working with biological data',
    author = 'Stephen Bond',
    author_email = 'biologyguy@gmail.com',
    url = 'https://github.com/biologyguy/BuddySuite',
    download_url = 'https://github.com/biologyguy/BuddySuite/archive/SeqBuddy_v1.0.tar.gz', # git tag 1.alpha -m "Adds a tag so that we can put this on PyPI."
    keywords = ['testing', 'logging', 'example'], # arbitrary keywords
    classifiers = [],
    install_requires=['biopython',
                      ],
)