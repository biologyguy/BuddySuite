from distutils.core import setup
setup(
    name = 'BuddySuite',
    packages = ['BuddySuite'], # this must be the same as the name above
    version = '1.alpha',
    description = 'A python toolkit for working with biological data',
    author = 'Stephen Bond',
    author_email = 'biologyguy@gmail.com',
    url = 'https://github.com/biologyguy/BuddySuite',
    download_url = 'https://github.com/biologyguy/BuddySuite/tarball/1.alpha', # git tag 1.alpha -m "Adds a tag so that we can put this on PyPI."
    keywords = ['testing', 'logging', 'example'], # arbitrary keywords
    classifiers = [],
)