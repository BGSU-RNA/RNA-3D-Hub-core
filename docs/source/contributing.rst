Contributing
============

Overview
--------

In general we try to follow the same guidelines that numpy does. You should
read `them https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt`_.

Here is a quick summary:

1. We try to follow `pep8 http://pep8.org/`_.

2. Always use a linter. Below is a useful list:

   - `pylint <http://www.logilab.org/857>`_
   - `pyflakes <https://pypi.python.org/pypi/pyflakes>`_
   - `pep8.py <http://svn.browsershots.org/trunk/devtools/pep8/pep8.py>`_
   - `flake8 <https://pypi.python.org/pypi/flake8>`_

   I've used flake8 a lot. Many editors have plugins to make using these tools
   easier. For example, `vim-flake8 <https://github.com/nvie/vim-flake8>`_. Please
   use one.

3. Write tests. We use `py.test <>`.

Testing
-------

Testing the pipeline is essential to developing it. However, it can be fairly
tricky
