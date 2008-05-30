Orange Extensions
V1.8 - Apr 2003

This library is copyleft 2001-2003 by Aleks Jakulin (http://ai.fri.uni-lj.si/~aleks).

It requires Orange 0.7 (http://magix.fri.uni-lj.si/orange) and Python 2.2.

The new version of Orange broke the existing MarginMetaLearner. Surprisingly
nobody complained. This is fixed, and the library is now ported to the new
version of SWIG.

With the help of Daniel Rubin, a Solaris version (which should also compile on
Linux) has been prepared.

BUGS:
- memory leaks in logistic regression and SVM
- no handling for out-of-memory situations (crash)
