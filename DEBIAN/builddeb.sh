#!/bin/bash

cd ../
#rm -rf web
#rm -f TODO
rm -rf debian
mkdir -p debian/usr/
mkdir -p debian/usr/share/doc/macs
mkdir -p debian/DEBIAN/
cp DEBIAN/control debian/DEBIAN/
cp DEBIAN/copyright debian/usr/share/doc/macs
mkdir -p debian/usr/share/man/man1/
makeinfo DEBIAN/macs.texi
info2man DEBIAN/macs.info
cp DEBIAN/macs.man debian/usr/share/man/man1/macs2.1
gzip --best debian/usr/share/man/man1/macs2.1
gzip --best -c DEBIAN/changelog > debian/usr/share/doc/macs/changelog.gz
echo "MACS Debian maintainer and upstream author are identical.\nTherefore see also normal changelog file for Debian changes." > debian/usr/share/doc/macs/changelog.Debian
gzip --best debian/usr/share/doc/macs/changelog.Debian
python2.6 setup.py install --prefix=debian/usr/ --install-layout=deb
#for f in `ls debian/usr/bin/*.py`;do mv ${f} ${f/.py/};done
for f in `find debian/usr/lib/ -name '*.pyc'`;do rm -f ${f};done
fakeroot dpkg-deb --build debian
mv debian.deb macs_2.0.deb
echo "finished!"