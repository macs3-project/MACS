docker run --rm -t macs_ppc64le bash -cx "pytest --runxfail && cd test && ./cmdlinetest macs2" &> ../temp/ppc64le.txt
docker run --rm -t macs_i386 bash -cx "pytest --runxfail && cd test && ./cmdlinetest macs2" &> ../temp/i386.txt
docker run --rm -t macs_arm32v7 bash -cx "pytest --runxfail && cd test && ./cmdlinetest macs2" &> ../temp/arm32v7.txt
docker run --rm -t macs_arm64v8 bash -cx "pytest --runxfail && cd test && ./cmdlinetest macs2" &> ../temp/arm64v8.txt
docker run --rm -t macs_s390x bash -cx "pytest --runxfail && cd test && ./cmdlinetest macs2" &> ../temp/s390x.txt
