# This script contains commands to create docker images for multi-arch
# testing of MACS3. It will take the 'base' images created by
# make_docker_base.sh, then copy the source code at current directory
# and install MACS3 into the container to make a new image.

# base images are from make_docker_base.sh, please run it first

# create Docker image for testing MACS3 in ppc64le
docker run --rm --privileged multiarch/qemu-user-static --reset -p yes &&
    docker build --rm --build-arg BASE_IMAGE="macs_ppc64le_base" -t macs_ppc64le -f test/Dockerfile.multiarch.install.macs3 .

# create Docker image for testing MACS3 in i386
docker run --rm --privileged multiarch/qemu-user-static --reset -p yes &&
    docker build --rm --build-arg BASE_IMAGE="macs_i386_base" -t macs_i386 -f test/Dockerfile.multiarch.install.macs3 .

# create Docker image for testing MACS3 in arm32v7pp
docker run --rm --privileged multiarch/qemu-user-static --reset -p yes &&
    docker build --rm --build-arg BASE_IMAGE="macs_arm32v7_base" -t macs_arm32v7 -f test/Dockerfile.multiarch.install.macs3 .

# create Docker image for testing MACS3 in arm64v8
docker run --rm --privileged multiarch/qemu-user-static --reset -p yes &&
    docker build --rm --build-arg BASE_IMAGE="macs_arm64v8_base" -t macs_arm64v8 -f test/Dockerfile.multiarch.install.macs3 .

# create Docker image for testing MACS3 in s390x
docker run --rm --privileged multiarch/qemu-user-static --reset -p yes &&
    docker build --rm --build-arg BASE_IMAGE="macs_s390x_base" -t macs_s390x -f test/Dockerfile.multiarch.install.macs3 .

# The final images are macs_ppc64le, macs_i386, macs_arm32v7, macs_arm64v8, and macs_s390x
