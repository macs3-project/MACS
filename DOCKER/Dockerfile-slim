# Use an official Python runtime as a parent image
FROM python:2.7 AS build

# install numpy and MACS2
RUN pip install --trusted-host pypi.python.org numpy
RUN pip install --trusted-host pypi.python.org MACS2

FROM python:2.7-slim

# copy compiled files over to python-slim
COPY --from=build /usr/local/bin/macs2 /usr/local/bin/
COPY --from=build /usr/local/lib/python2.7/site-packages/ /usr/local/lib/python2.7/site-packages/

# Define environment variable
ENTRYPOINT ["macs2"]
