# Use an official Python runtime as a parent image
FROM python:2.7.15

# install numpy and MACS2
RUN pip install --trusted-host pypi.python.org numpy
RUN pip install --trusted-host pypi.python.org MACS2

# Define environment variable
ENTRYPOINT ["macs2"]
