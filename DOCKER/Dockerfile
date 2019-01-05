# Use an official Python runtime as a parent image
FROM python:2.7

# Set the working directory to /app
WORKDIR /

RUN pip install --trusted-host pypi.python.org numpy
RUN pip install --trusted-host pypi.python.org MACS2

# Define environment variable
ENV NAME MACS2

CMD ["macs2", "-h"]
