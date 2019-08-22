# Use an official Python runtime as a parent image
FROM python:2.7.15

# install numpy and MACS2
RUN git clone https://github.com/taoliu/MACS.git /MACS && cd /MACS
RUN pip install --trusted-host pypi.python.org -r requirements.txt
RUN python setup_w_cython.py install

# Define environment variable
ENTRYPOINT ["macs2"]
