# Use an official Python runtime as a parent image
FROM python:2.7.15

# install numpy and MACS2
RUN git clone https://github.com/taoliu/MACS.git /MACS
RUN pip install --trusted-host pypi.python.org --upgrade pip && pip install --trusted-host pypi.python.org -r /MACS/requirements.txt
RUN cd /MACS && python setup_w_cython.py install

# Define environment variable
ENTRYPOINT ["macs2"]
