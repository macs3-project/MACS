# Use an official Python runtime as a parent image
FROM python:3.7

# install numpy and MACS2
RUN git clone https://github.com/taoliu/MACS.git /MACS
RUN pip install --trusted-host pypi.python.org --upgrade pip && pip install --trusted-host pypi.python.org -r /MACS/requirements.txt
RUN cd /MACS && python setup.py install

# Define environment variable
ENTRYPOINT ["macs2"]
