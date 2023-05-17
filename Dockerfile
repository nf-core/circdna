FROM python:3.7

# Set up the working directory
WORKDIR ./

# Copy the requirements file
COPY requirements.txt .

COPY bin/ .

# Install Python packages
RUN pip install --no-cache-dir -r requirements.txt
