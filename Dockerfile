# Start with a base Ubuntu image and install Python
FROM ubuntu:22.04

# Avoid prompts from apt
ARG DEBIAN_FRONTEND=noninteractive

# Install Python and other dependencies, and apply updates
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y python3-pip python3-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*


# Set the working directory in the container
WORKDIR /usr/src/app

# Install SigProfilerMatrixGenerator using pip
RUN pip3 install SigProfilerPlotting==1.3.21

# Create a non-root user named 'spm_user'
RUN useradd -m -s /bin/bash spm_user

# Change the ownership of the /usr/src/app directory and its contents to the new non-root user
RUN chown -R spm_user:spm_user /usr/src/app

# Switch to the non-root user for subsequent commands and when running the container
USER spm_user