FROM gcr.io/broad-getzlab-workflows/base_image:v0.0.2

WORKDIR build
# build steps go here
# remember to clear the build directory!

WORKDIR /app
ENV PATH=$PATH:/app
COPY <your script> .
