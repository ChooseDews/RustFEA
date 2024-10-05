FROM rust:bullseye

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    liblapacke-dev \
    libopenblas-dev \
    libsuitesparse-dev \
    libmumps-seq-dev \
    bash    # Ensure bash is installed

WORKDIR /app

# Avoid copying the src directory, since we will mount it
COPY ./examples ./examples
COPY ./Cargo.toml ./Cargo.toml
COPY ./Cargo.lock ./Cargo.lock
#COPY ./src ./src

RUN cargo build --release
