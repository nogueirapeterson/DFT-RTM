FROM docker://nvcr.io/nvidia/nvhpc:22.5-devel-cuda_multi-ubuntu20.04

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip

RUN pip install meson ninja
