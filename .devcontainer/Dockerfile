FROM mcr.microsoft.com/devcontainers/base:ubuntu
RUN apt-get update && \
    apt-get install -y \
    xorg-dev mesa-utils xvfb libgl1 freeglut3-dev libxrandr-dev libxinerama-dev libxcursor-dev \
    libxi-dev libxext-dev libcairo2-dev libfreetype6-dev libffi-dev libjpeg-dev libpng-dev libz-dev \
    libssl-dev curl
USER vscode
SHELL [ "/bin/bash", "-c" ]
ENV SHELL=/bin/bash
RUN curl -fsSL https://install.julialang.org | sh -s -- --yes && \
    ~/.juliaup/bin/juliaup add 1.10.4 && \
    ~/.juliaup/bin/juliaup default 1.10.4