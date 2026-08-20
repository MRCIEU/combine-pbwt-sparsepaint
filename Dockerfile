FROM alpine:3.23

RUN apk add --update --no-cache \
        curl \
        gcc \
        musl-dev \
        procps

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs/ | sh -s -- -y

COPY . combine/

WORKDIR /combine

RUN source $HOME/.cargo/env && \
    cargo build --release && \
    mv target/release/combine /usr/local/bin/combine && \
    rm -rf /combine

CMD ["combine", "-h"]
