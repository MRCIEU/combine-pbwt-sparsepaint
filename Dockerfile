FROM alpine:3.23 AS build

RUN apk add --update --no-cache \
        build-base \
        cmake \
        curl

RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs/ | sh -s -- -y

COPY . combine/

WORKDIR /combine

RUN source $HOME/.cargo/env && \
    cargo build --release

# FROM alpine:3.23

# RUN apk add --update --no-cache \
#         procps

FROM scratch

COPY --from=build /combine/target/release/combine /usr/local/bin/combine

CMD ["combine", "-h"]
