# Build stage: Alpine (musl), so the binary is statically linked and can run from scratch.
# If the Alpine version changes, check the musl version noted in notice.hbs.
FROM rust:1.94-alpine3.23 AS build

# build-base and cmake are needed to compile the C compression libraries (libdeflate, zlib-ng)
RUN apk add --no-cache \
        build-base \
        cmake

WORKDIR /combine

# NOTICE is embedded in the binary (include_str!), so it is needed at build time
COPY Cargo.toml Cargo.lock LICENSE NOTICE ./
COPY src ./src

RUN cargo build --release --locked

# Final stage: just the binary, plus its licence and third-party notices
FROM scratch

LABEL org.opencontainers.image.title="combine" \
      org.opencontainers.image.description="Combine per-chromosome pbwt -paintSparse output into a single sparse matrix in map units" \
      org.opencontainers.image.source="https://github.com/mrcieu/combine-pbwt-sparsepaint" \
      org.opencontainers.image.licenses="GPL-3.0-only"

COPY --from=build /combine/target/release/combine /usr/local/bin/combine
COPY --from=build /combine/LICENSE /combine/NOTICE /usr/share/licenses/combine/

CMD ["combine", "-h"]

# Alternative final stage with a shell (e.g. for interactive use), instead of scratch:
#
# FROM alpine:3.23
# RUN apk add --no-cache procps
# COPY --from=build /combine/target/release/combine /usr/local/bin/combine
# COPY --from=build /combine/LICENSE /combine/NOTICE /usr/share/licenses/combine/
# CMD ["combine", "-h"]
