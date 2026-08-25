# ==========================================
# Stage 1: Build the SHEsis C++ Core Binary
# ==========================================
FROM ubuntu:16.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    python3 \
    libboost-all-dev \
    libarmadillo-dev \
    libmlpack-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build
COPY . .

# Run external patch script
RUN python3 patch.py

# Compile SHEsis binary
RUN g++ -O3 -std=c++11 -include iostream \
    -D__STDC_LIMIT_MACROS -D__STDC_FORMAT_MACROS \
    -Wno-parentheses -Wno-write-strings -Wno-conversion-null \
    -I/usr/include/libxml2 \
    main.cpp HaplotypeLD.cpp SHEsisData.cpp fisher.cpp utility.cpp \
    AssociationTest.cpp HWETest.cpp LDTest.cpp QTL.cpp HaplotypeBase.cpp \
    Haplotype.cpp HaplotypeEM.cpp IndexingVariables.cpp ArrayStorage.cpp \
    System.cc Solver.cc Options.cc BMP.cpp font.cpp minifont.cpp \
    CreatHtmlTable.cpp linear.cpp regression.cpp logistic.cpp \
    MarkerRegression.cpp GeneInteractionQTL.cpp GeneInteractionBinary.cpp GeneInteraction.cpp \
    -o SHEsis -lmlpack -larmadillo -lboost_program_options -lxml2

# ==========================================
# Stage 2: Runtime Environment (Web UI & API)
# ==========================================
FROM ubuntu:16.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    ca-certificates \
    xz-utils \
    libarmadillo6 \
    libmlpack2 \
    libboost-program-options1.58.0 \
    libxml2 \
    && rm -rf /var/lib/apt/lists/*

# Install Node.js 10 runtime
RUN curl -fsSL https://nodejs.org/dist/v10.24.1/node-v10.24.1-linux-x64.tar.xz | tar -xJ -C /usr/local --strip-components=1

WORKDIR /app/SHEsisWebServer

# Install Node dependencies
COPY SHEsisWebServer/package.json .
RUN npm install && npm install shelljs@^0.8.0 mongodb@^2.2.36 --save

# Copy patched web server files from builder stage
COPY --from=builder /build/SHEsisWebServer ./

# Copy compiled SHEsis binary into bin/
COPY --from=builder /build/SHEsis ./bin/SHEsis
RUN chmod +x ./bin/SHEsis

# Ensure working directory permissions
RUN mkdir -p public/tmp && chmod -R 777 public/tmp

EXPOSE 5903

CMD ["node", "SHEsisServer.js"]
