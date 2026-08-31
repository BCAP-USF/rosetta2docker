# rosetta2docker

Using Rosetta in a container, managed with Docker Compose. Your data lives in the `./data` directory and is mounted into each container at `/data`.

## Prerequisites

- [Docker](https://docs.docker.com/get-docker/) and [Docker Compose](https://docs.docker.com/compose/install/)
- The `rosettacommons/rosetta` images are pulled from [Docker Hub](https://hub.docker.com/r/rosettacommons/rosetta) the first time you start a container.

## Getting started

The compose file defines four services. Start them with:

```bash
docker compose up -d
```

Start a single service with:

```bash
docker compose up -d <service>
```

Place your input files in the `./data` directory. They appear in the container at `/data`.

## Data sharing

All services mount the local `./data` folder at `/data` inside the container, so anything you write there (results, outputs) is available on your host immediately.

## Services

### rosetta (serial)

Single-process Rosetta commands.

```bash
docker exec -it rosetta bash
```

### rosetta-mpi (MPI)

MPI-enabled build, for running parallel jobs with `mpirun`.

```bash
docker exec -it rosetta-mpi bash
```

### rosetta-jupyter (Jupyter)

Jupyter notebook server. Access it at <http://localhost:8888> (port `8888` is published).

```bash
docker compose up -d rosetta-jupyter
```

### rosetta-ml (Machine Learning)

Rosetta build with Python machine-learning dependencies (e.g. PyRosetta, TensorFlow).

```bash
docker exec -it rosetta-ml bash
```

## Stopping the containers

```bash
docker compose down
```

To also remove the containers after they stop:

```bash
docker compose stop
```
