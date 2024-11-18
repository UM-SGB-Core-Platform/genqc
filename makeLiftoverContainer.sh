cd liftOver-Docker
docker build -t liftover .
docker save liftover -o liftover.tar
singularity build liftover.sif docker-archive://liftover.tar
docker image rm liftover
