Generates a Docker image and pushes it to a container registry.

Make sure you are logged into the container registry first:
```
# for Docker Hub, user ogdf:
docker login -u ogdf

# or, for a different container registry:
docker login MY_CONTAINER_REGISTRY
```

Then invoke the build script, using the registry and the tag of the base image
as arguments. Examples:
```
# for Docker Hub, user ogdf:
./build.sh ogdf gcc:10

# or, for a different container registry:
./build.sh MY_CONTAINER_REGISTRY/ogdf gcc:10
```

This will fetch the gcc:10 base image from Docker Hub, install the libraries
necessary for the OGDF on top of it, and publish the new image.

If the base image is of the form 'clang:{version number}', first build a clang
image and then the actual CI-image on top.
