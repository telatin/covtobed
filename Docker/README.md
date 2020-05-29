# Docker container

This directories contains files to generate the Docker/Singularity images that provide:

* `covtobed` as described in the [README](../README.md)
* `coverage` a legacy program with some more features that has been used for the EU project _MD-PAEDIGREE_.
 
## How to get it
 
The Docker image is distributed via [Docker Hub](https://hub.docker.com/r/andreatelatin/covtobed):
 
```bash
sudo docker pull andreatelatin/covtobed
```
 
The Singularity image can be pulled from Docker Hub as well. A snapshot (v. 0.3) is available from [Zenodo](https://zenodo.org/record/1063493):

```bash
singularity pull andreatelatin/covtobed
```
